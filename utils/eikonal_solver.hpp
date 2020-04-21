#include <math.h>
#include <iostream>

#define EK_EPS 1e-16
// #define EK_HIGHPREC

#ifdef EK_HIGHPREC
  #define EK_ACC 1e-4
#else
  #define EK_ACC 1e-4
#endif

#define EK_COEFSZ 16 // the size of the coefficients block, at least 14 is needed
#define EK_INF 1e100

#ifdef EK_HIGHPREC
#define ek_double long double
#define ek_sqrt sqrtl
#else
#define ek_double double
#define ek_sqrt sqrt
#endif


//#define EK_PRINT_STATS


template<class T>
struct ek_coeff
{
  // tet
  T a, b, c, d, e, f;

  // tri
  T a1, a2, a3, b1, b2, b3, c1, c2, c3;

  // line
  T l1, l2, l3;
};

template<class S>
struct fib_velo {
  S vf;
  S vs;
  S vn;
};


template<class T, class S>
struct ek_initdata
{
  mt_vector<T> nod;  ///< Activation nodes.
  mt_vector<S> val;  ///< Activation times.

  // default global velocities
  S vf;   ///< fiber velo
  S vs;   ///< sheet velo
  S vn;   ///< normal velo
  S vps;  ///< PS velo
  /// tag specific velocities
  std::map<T, fib_velo<S> > tag_velo;

  // PS related delays
  S antero_delay;
  S retro_delay;
};

// a, b vectors of dim 3
// M    3x3 matrix
// compute a^T * M * b
template<class T> inline
T comp_3x3_aMb_symm(const T* a, const T* M, const T* b)
{
  T ax = a[0], ay = a[1], az = a[2];
  T bx = b[0], by = b[1], bz = b[2];

  T M0 = M[0], M1 = M[1], M2 = M[2],
               M4 = M[3], M5 = M[4],
                          M8 = M[5];

  return (M0*ax + M1*ay + M2*az)*bx + (M1*ax + M4*ay + M5*az)*by + (M2*ax + M5*ay + M8*az)*bz;
}

// Assuming B \subsetneq A
// compute C := A \ B
template<class T>
void set_cut(const T* A, int asize, const T* B, int bsize, T* C)
{
  int widx = 0;

  for(int i=0; i<asize; i++){
    bool inside = false;
    T aval = A[i];
    for(int j=0; j<bsize; j++) {
      if(aval == B[j]) {
        inside = true;
        break;
      }
    }
    if(! inside) {
      C[widx] = aval;
      widx++;
    }
  }
}

template<class T>
T get_tet_base_and_vpos(const T* tet, const T & v, T* tet_base)
{
  int widx = 0;
  const int asize = 4;
  T vpos = -1;

  for(int i=0; i<asize; i++){
    if(tet[i] != v) {
      tet_base[widx] = tet[i];
      widx++;
    }
    else {
      vpos = i;
    }
  }
  return vpos;
}


template<class T, class S, class V>
class eikonal_solver
{
private:
  // e2n: graph connecting elements to nodes
  // n2e: graph connecting nodes to elements, thus e2n transposed
  // n2n: graph connecting nodes with nodes
  const mt_meshdata & _mesh;

  // initial nodes and their respective value
  const ek_initdata<T, V> & _init;

  size_t _numNodes;
  size_t _numElems;

  // flag specifiying if 1 or 2 vectors define the element fiber orientation
  bool _2lon;

  V _vf;
  V _vs;
  V _vn;
  V _vps;

  V _antero_inc;
  V _retro_inc;

  V _initval;

  mt_vector<T> _act_nod;
  size_t            _act_idx;

  mt_vector<bool> _fin_nod_map;
  mt_vector<bool> _act_nod_map;
  mt_vector<bool> _jun_nod_map;

  mt_vector<T>         _nod_upd;
  mt_vector<ek_double> _sca_coeff;


  inline void nodal_expand(mt_vector<V> & phi, const T v, const V vval)
  {
    bool active;
    mt_int start = _mesh.n2n_dsp[v], stop = start + _mesh.n2n_cnt[v];

    for(int j = start; j < stop; j++) {
      mt_int nidx = _mesh.n2n_con[j];

      if( (nidx != v) && (_fin_nod_map[nidx] == false) ) {
        #ifdef OPENMP
        #pragma omp atomic read
        #endif
        active = _act_nod_map[nidx];

        if( !active ) {
          V nval = phi[nidx];

          if(nval > vval) {
            V newval = solve_nbh(phi, nidx, true);
            if( (nval - newval) > (EK_ACC * (1 + newval)) )
            {
              #ifdef OPENMP
              #pragma omp atomic write
              #endif
              phi[nidx] = newval;

              _nod_upd[nidx]++;

              #ifdef OPENMP
              #pragma omp atomic write
              #endif
              _act_nod_map[nidx] = true;
            }
          }
        }
      }
    }
  }



  // compute the velocity information M based on the
  // two fiber directions in _lon
  inline void comp_M_2lon(ek_double* M, const size_t eidx)
  {
    const S* f = _mesh.lon.data() + eidx*6;
    const S* s = f+3;
    S f0 = f[0], f1 = f[1], f2 = f[2], s0 = s[0], s1 = s[1], s2 = s[2];

    // if current element tag is in the map fib_velo, take the respective velocities
    // else take default velocities
    S vf, vs, vn;
    T tag = _mesh.etags[eidx];
    if(_init.tag_velo.count(tag)) {
      const struct fib_velo<V> & velo = _init.tag_velo.find(tag)->second;
      vf = velo.vf; vf = 1. / vf / vf;
      vs = velo.vs; vs = 1. / vs / vs;
      vn = velo.vn; vn = 1. / vn / vn;
    }
    else {
      vf = _vf;
      vs = _vs;
      vn = _vn;
    }

    #if 0
    M[0] = f0*f0*_vf + s0*s0*_vs - (f0*f0 + s0*s0 - 1)*_vn;
    M[1] = f0*f1*_vf + s0*s1*_vs - (f0*f1 + s0*s1)    *_vn;
    M[2] = f0*f2*_vf + s0*s2*_vs - (f0*f2 + s0*s2)    *_vn;
    //M[3] = f0*f1*_vf + s0*s1*_vs - (f0*f1 + s0*s1)    *_vn;
    M[4] = f1*f1*_vf + s1*s1*_vs - (f1*f1 + s1*s1 - 1)*_vn;
    M[5] = f1*f2*_vf + s1*s2*_vs - (f1*f2 + s1*s2)    *_vn;
    //M[6] = f0*f2*_vf + s0*s2*_vs - (f0*f2 + s0*s2)    *_vn;
    //M[7] = f1*f2*_vf + s1*s2*_vs - (f1*f2 + s1*s2)    *_vn;
    M[8] = f2*f2*_vf + s2*s2*_vs - (f2*f2 + s2*s2 - 1)*_vn;
    #endif

    M[0] = f0*f0*vf + s0*s0*vs - (f0*f0 + s0*s0 - 1)*vn;
    M[1] = f0*f1*vf + s0*s1*vs - (f0*f1 + s0*s1)    *vn;
    M[2] = f0*f2*vf + s0*s2*vs - (f0*f2 + s0*s2)    *vn;
    M[3] = f1*f1*vf + s1*s1*vs - (f1*f1 + s1*s1 - 1)*vn; // M4
    M[4] = f1*f2*vf + s1*s2*vs - (f1*f2 + s1*s2)    *vn; // M5
    M[5] = f2*f2*vf + s2*s2*vs - (f2*f2 + s2*s2 - 1)*vn; // M8
  }


  // compute the velocity information M based on the
  // one fiber direction in _lon
  inline void comp_M_1lon(ek_double* M, const size_t eidx)
  {
    const S* f = _mesh.lon.data() + eidx*3;
    S f0 = f[0], f1 = f[1], f2 = f[2];

    // if current element tag is in the map fib_velo, take the respective velocities
    // else take default velocities
    S vf, vs;
    T tag = _mesh.etags[eidx];
    if(_init.tag_velo.count(tag)) {
      const struct fib_velo<V> & velo = _init.tag_velo.find(tag)->second;
      vf = velo.vf; vf = 1. / vf / vf;
      vs = velo.vs; vs = 1. / vs / vs;
    }
    else {
      vf = _vf;
      vs = _vs;
    }

    #if 0
    M[0] = f0*f0*_vf - (f0*f0 - 1)*_vs;
    M[1] = f0*f1*_vf - f0*f1      *_vs;
    M[2] = f0*f2*_vf - f0*f2      *_vs;
    M[3] = f0*f1*_vf - f0*f1      *_vs; // M[3] = M[1];
    M[4] = f1*f1*_vf - (f1*f1 - 1)*_vs;
    M[5] = f1*f2*_vf - f1*f2      *_vs;
    M[6] = f0*f2*_vf - f0*f2      *_vs; // M[6] = M[2];
    M[7] = f1*f2*_vf - f1*f2      *_vs; // M[7] = M[5];
    M[8] = f2*f2*_vf - (f2*f2 - 1)*_vs;
    #endif

    M[0] = f0*f0*vf - (f0*f0 - 1)*vs;
    M[1] = f0*f1*vf - f0*f1      *vs;
    M[2] = f0*f2*vf - f0*f2      *vs;
    M[3] = f1*f1*vf - (f1*f1 - 1)*vs; // M4
    M[4] = f1*f2*vf - f1*f2      *vs; // M5
    M[5] = f2*f2*vf - (f2*f2 - 1)*vs; // M8
  }

  // compute velocity information for a purkinje fiber
  inline void comp_M_PS(ek_double* M, const size_t eidx, const T a, const T b)
  {
    const S* xyz = _mesh.xyz.data(), *xp;
    xp = xyz + a*3;
    S ax = *xp++, ay = *xp++, az = *xp++;
    xp = xyz + b*3;
    S bx = *xp++, by = *xp++, bz = *xp++;
    S ex = bx - ax, ey = by - ay, ez = bz - az;

    S l = sqrt(ex*ex + ey*ey + ez*ez);
    ex /= l, ey /= l, ez /= l;

    V vps;
    T tag = _mesh.etags[eidx];
    if(_init.tag_velo.count(tag)) {
      const struct fib_velo<V> & velo = _init.tag_velo.find(tag)->second;
      vps = velo.vf * 1e3; vps = 1. / vps / vps;
    }
    else
      vps = _vps;

    M[0] = vps * ex*ex;
    M[1] = vps * ex*ey;
    M[2] = vps * ex*ez;
    M[3] = vps * ey*ey;
    M[4] = vps * ey*ez;
    M[5] = vps * ez*ez;
  }



  // solve tetrahedra for vertex v4
  inline V solve_tetra(const V phi1, const V phi2, const V phi3,
                const ek_double a, const ek_double b, const ek_double c,
                const ek_double d, const ek_double e, const ek_double f)
  {
    // phi
    const V phi23 = phi3 - phi2;
    const V phi13 = phi3 - phi1;

    bool phi23_bad = fabs(phi23) == EK_INF;
    bool phi13_bad = fabs(phi13) == EK_INF;

    if(phi23_bad || phi13_bad)
      return V(EK_INF);

    // solve quadratic system

    ek_double sqrt_val = -e*e*c*d*d + a*d*d*d*d - (b*b*c - a*c*c)*f*f + ((b*b - a*c)*f*f -
                          (e*e*c - a*d*d - 2*(e*c - b*d)*e)*f)*phi13*phi13 - 2*(e*e*c*d -
                          2*e*b*d*d + a*d*d*d + (b*b - a*c)*d*f)*phi13*phi23 + (e*e*c*c -
                          2*e*b*c*d + a*c*d*d + (b*b*c - a*c*c)*f)*phi23*phi23 + 2*(e*c*d*d -
                          b*d*d*d)*e + (e*e*c*c + (b*b - 2*a*c)*d*d - 2*(e*c*c - b*c*d)*e)*f;

    if(sqrt_val < 0)
      return V(EK_INF);

    sqrt_val = ek_sqrt(sqrt_val);

    const ek_double p1 = e*c*d*d - b*d*d*d + (e*c - b*d)*f*phi13*phi13 - 2*(e*c*d - b*d*d)*phi13*phi23 +
                         (e*c*c - b*c*d)*phi23*phi23 - (e*c*c - b*c*d)*f;

    const ek_double p2 = (d*d*d*d - 2*c*d*d*f + c*c*f*f + (d*d*f - c*f*f)*phi13*phi13 - 2*(d*d*d - c*d*f)*phi13*phi23
                         + (c*d*d - c*c*f)*phi23*phi23);

    const ek_double l2_1 = (p1 + sqrt_val*(d*phi13 - c*phi23)) / p2;
    const ek_double l2_2 = (p1 - sqrt_val*(d*phi13 - c*phi23)) / p2;

    const ek_double l1_1 = -((f*l2_1 + e)*phi13 - (d*l2_1 + b)*phi23)/(d*phi13 - c*phi23);
    const ek_double l1_2 = -((f*l2_2 + e)*phi13 - (d*l2_2 + b)*phi23)/(d*phi13 - c*phi23);

    // compute phi4
    V phi4_1 = EK_INF, phi4_2 = EK_INF;

    V l3 = 1 - l1_1 - l2_1;
    if( (l1_1 > 0. && l1_1 < 1.) && (l2_1 > 0. && l2_1 < 1.) && l3 > 0 )
    {
      phi4_1 = l1_1 * phi1 + l2_1 * phi2 + l3*phi3 +
               ek_sqrt(a + 2*l1_1*b + l1_1*l1_1*c + 2*l1_1*l2_1*d + 2*l2_1*e + l2_1*l2_1*f);
    }
    l3 = 1 - l1_2 - l2_2;
    if( (l1_2 > 0. && l1_2 < 1.) && (l2_2 > 0. && l2_2 < 1.) && l3 > 0 )
    {
      phi4_2 = l1_2 * phi1 + l2_2 * phi2 + l3*phi3 +
               ek_sqrt(a + 2*l1_2*b + l1_2*l1_2*c + 2*l1_2*l2_2*d + 2*l2_2*e + l2_2*l2_2*f);
    }

    V minsol = phi4_1;
    if(minsol > phi4_2) minsol = phi4_2;

    return minsol;
  }


  // solve line for v2
  inline V solve_line(const ek_double* M,
                      const mt_vector<V> & phi,
                      const T v1, const T v2)
  {
     // we solve travel time from x1 to x2
    ek_double e12[3];
    const S* xyz = _mesh.xyz.data(), *xp;
    xp = xyz + v1*3;
    S x1_x = *xp++, x1_y = *xp++, x1_z = *xp++;
    xp = xyz + v2*3;
    S x2_x = *xp++, x2_y = *xp++, x2_z = *xp++;

    e12[0] = x2_x - x1_x; e12[1] = x2_y - x1_y; e12[2] = x2_z - x1_z;

    V phi1;
    #ifdef OPENMP
    #pragma omp atomic read
    #endif
    phi1 = phi[v1];

    V sol = phi1 + sqrt(comp_3x3_aMb_symm(e12, M, e12));
    return sol;
  }



  // The solve_triangle version using precomputed coeffs
  inline V solve_triangle(const V f1,
                   const V f2,
                   const ek_double a,
                   const ek_double b,
                   const ek_double c)
  {
    // phi
    V f12 = f2 - f1;

    //bool f12_bad = (fabs(f12) < EK_EPS) || (fabs(f12) == EK_INF);
    bool f12_bad = fabs(f12) == EK_INF;
    if(f12_bad)
      return EK_INF;

    // compute edges
    ek_double sqrt_val, p1, p2, l1, l2;
    sqrt_val = -b*b*c + a*c*c + (b*b - a*c)*f12*f12;
    if(sqrt_val < 0)
      return EK_INF;

    sqrt_val = ek_sqrt(sqrt_val);
    p1 = b*f12*f12 - b*c;
    p2 = c*f12*f12 - c*c;

    l1 = -( p1 - sqrt_val*f12 ) / p2;
    l2 = -( p1 + sqrt_val*f12 ) / p2;

    V ret1 = EK_INF, ret2 = EK_INF;

    if( (l1 > 0) && (l1 < 1) )
      ret1 = l1 * f1 + (1-l1) * f2 + ek_sqrt(a + 2*b*l1 + l1*l1*c);

    if( (l2 > 0) && (l2 < 1) )
      ret2 = l2 * f1 + (1-l2) * f2 + ek_sqrt(a + 2*b*l2 + l2*l2*c);

    V minsol = ret1;
    if(minsol > ret2) minsol = ret2;

    return minsol;
  }
  // The solve_triangle version computing the products on-the-fly
  inline V solve_triangle(const ek_double* M,
                   const mt_vector<V> & phi,
                   const T v1, const T v2, const T v3)
  {
    // phi
    V f1, f2;

    #ifdef OPENMP
    #pragma omp atomic read
    #endif
    f1 = phi[v1];
    #ifdef OPENMP
    #pragma omp atomic read
    #endif
    f2 = phi[v2];

    V f12 = f2 - f1;

    bool f12_bad = (fabs(f12) < EK_EPS) || (fabs(f12) == EK_INF);
    if(f12_bad)
      return EK_INF;

    // compute edges
    const S* xyz = _mesh.xyz.data(), *xp;
    ek_double e12[3], e23[3];

    xp = xyz + v1*3;
    S x1_x = *xp++, x1_y = *xp++, x1_z = *xp++;
    xp = xyz + v2*3;
    S x2_x = *xp++, x2_y = *xp++, x2_z = *xp++;
    xp = xyz + v3*3;
    S x3_x = *xp++, x3_y = *xp++, x3_z = *xp++;

    e12[0] = x2_x - x1_x; e12[1] = x2_y - x1_y; e12[2] = x2_z - x1_z;
    e23[0] = x3_x - x2_x; e23[1] = x3_y - x2_y; e23[2] = x3_z - x2_z;

    long double a = comp_3x3_aMb_symm(e23, M, e23);
    long double b = comp_3x3_aMb_symm(e23, M, e12);
    long double c = comp_3x3_aMb_symm(e12, M, e12);

    ek_double sqrt_val, p1, p2, l1, l2;
    sqrt_val = -b*b*c + a*c*c + (b*b - a*c)*f12*f12;
    if(sqrt_val < 0)
      return EK_INF;

    sqrt_val = ek_sqrt(sqrt_val);
    p1 = b*f12*f12 - b*c;
    p2 = c*f12*f12 - c*c;

    l1 = -( p1 - sqrt_val*f12 ) / p2;
    l2 = -( p1 + sqrt_val*f12 ) / p2;

    V ret1 = EK_INF, ret2 = EK_INF;

    if( (l1 > 0) && (l1 < 1) )
      ret1 = l1 * f1 + (1-l1) * f2 + ek_sqrt(a + 2*b*l1 + l1*l1*c);

    if( (l2 > 0) && (l2 < 1) )
      ret2 = l2 * f1 + (1-l2) * f2 + ek_sqrt(a + 2*b*l2 + l2*l2*c);

    V minsol = ret1;
    if(minsol > ret2) minsol = ret2;

    return minsol;
  }


  // get 1D minimal phi4 at v4 for one tetra defined by v1 to v4
  inline V get_min_1D(const V phi1, const V phi2, const V phi3, struct ek_coeff<ek_double> & c)
  {
    V sol1, sol2, sol3, minsol;

    // line (1,4)
    sol1 = phi1 + sqrt(c.l1);
    // line (2,4)
    sol2 = phi2 + sqrt(c.l2);
    // line (3,4)
    sol3 = phi3 + sqrt(c.l3);

    minsol = sol1;
    if(minsol > sol2) minsol = sol2;
    if(minsol > sol3) minsol = sol3;

    return minsol;
  }

  // get 2D minimal phi4 at v4 for one tetra defined by v1 to v4
  inline V get_min_2D(const V phi1, const V phi2, const V phi3, struct ek_coeff<ek_double> & c)
  {
    V sol1, sol2, sol3, minsol;

    // triangle {1,2,4}
    sol1 = solve_triangle(phi1, phi2, c.a1, c.b1, c.c1);
    // triangle {2,3,4}
    sol2 = solve_triangle(phi2, phi3, c.a2, c.b2, c.c2);
    // triangle {1,3,4}
    sol3 = solve_triangle(phi1, phi3, c.a3, c.b3, c.c3);

    minsol = sol1;
    if(minsol > sol2) minsol = sol2;
    if(minsol > sol3) minsol = sol3;

    return minsol;
  }





  // coeff[0 ] = comp_3x3_aMb_symm(e12, M, e12);
  // coeff[1 ] = comp_3x3_aMb_symm(e12, M, e23);
  // coeff[2 ] = comp_3x3_aMb_symm(e12, M, e24);
  // coeff[3 ] = comp_3x3_aMb_symm(e13, M, e13);
  // coeff[4 ] = comp_3x3_aMb_symm(e13, M, e23);
  // coeff[5 ] = comp_3x3_aMb_symm(e13, M, e34);
  // coeff[6 ] = comp_3x3_aMb_symm(e14, M, e14);
  // coeff[7 ] = comp_3x3_aMb_symm(e14, M, e24);
  // coeff[8 ] = comp_3x3_aMb_symm(e14, M, e34);
  // coeff[9 ] = comp_3x3_aMb_symm(e23, M, e23);
  // coeff[10] = comp_3x3_aMb_symm(e23, M, e34);
  // coeff[11] = comp_3x3_aMb_symm(e24, M, e24);
  // coeff[12] = comp_3x3_aMb_symm(e24, M, e34);
  // coeff[13] = comp_3x3_aMb_symm(e34, M, e34);
  inline void get_coeffs(const size_t eidx, const T vpos, struct ek_coeff<ek_double> & c)
  {
    ek_double* coeffs = _sca_coeff.data() + eidx * EK_COEFSZ;

    ek_double C0 = coeffs[0], C1 = coeffs[1], C2 = coeffs[2], C3 = coeffs[3],
              C4 = coeffs[4], C5 = coeffs[5], C6 = coeffs[6], C7 = coeffs[7],
              C8 = coeffs[8], C9 = coeffs[9], C10 = coeffs[10], C11 = coeffs[11],
              C12 = coeffs[12], C13 = coeffs[13];

    switch(vpos)
    {
      case 0:
        // edge -> orig_edge
        // (1,2) -> (2,3), (1,3) ->  (2,4), (1,4) -> -(1,2)
        // (2,3) -> (3,4), (2,4) -> -(1,3), (3,4) -> -(1,4)

        // tet
        //  a =  (e14, M, e14), b = -(e14, M, e24), c = (e24, M, e24)
        //  d =  (e34, M, e24), e = -(e34, M, e14), f = (e34, M, e34)
        c.a = C6;  c.b = -C7; c.c = C11;
        c.d = C12; c.e = -C8; c.f = C13;
        // tri
        // a1 =  (e13, M, e13), a2 =  (e14, M, e14), a3 =  (e14, M, e14)
        // b1 = -(e13, M, e23), b2 = -(e14, M, e34), b3 = -(e14, M, e24)
        // c1 =  (e23, M, e23), c2 =  (e34, M, e34), c3 =  (e24, M, e24)
        c.a1 = C3; c.b1 = -C4; c.c1 = C9;
        c.a2 = C6; c.b2 = -C8; c.c2 = C13;
        c.a3 = C6; c.b3 = -C7; c.c3 = C11;
        // line
        // l1 = (e12, M, e12), l2 = (e13, M, e13), l3 = (e14, M, e14)
        c.l1 = C0; c.l2 = C3; c.l3 = C6;
        break;

      case 1:
        // edge -> orig_edge
        // (1,2) -> (1,3), (1,3) ->  (1,4), (1,4) ->  (1,2)
        // (2,3) -> (3,4), (2,4) -> -(2,3), (3,4) -> -(2,4)

        // tet
        //  a = (e24, M, e24), b = -(e24, M, e14), c = (e14, M, e14)
        //  d = (e34, M, e14), e = -(e34, M, e24), f = (e34, M, e34)
        c.a = C11; c.b = -C7;  c.c = C6;
        c.d = C8;  c.e = -C12; c.f = C13;
        // tri
        // a1 =  (e23, M, e23) a2 =  (e24, M, e24) a3 =  (e24, M, e24)
        // b1 = -(e23, M, e13) b2 = -(e24, M, e34) b3 = -(e24, M, e14)
        // c1 =  (e13, M, e13) c2 =  (e34, M, e34) c3 =  (e14, M, e14)
        c.a1 = C9;  c.b1 = -C4;  c.c1 = C3;
        c.a2 = C11; c.b2 = -C12; c.c2 = C13;
        c.a3 = C11; c.b3 = -C7;  c.c3 = C6;
        // line
        // l1 = (e12, M, e12), l2 = (e23, M, e23), l3 = (e24, M, e24)
        c.l1 = C0; c.l2 = C9; c.l3 = C11;
        break;

      case 2:
        // edge -> orig_edge
        // (1,2) -> (1,2), (1,3) -> (1,4), (1,4) -> (1,3)
        // (2,3) -> (2,4), (2,4) -> (2,3), (3,4) -> -(3,4)

        // tet
        //  a = (e34, M, e34), b = -(e34, M, e14), c = (e14, M, e14)
        //  d = (e24, M, e14), e = -(e24, M, e34), f = (e24, M, e24)
        c.a = C13; c.b = -C8; c.c = C6;
        c.d = C7; c.e = -C12; c.f = C11;
        // tri
        // a1 = (e23, M, e23) a2 =  (e34, M, e34) a3 =  (e34, M, e34)
        // b1 = (e23, M, e12) b2 = -(e34, M, e24) b3 = -(e34, M, e14)
        // c1 = (e12, M, e12) c2 =  (e24, M, e24) c3 =  (e14, M, e14)
        c.a1 = C9;  c.b1 =  C1 ; c.c1 = C0 ;
        c.a2 = C13; c.b2 = -C12; c.c2 = C11;
        c.a3 = C13; c.b3 = -C8;  c.c3 = C6;
        // line
        // l1 = (e13, M, e13), l2 = (e23, M, e23), l3 = (e34, M, e34)
        c.l1 = C3; c.l2 = C9; c.l3 = C13;
        break;

      case 3:
        // This is the original configuration. No mapping needed.

        // tet
        //  a = (e34, M, e34), b = (e34, M, e13), c = (e13, M, e13)
        //  d = (e23, M, e13), e = (e23, M, e34), f = (e23, M, e23)
        c.a = C13; c.b = C5;  c.c = C3 ;
        c.d = C4;  c.e = C10; c.f = C9;
        // tri
        // a1 = (e24, M, e24) a2 = (e34, M, e34) a3 = (e34, M, e34)
        // b1 = (e24, M, e12) b2 = (e34, M, e23) b3 = (e34, M, e13)
        // c1 = (e12, M, e12) c2 = (e23, M, e23) c3 = (e13, M, e13)
        c.a1 = C11; c.b1 = C2 ; c.c1 = C0 ;
        c.a2 = C13; c.b2 = C10; c.c2 = C9;
        c.a3 = C13; c.b3 = C5;  c.c3 = C3 ;
        // line
        // l1 = (e14, M, e14), l2 = (e24, M, e24), l3 = (e34, M, e34)
        c.l1 = C6; c.l2 = C11; c.l3 = C13;
        break;

      default:
        std::cerr << "get_coeffs: ERROR, bad vertex position index!" << std::endl;
        break;
    }
  }


  // solve tetra with index eidx, for vertex with index v
  inline V solve_3D(T* tet_base, const mt_vector<V> & phi, const size_t eidx, const T v, const bool check1D)
  {
    const T* tet = _mesh.e2n_con.data() + _mesh.e2n_dsp[eidx];
    T vpos = get_tet_base_and_vpos(tet, v, tet_base);

    V phi1, phi2, phi3;
    #ifdef OPENMP
    #pragma omp atomic read
    #endif
    phi1 = phi[tet_base[0]];

    #ifdef OPENMP
    #pragma omp atomic read
    #endif
    phi2 = phi[tet_base[1]];

    #ifdef OPENMP
    #pragma omp atomic read
    #endif
    phi3 = phi[tet_base[2]];

    if( (phi1 == EK_INF) && (phi2 == EK_INF) && (phi3 == EK_INF) )
      return EK_INF;

    struct ek_coeff<ek_double> c;
    memset(&c, 0, sizeof(c));
    get_coeffs(eidx, vpos, c);

    if(check1D) {
      #ifdef EK_FAST
      V sol1, sol2, minsol;
      sol1 = solve_tetra(phi1, phi2, phi3, c.a, c.b, c.c, c.d, c.e, c.f);
      if(sol1 == EK_INF) {
        sol1 = get_min_2D(phi1, phi2, phi3, c);
      }
      sol2 = get_min_1D(phi1, phi2, phi3, c);

      if(sol1 < sol2) minsol = sol1;
      else minsol = sol2;
      #else
      V sol1, sol2, sol3, minsol;
      sol1 = solve_tetra(phi1, phi2, phi3, c.a, c.b, c.c, c.d, c.e, c.f);
      sol2 = get_min_2D(phi1, phi2, phi3, c);
      sol3 = get_min_1D(phi1, phi2, phi3, c);

      minsol = sol1;
      if(sol2 < minsol) minsol = sol2;
      if(sol3 < minsol) minsol = sol3;
      #endif
      return minsol;
    }
    else {
      // solve 3D
      V sol = solve_tetra(phi1, phi2, phi3, c.a, c.b, c.c, c.d, c.e, c.f);
      if(sol == EK_INF) {
        // solve 2D
        sol = get_min_2D(phi1, phi2, phi3, c);
        if (sol == EK_INF) {
          // solve 1D
          sol = get_min_1D(phi1, phi2, phi3, c);
        }
      }
      return sol;
    }
  }

  inline V solve_2D(T* base, const mt_vector<V> & phi, const T eidx, const T v)
  {
    const T* tri = _mesh.e2n_con.data() + _mesh.e2n_dsp[eidx];
    set_cut(tri, 3, &v, 1, base);
    T v1 = base[0], v2 = base[1];
    ek_double M[6];

    V phi1, phi2;
    #ifdef OPENMP
    #pragma omp atomic read
    #endif
    phi1 = phi[v1];
    #ifdef OPENMP
    #pragma omp atomic read
    #endif
    phi2 = phi[v2];

    if( (phi1 == EK_INF) && (phi2 == EK_INF) )
      return EK_INF;

    comp_M_1lon(M, eidx);

    // solve 2D
    V sol = solve_triangle(M, phi, v1, v2, v);
    if (sol == EK_INF) {
      // solve 1D
      V sol1, sol2;

      sol1 = solve_line(M, phi, v1, v);
      sol2 = solve_line(M, phi, v2, v);

      sol = sol1;
      if(sol > sol2) sol = sol2;
    }
    return sol;
  }


  // solve line element of Purkinje System
  inline V solve_PS(const mt_vector<V> & phi, const size_t eidx, const T v)
  {
    ek_double M[6];
    const T* line = _mesh.e2n_con.data() + _mesh.e2n_dsp[eidx];

    // determine start and endpoint of line, v1 = start, v2 = end
    T v_start;
    if(line[0] == v) v_start = line[1];
    else v_start = line[0];

    // the line goes from tet_base[0] to v
    comp_M_PS(M, eidx, v_start, v);
    V sol = solve_line(M, phi, v_start, v);
    return sol;
  }

  // solve for vertex v
  // this involves computing 3D-1D solutions for all elements connected to v,
  // then taking the minimal one
  inline V solve_nbh(const mt_vector<V> & phi, const T v, const bool & check1D)
  {
    T base[3];
    V ainc = 0.0, rinc = 0.0, sol = EK_INF, minsol = EK_INF;

    if(_jun_nod_map[v]) {
      ainc = _antero_inc;
      rinc = _retro_inc;
    }

    // iterate over tetrahedra
    T start = _mesh.n2e_dsp[v], end = start + _mesh.n2e_cnt[v];

    for(T i=start, sidx=0; i<end; i++, sidx++)
    {
      T eidx  = _mesh.n2e_con[i];
      T esize = _mesh.e2n_cnt[eidx];
      switch(esize)
      {
        case 2:
          sol = solve_PS(phi, eidx, v) + ainc;
          break;
        case 3:
          sol = solve_2D(base, phi, eidx, v) + rinc;
          break;
        case 4:
          sol = solve_3D(base, phi, eidx, v, check1D) + rinc;
          break;
      }
      if(minsol > sol) minsol = sol;
    }

    return minsol;
  }

  inline void precompute_coeffs()
  {
    _sca_coeff.resize(_numElems * EK_COEFSZ);
    ek_double* coeff = _sca_coeff.data();

    ek_double M[6];
    ek_double e12[3], e13[3], e14[3], e23[3], e24[3], e34[3];
    const S* xyz = _mesh.xyz.data(), *xp;

    for(size_t idx=0; idx < _numElems; idx++)
    {
      if(_2lon) comp_M_2lon(M, idx);
      else comp_M_1lon(M, idx);

      const T* tet = _mesh.e2n_con.data() + _mesh.e2n_dsp[idx];

      xp = xyz + tet[0]*3;
      S x1_x = *xp++, x1_y = *xp++, x1_z = *xp++;
      xp = xyz + tet[1]*3;
      S x2_x = *xp++, x2_y = *xp++, x2_z = *xp++;
      xp = xyz + tet[2]*3;
      S x3_x = *xp++, x3_y = *xp++, x3_z = *xp++;
      xp = xyz + tet[3]*3;
      S x4_x = *xp++, x4_y = *xp++, x4_z = *xp++;

      e12[0] = x2_x - x1_x; e12[1] = x2_y - x1_y; e12[2] = x2_z - x1_z;
      e13[0] = x3_x - x1_x; e13[1] = x3_y - x1_y; e13[2] = x3_z - x1_z;
      e14[0] = x4_x - x1_x; e14[1] = x4_y - x1_y; e14[2] = x4_z - x1_z;
      e23[0] = x3_x - x2_x; e23[1] = x3_y - x2_y; e23[2] = x3_z - x2_z;
      e24[0] = x4_x - x2_x; e24[1] = x4_y - x2_y; e24[2] = x4_z - x2_z;
      e34[0] = x4_x - x3_x; e34[1] = x4_y - x3_y; e34[2] = x4_z - x3_z;

      coeff[0 ] = comp_3x3_aMb_symm(e12, M, e12);
      coeff[1 ] = comp_3x3_aMb_symm(e12, M, e23);
      coeff[2 ] = comp_3x3_aMb_symm(e12, M, e24);
      coeff[3 ] = comp_3x3_aMb_symm(e13, M, e13);
      coeff[4 ] = comp_3x3_aMb_symm(e13, M, e23);
      coeff[5 ] = comp_3x3_aMb_symm(e13, M, e34);
      coeff[6 ] = comp_3x3_aMb_symm(e14, M, e14);
      coeff[7 ] = comp_3x3_aMb_symm(e14, M, e24);
      coeff[8 ] = comp_3x3_aMb_symm(e14, M, e34);
      coeff[9 ] = comp_3x3_aMb_symm(e23, M, e23);
      coeff[10] = comp_3x3_aMb_symm(e23, M, e34);
      coeff[11] = comp_3x3_aMb_symm(e24, M, e24);
      coeff[12] = comp_3x3_aMb_symm(e24, M, e34);
      coeff[13] = comp_3x3_aMb_symm(e34, M, e34);

      coeff += EK_COEFSZ;
    }
  }


public:
  eikonal_solver(const struct mt_meshdata & mesh,
                 const struct ek_initdata<T, V> & init):
                 _mesh(mesh),
                 _init(init),
                 _numNodes(mesh.n2e_cnt.size()),
                 _numElems(mesh.e2n_cnt.size()),
                 _antero_inc(init.antero_delay),
                 _retro_inc(init.retro_delay)
  {
    _fin_nod_map.resize(_numNodes, false);
    _act_nod_map.resize(_numNodes, false);
    _jun_nod_map.resize(_numNodes, false);

    _nod_upd.resize(_numNodes, 0);
    _act_nod.resize(_numNodes);

    _2lon = (_mesh.lon.size() / _numElems) == 6;

    // find junction nodes in mesh (nodes connecting tetrahedral and line elements)
    for(size_t i=0; i<_numNodes; i++)
    {
      T start = _mesh.n2e_dsp[i], end = start + _mesh.n2e_cnt[i];
      T eidx  = _mesh.n2e_con[start];
      T tval  = _mesh.e2n_cnt[eidx];

      for(int j = start+1; j < end; j++)
      {
        eidx = _mesh.n2e_con[j];
        T esize = _mesh.e2n_cnt[eidx];
        if(esize != tval) {
          _jun_nod_map[i] = true;
          break;
        }
      }
    }

    _vf  = init.vf;
    _vs  = init.vs;
    _vn  = init.vn;
    // convert default init velocities (which are m/s) to um/ms
    // _vf = init.vf * 1e3;
    // _vs = init.vs * 1e3;
    // _vn = init.vn * 1e3;
    // _vps = init.vps * 1e3;

    // _vf = 1.0 / _vf / _vf;
    // _vs = 1.0 / _vs / _vs;
    // _vn = 1.0 / _vn / _vn;
    // _vps = 1.0 / _vps / _vps;

    _initval = EK_INF;

    // precompute local solver coefficients
    precompute_coeffs();
  }


  inline void operator()(mt_vector<V> & phi)
  {
    phi.assign(_numNodes, _initval);

    // initialize active nodes and phi
    for(size_t i=0; i<_init.nod.size(); i++)
    {
      T nidx = _init.nod[i];
      _fin_nod_map[nidx] = true;
      phi[nidx] = _init.val[i];
    }
    for(const T nidx : _init.nod) {
      mt_int start = _mesh.n2n_dsp[nidx], end = start + _mesh.n2n_cnt[nidx];
      for(int j=start; j<end; j++) {
        T idx = _mesh.n2n_con[j];
        if(_fin_nod_map[idx] == false)
          _act_nod_map[idx] = true;
      }
    }

    // set up active list
    _act_idx = 0;
    for(size_t i=0; i<_numNodes; i++)
      if(_act_nod_map[i]) _act_nod[_act_idx++] = i;

    // solve for active nodes ---------------------------------------------
    size_t loopidx = 0;
    size_t numConv;
    PROGRESS<size_t> prg(_numNodes, "Eikonal solve progress: ");

    while( _act_idx > 0 )
    {
      loopidx ++;
      numConv = 0;

      // if( !(loopidx % 10) )
      //   std::cout << "Iter: " << loopidx << ", active list size: " << _act_idx << ", "
      //             << ((float)_act_idx) / _numNodes * 100 << " % of all nodes" << std::endl;

      size_t asize = _act_idx;
      // compute phi
      #ifdef OPENMP
      #pragma omp parallel for schedule(guided)
      #endif
      for(size_t i=0; i<asize; i++)
      {
        T actidx = _act_nod[i];
        V oldval, newval;

        oldval = phi[actidx];
        newval = solve_nbh(phi, actidx, false);

        if(newval < oldval) {
          phi[actidx] = newval;
          _nod_upd[actidx]++;
        }
        else
          newval = oldval;


        if( fabs(newval - oldval) < (EK_ACC * (1 + (newval + oldval)*0.5)) )
        {
          nodal_expand(phi, actidx, newval);
          _act_nod_map[actidx] = false;
          numConv++;
        }
      }

      _act_idx = 0;
      for(size_t i=0; i<_numNodes; i++)
        if(_act_nod_map[i]) _act_nod[_act_idx++] = i;

      // we estimate that 80% of the nodes removed from the active list
      // won't enter the list again. Therefore this amount makes up the
      // progress of the current iteration.
      prg.increment(size_t(numConv * 0.8f));
    }
    prg.finish();
  }

  inline const mt_vector<T> & get_upd() const
  {
    return _nod_upd;
  }

};


