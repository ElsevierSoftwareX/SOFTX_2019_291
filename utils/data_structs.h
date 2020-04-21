/**
* @file data_structs.h
* @brief The main meshtool data structs.
* @author Aurel Neic
* @version
* @date 2016-12-13
*/

#ifndef _DATA_STRUCTS_H
#define _DATA_STRUCTS_H

#include <cmath>
#include <cstring>
#include <cstdlib>

#include "mt_vector.h"
#include "hashmap.hpp"

typedef MT_INT mt_int;
typedef MT_REAL mt_real;

typedef enum
{
  Tetra,
  Hexa,
//  Octa,
  Pyramid,
  Prism,
  Quad,
  Tri,
  Line
} elem_t;

inline short
get_nnodes(elem_t type)
{
  switch(type) {
    case Line:    return 2;
    case Tri:     return 3;
    case Quad:
    case Tetra:   return 4;
    case Pyramid: return 5;
    case Prism:   return 6;
    case Hexa:    return 8;
  }

  return -1;
}

/**
* @brief The graph of a (sub)mesh
*/
struct mt_meshgraph
{
  mt_vector<elem_t> etype;
  mt_vector<mt_int> e2n_cnt;
  mt_vector<mt_int> e2n_con;
  mt_vector<mt_int> eidx;
};

/**
* @brief Struct needed for processing data for Neumann BC files.
*/
struct nbc_data {
  mt_vector<mt_int> eidx;
  mt_vector<mt_int> sp_vtx;
  mt_vector<mt_int> tag;
};


/**
* @brief The main mesh struct.
*
* Note the following naming convention:
*  e2n: graph connecting elements to nodes
*  n2e: graph connecting nodes to elements, thus e2n transposed
*  n2n: graph connecting nodes with nodes
*/
struct mt_meshdata
{
  mt_vector<mt_int> e2n_cnt;
  mt_vector<mt_int> e2n_dsp;
  mt_vector<mt_int> e2n_con;

  mt_vector<mt_int> n2e_cnt;
  mt_vector<mt_int> n2e_dsp;
  mt_vector<mt_int> n2e_con;

  mt_vector<mt_int> n2n_cnt;
  mt_vector<mt_int> n2n_dsp;
  mt_vector<mt_int> n2n_con;

  mt_vector<elem_t> etype;
  mt_vector<mt_int> etags;
  mt_vector<mt_real> xyz;
  mt_vector<mt_real> lon;
};

class mt_manifold {
  public:
  mt_meshdata mesh;
  mt_mask     on_mnfld;

  inline void set_mask() {
    on_mnfld.insert(mesh.e2n_con.begin(), mesh.e2n_con.end());
  }
};


/**
* @brief Clamp a value into an interval [start, end]
*
* @tparam V   Value type
* @tparam W   Interval boundary type
* @param val    The value we clamp
* @param start  The interval start value
* @param end    The interval end value
*
* @return The clamped value
*/
template<typename V, typename W> inline
V clamp(const V val, const W start, const W end, bool warn = false)
{
  bool outside_clamp_region = false;
  V ret = val;

  if(val < start) {
    ret = start;
    outside_clamp_region = true;
  }
  else if(val > end) {
    ret = end;
    outside_clamp_region = true;
  }
  static size_t cnt = 0;
  if(warn && outside_clamp_region) {
    if(cnt < 20)
      fprintf(stderr, "Warning: Clamped value %g outside clamping region [%g, %g]!\n", static_cast<mt_real>(val), static_cast<mt_real>(start), static_cast<mt_real>(end));
    else if(cnt == 20)
      fprintf(stderr, "Warning: Further clamping warnings will be suppressed!\n");
    cnt++;
  }
  return ret;
}

/**
* @brief Get an estimate for the decimal power of a value.
*
* This may be useful for estimate how to scale a mesh etc..
*
*/
template<typename V> inline
int get_dec_power_estimate(V val)
{
  double dec_log = log(val) / log(10);

  // we want to round towards negative infinity also for negative values
  int ret = dec_log < 0 ? (dec_log - 1.0) : dec_log;

  return ret;
}

/**
* @brief A class that stores the graph of a mapping a -> b,
* its transposed (backward) mapping b -> a, a -> a and b -> b.
*
* @tparam T  Index data type.
*/
template<class T>
struct mt_mapping
{
  mt_vector<T> fwd_cnt;
  mt_vector<T> fwd_dsp;
  mt_vector<T> fwd_con;

  mt_vector<T> bwd_cnt;
  mt_vector<T> bwd_dsp;
  mt_vector<T> bwd_con;

  mt_vector<T> a2a_cnt;
  mt_vector<T> a2a_dsp;
  mt_vector<T> a2a_con;

  mt_vector<T> b2b_cnt;
  mt_vector<T> b2b_dsp;
  mt_vector<T> b2b_con;

  /// create b -> a from a -> b
  void transpose()
  {
    transpose_connectivity(fwd_cnt, fwd_con, bwd_cnt, bwd_con, false);
  }

  void setup_a2a()
  {
    if(bwd_cnt.size() == 0)
      transpose();

    multiply_connectivities(fwd_cnt, fwd_con, bwd_cnt, bwd_con, a2a_cnt, a2a_con);
  }
  void setup_b2b()
  {
    if(bwd_cnt.size() == 0)
      transpose();

    multiply_connectivities(bwd_cnt, bwd_con, fwd_cnt, fwd_con, b2b_cnt, b2b_con);
  }

  /// update displacement vectors
  void setup_dsp()
  {
    fwd_dsp.resize(fwd_cnt.size());
    bucket_sort_offset(fwd_cnt, fwd_dsp);

    bwd_dsp.resize(bwd_cnt.size());
    bucket_sort_offset(bwd_cnt, bwd_dsp);

    a2a_dsp.resize(a2a_cnt.size());
    bucket_sort_offset(a2a_cnt, a2a_dsp);

    b2b_dsp.resize(b2b_cnt.size());
    bucket_sort_offset(b2b_cnt, b2b_dsp);
  }

  void resize(size_t n)
  {
    if(fwd_cnt.size()) fwd_cnt.resize(n), fwd_cnt.reallocate();
    if(fwd_dsp.size()) fwd_dsp.resize(n), fwd_dsp.reallocate();
    if(fwd_con.size()) fwd_con.resize(n), fwd_con.reallocate();
    if(bwd_cnt.size()) bwd_cnt.resize(n), bwd_cnt.reallocate();
    if(bwd_dsp.size()) bwd_dsp.resize(n), bwd_dsp.reallocate();
    if(bwd_con.size()) bwd_con.resize(n), bwd_con.reallocate();
    if(a2a_cnt.size()) a2a_cnt.resize(n), a2a_cnt.reallocate();
    if(a2a_dsp.size()) a2a_dsp.resize(n), a2a_dsp.reallocate();
    if(a2a_con.size()) a2a_con.resize(n), a2a_con.reallocate();
    if(b2b_cnt.size()) b2b_cnt.resize(n), b2b_cnt.reallocate();
    if(b2b_dsp.size()) b2b_dsp.resize(n), b2b_dsp.reallocate();
    if(b2b_con.size()) b2b_con.resize(n), b2b_con.reallocate();
  }
};

struct mt_pscable
{
  mt_vector<mt_real> pts;
  mt_int par1, par2;
  mt_int br1, br2;
  mt_real size, rgj, cond;
};

struct mt_psdata {

  mt_vector<mt_pscable> cables;

};


template<class T>
struct tuple {
  T v1;
  T v2;
};
template<class T, class S>
struct mixed_tuple
{
  T v1;
  S v2;
};


template<class T>
struct triple {
  T v1;
  T v2;
  T v3;
};

template<class T, class S, class V>
struct mixed_triple {
  T v1;
  S v2;
  V v3;
};

template<class T>
struct quadruple {
  T v1;
  T v2;
  T v3;
  T v4;
};

template<class T>
bool operator<(const struct tuple<T> & lhs, const struct tuple<T> & rhs)
{
  if(lhs.v1 != rhs.v1)
    return lhs.v1 < rhs.v1;
  else
    return lhs.v2 < rhs.v2;
}
template<class T>
bool operator==(const struct tuple<T> & lhs, const struct tuple<T> & rhs)
{
    return lhs.v1 == rhs.v1 && lhs.v2 == rhs.v2;
}

template<class T, class S>
bool operator<(const struct mixed_tuple<T, S> & lhs, const struct mixed_tuple<T, S> & rhs)
{
  if(lhs.v1 != rhs.v1)
    return lhs.v1 < rhs.v1;
  else
    return lhs.v2 < rhs.v2;
}

template<class T>
bool operator<(const struct triple<T> & lhs, const struct triple<T> & rhs)
{
  if(lhs.v1 != rhs.v1)
    return lhs.v1 < rhs.v1;
  else if(lhs.v2 != rhs.v2)
    return lhs.v2 < rhs.v2;
  else
    return lhs.v3 < rhs.v3;
}

template<class T>
bool operator<(const struct quadruple<T> & lhs, const struct quadruple<T> & rhs)
{
  if(lhs.v1 != rhs.v1)
    return lhs.v1 < rhs.v1;
  else if(lhs.v2 != rhs.v2)
    return lhs.v2 < rhs.v2;
  else if(lhs.v3 != rhs.v3)
    return lhs.v3 < rhs.v3;
  else
    return lhs.v4 < rhs.v4;
}

template<class T>
short common_components(const triple<T> & lhs, const triple<T> & rhs)
{
  short c = 0;
  c += short(lhs.v1 == rhs.v1 || lhs.v1 == rhs.v2 || lhs.v1 == rhs.v3);
  c += short(lhs.v2 == rhs.v1 || lhs.v2 == rhs.v2 || lhs.v2 == rhs.v3);
  c += short(lhs.v3 == rhs.v1 || lhs.v3 == rhs.v2 || lhs.v3 == rhs.v3);

  return c;
}

// expand hashing for custom triple struct
namespace hashmap {

template<typename T>
struct hash_ops< tuple<T> >
{
  static inline bool cmp(tuple<T> a, tuple<T> b) {
    return a.v1 == b.v1 && a.v2 == b.v2;
  }

  static inline hm_uint hash(tuple<T> a) {
    hm_uint h = mkhash_init;
    h = mkhash(h, hash_ops<T>::hash(a.v1));
    h = mkhash(h, hash_ops<T>::hash(a.v2));
    return h;
  }
};

template<typename T>
struct hash_ops< triple<T> >
{
  static inline bool cmp(triple<T> a, triple<T> b) {
    return a.v1 == b.v1 && a.v2 == b.v2 && a.v3 == b.v3;
  }

  static inline hm_uint hash(triple<T> a) {
    hm_uint h = mkhash_init;
    h = mkhash(h, hash_ops<T>::hash(a.v1));
    h = mkhash(h, hash_ops<T>::hash(a.v2));
    h = mkhash(h, hash_ops<T>::hash(a.v3));
    return h;
  }
};
// expand hashing for custom triple struct
template<typename T>
struct hash_ops< quadruple<T> >
{
  static inline bool cmp(quadruple<T> a, quadruple<T> b) {
    return a.v1 == b.v1 && a.v2 == b.v2 && a.v3 == b.v3 && a.v4 == b.v4;
  }

  static inline hm_uint hash(quadruple<T> a) {
    hm_uint h = mkhash_init;
    h = mkhash(h, hash_ops<T>::hash(a.v1));
    h = mkhash(h, hash_ops<T>::hash(a.v2));
    h = mkhash(h, hash_ops<T>::hash(a.v3));
    h = mkhash(h, hash_ops<T>::hash(a.v4));
    return h;
  }
};
}

struct edgeele {
  mt_int v1;
  mt_int v2;
  mt_int idx;
  mt_real length;
  float split_at;
};

struct edge_len_dsc {
  bool operator() (const edgeele & lhs, const edgeele & rhs) const
  {
    if (lhs.length != rhs.length) return lhs.length > rhs.length;
    else                          return lhs.idx > rhs.idx;
  }
};
struct edge_len_asc {
  bool operator() (const edgeele & lhs, const edgeele & rhs) const
  {
    if (lhs.length != rhs.length) return lhs.length < rhs.length;
    else                          return lhs.idx < rhs.idx;
  }
};


struct tri_sele {
  mt_int v1;
  mt_int v2;
  mt_int v3;
  mt_int eidx;
};

struct quad_sele {
  mt_int v1;
  mt_int v2;
  mt_int v3;
  mt_int v4;
  mt_int eidx;
};

/**
* @brief A point class to simplify geometric computations.
*
* Only use this for local computations. Do not store an array of mt_points.
*
* @tparam S floating point type.
*/
template<class S>
class mt_point
{
public:
  S x, y, z;

  mt_point(): x(S()), y(S()), z(S())
  {}

  mt_point(S ix, S iy, S iz): x(ix), y(iy), z(iz)
  {}

  template<class V>
  mt_point(const mt_point<V> & vec): x(vec.x), y(vec.y), z(vec.z)
  {}

  template<class V>
  mt_point(const V* pts)
  {
    x = pts[0];
    y = pts[1];
    z = pts[2];
  }

  template<class V>
  void operator= (const mt_point<V> & pt)
  {
    x = pt.x;
    y = pt.y;
    z = pt.z;
  }

  void operator+= (const mt_point<S> & pt)
  {
    x += pt.x;
    y += pt.y;
    z += pt.z;
  }
  void operator-= (const mt_point<S> & pt)
  {
    x -= pt.x;
    y -= pt.y;
    z -= pt.z;
  }
  void operator*= (S s)
  {
    x *= s;
    y *= s;
    z *= s;
  }
  void operator/= (S s)
  {
    x /= s;
    y /= s;
    z /= s;
  }

  /// return vector length
  S length() const
  {
    S len = sqrt(x*x + y*y + z*z);
    return len;
  }
  /// return squared vector length
  S length2() const
  {
    return (x*x + y*y + z*z);
  }

  mt_point<S> crossProd(const mt_point<S> & v) const
  {
    return {(y * v.z) - (z * v.y), (z * v.x) - (x * v.z), (x * v.y) - (y * v.x)};
  }

  S scaProd(const mt_point<S> & v) const
  {
    return x*v.x + y*v.y + z*v.z;
  }

  void normalize()
  {
    S len = this->length();
    x /= len;
    y /= len;
    z /= len;
  }

  template<class V>
  void get(const V* loc)
  {
    x = loc[0];
    y = loc[1];
    z = loc[2];
  }

  template<class V>
  void set(V* loc) const
  {
    loc[0] = x;
    loc[1] = y;
    loc[2] = z;
  }
};

template<class S>
mt_point<S> operator+ (const mt_point<S> & lhs, const mt_point<S> & rhs)
{
  return mt_point<S>(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
}
template<class S>
mt_point<S> operator- (const mt_point<S> & lhs, const mt_point<S> & rhs)
{
  return mt_point<S>(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
}

template<class S>
mt_point<S> operator* (const mt_point<S> & lhs, const S s)
{
  return mt_point<S>(lhs.x * s, lhs.y * s, lhs.z * s);
}

template<class S>
mt_point<S> operator/ (const mt_point<S> & lhs, const S s)
{
  return mt_point<S>(lhs.x / s, lhs.y / s, lhs.z / s);
}

template<class S, class V>
void array_to_points(const mt_vector<S> & xyz, mt_vector<mt_point<V> > & p)
{
  p.resize(xyz.size()/3);
  for(size_t i=0; i<p.size(); i++)
  {
    p[i].x = xyz[i*3+0];
    p[i].y = xyz[i*3+1];
    p[i].z = xyz[i*3+2];
  }
}
template<class S, class V>
void points_to_array(const mt_vector<mt_point<S> > & p, mt_vector<V> & xyz)
{
  xyz.resize(p.size()*3);
  for(size_t i=0; i<p.size(); i++)
  {
    xyz[i*3+0] = p[i].x;
    xyz[i*3+1] = p[i].y;
    xyz[i*3+2] = p[i].z;
  }
}

template<class S>
mt_point<S> unit_vector(const mt_point<S> & v)
{
  S len = v.length();
  return mt_point<S>(v.x / len, v.y / len, v.z / len);
}

template<class T, class V>
inline mt_point<T> cast_vec(const V & vec)
{
  mt_point<T> ret;
  ret.x = vec.x, ret.y = vec.y, ret.z = vec.z;

  return ret;
}


typedef mt_point<float>   vec3f;
typedef mt_point<mt_real> vec3r;
typedef mt_point<mt_int>  vec3i;
// wrapper routine to switch between fast inverse sqrt algorithm (for floats) or normal
float inv_sqrtf(const float number );
// wrapper routine to switch between fast inverse sqrt algorithm (for doubles) or normal
double inv_sqrtd(const double number);

//Explict specialization for float values with fast inverse sq algorithm
template<> inline
void mt_point<float>::normalize() {
  const float invsq = inv_sqrtf(this->length2());//fast_inv_sqrtf(lensq);
  this->x *= invsq;
  this->y *= invsq;
  this->z *= invsq;
}

//Explict specialization for double values with fast inverse sq algorithm
template<> inline
void mt_point<double>::normalize() {
  const double invsq = inv_sqrtd(this->length2());
  this->x *= invsq;
  this->y *= invsq;
  this->z *= invsq;
}

//Explicit specialization for floats
template<> inline
mt_point<float> unit_vector(const mt_point<float> & v)
{
  const float invsq = inv_sqrtf(v.length2());
  return mt_point<float>(v.x * invsq, v.y * invsq, v.z * invsq);
}

template<> inline
mt_point<double> unit_vector(const mt_point<double> & v)
{
  const double invsq = inv_sqrtd(v.length2());
  return mt_point<double>(v.x * invsq, v.y * invsq, v.z * invsq);
}

#endif
