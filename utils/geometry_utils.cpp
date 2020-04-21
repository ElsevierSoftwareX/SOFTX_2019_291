/**
* @file geometry_utils.cpp
* @brief Utility functions for geometric operations.
* @author Aurel Neic
* @version
* @date 2018-10-28
*/

#include "mt_utils_base.h"
#include "geometry_utils.h"
#include "dense_mat.hpp"
#include <limits>
#include <vector>
#include <cmath>

//functions are defined in utils/tetgen/predicates.cxx
//CAREFUL the last parameter is the base point for calculating the volume
//So in our orientation the p0 point needs to be always in the last slot
extern double orient3d(double*, double*, double*, double*);
extern double orient2d(double*, double*, double*);

//Compute angle
mt_real angle_between_line(const vec3r & A, const vec3r & B, const vec3r & C)
{
  mt_real a = (B-A).length();
  mt_real b = (C-A).length();
  const mt_real c = (C-B).length();

  const mt_real tmp = a;
  if(b > a) {
   a = b;
   b = tmp;
  }
  mt_real mu = 0.0;
  if(b >= c && c >= 0.0) {
    mu = c - (a - b);
  } else if(c > b && b >= 0.0) {
    mu = b - (a - c);
  } else {
    fprintf(stderr, "%s error: invalid triangle provided!\n", __func__);
    exit(EXIT_FAILURE);
  }
  const mt_real fac = sqrt(((a-b)+c)*mu / ((a+(b+c))*((a-c)+b)));
  return 2.0*atan(fac);
}


vec3r triangle_normal(const mt_int v1, const mt_int v2, const mt_int v3,
                      const mt_vector<mt_real> & xyz)
{
  vec3r e1 = vec3r(xyz.data() + v2*3) - vec3r(xyz.data() + v1*3);
  vec3r e2 = vec3r(xyz.data() + v3*3) - vec3r(xyz.data() + v1*3);

  vec3r ret = e1.crossProd(e2) * mt_real(-1.0);
  ret.normalize();

  return ret;
}

vec3r triangle_normal(const vec3r & p0, const vec3r & p1, const vec3r & p2)
{
  vec3r e1 = p1 - p0;
  vec3r e2 = p2 - p0;

  vec3r ret = e1.crossProd(e2) * mt_real(-1.0);
  ret.normalize();

  return ret;
}

vec3r triangle_centerpoint(const mt_int v1, const mt_int v2, const mt_int v3,
                           const mt_vector<mt_real> & xyz)
{
  vec3r p1(xyz.data() + v1*3), p2(xyz.data() + v2*3), p3(xyz.data() + v3*3);

  vec3r ret = p1 + p2 + p3;
  ret /= mt_real(3.0);

  return ret;
}


void fit_sphere(const mt_vector<vec3r> & vert, mt_real & R, vec3r & pos)
{
  int nvert = vert.size();
  vec3r ctr(0,0,0);

  // compute centroid of vertex set
  for(const vec3r & v : vert) ctr += v;
  ctr /= mt_real(nvert);

  dmat<mt_real> A(3,3), B(3,1);
  A.assign(0.0), B.assign(0.0);

  // assemble A and B matrices
  for(int i=0; i<nvert; i++) {
    A[0][0] += vert[i].x * (vert[i].x - ctr.x);
    A[0][1] += vert[i].x * (vert[i].y - ctr.y);
    A[0][2] += vert[i].x * (vert[i].z - ctr.z);

    A[1][0] += vert[i].y * (vert[i].x - ctr.x);
    A[1][1] += vert[i].y * (vert[i].y - ctr.y);
    A[1][2] += vert[i].y * (vert[i].z - ctr.z);

    A[2][0] += vert[i].z * (vert[i].x - ctr.x);
    A[2][1] += vert[i].z * (vert[i].y - ctr.y);
    A[2][2] += vert[i].z * (vert[i].z - ctr.z);

    B[0][0] += vert[i].length2() * (vert[i].x - ctr.x);
    B[1][0] += vert[i].length2() * (vert[i].y - ctr.y);
    B[2][0] += vert[i].length2() * (vert[i].z - ctr.z);
  }

  A *= 2.0 / mt_real(nvert);
  B *= 1.0 / mt_real(nvert);

  #if 1
  // compute A^T
  dmat<mt_real> At = A; At.transpose();

  // compute (A^T . A) ^ -1
  mt_real det;
  dmat<mt_real> AtA_inv = At * A;
  invert_3x3(AtA_inv, det);

  // [x,y,z] = ((A^T . A) ^ -1) . (A^T . B)
  dmat<mt_real> rhs = At * B, sol = AtA_inv * rhs;
  pos.x = sol[0][0], pos.y = sol[1][0], pos.z = sol[2][0];
  #else
  // just compute [x,y,z] = A^-1 . B
  mt_real det;
  invert_3x3(A, det);
  dmat<mt_real> sol = A * B;
  pos.x = sol[0][0], pos.y = sol[1][0], pos.z = sol[2][0];
  #endif

  // compute R = sqrt( frac{sum_(i = 0)^(n-1) <vert[i] - pos, vert[i] - pos> } {nvert} )
  for(int i=0; i<nvert; i++) {
    vec3r d = vert[i] - pos;
    R += d.length2();
  }
  R = sqrt(R / mt_real(nvert));
}

mt_real circumsphere(const elem_t type, const mt_int* con, const mt_real* xyz, vec3r & ctr)
{
  mt_vector<vec3r> pts(get_nnodes(type));

  for(size_t i=0; i<pts.size(); i++)
  {
    mt_int idx = con[i];
#ifdef OPENMP
#pragma omp atomic read
#endif
    pts[i].x = xyz[idx*3 + 0];
#ifdef OPENMP
#pragma omp atomic read
#endif
    pts[i].y = xyz[idx*3 + 1];
#ifdef OPENMP
#pragma omp atomic read
#endif
    pts[i].z = xyz[idx*3 + 2];
  }

  switch(type)
  {
    case Line:  return line_circumcircle(pts.data(), ctr);
    case Tri:   return tri_circumcircle(pts.data(), ctr);
    case Tetra: return tet_circumsphere(pts.data(), ctr);
    default:
    {
      fprintf(stderr, "ERROR in %s: Unsupported element type. Aborting!\n", __func__);
      exit(EXIT_FAILURE);
    }
  }
}

mt_real volume(const elem_t type, const mt_int* con, const mt_real* xyz)
{
  mt_vector<vec3r> pts(get_nnodes(type));

  for(size_t i=0; i<pts.size(); i++)
  {
    mt_int idx = con[i];
#ifdef OPENMP
#pragma omp atomic read
#endif
    pts[i].x = xyz[idx*3 + 0];
#ifdef OPENMP
#pragma omp atomic read
#endif
    pts[i].y = xyz[idx*3 + 1];
#ifdef OPENMP
#pragma omp atomic read
#endif
    pts[i].z = xyz[idx*3 + 2];
  }

  switch(type)
  {
    case Line:    return (pts[0]-pts[1]).length();
    case Tri:     return tri_surf    (pts.data());
    case Quad:    return quad_surf   (pts.data());
    case Tetra:   return tet_volume  (pts.data());
    case Pyramid: return pyr_volume  (pts.data());
    case Prism:   return prism_volume(pts.data());
    case Hexa:    return hex_volume  (pts.data());
    default:
      fprintf(stderr, "%s error, unsupported element type. Aborting!\n", __func__);
      exit(EXIT_FAILURE);
  }
}

vec3r barycenter(const mt_int npts, const mt_int* con, const mt_real* xyz)
{
  vec3r bary(0,0,0), p;

  for(mt_int i=0; i < npts; i++) {
    mt_int c = con[i];
#ifdef OPENMP
#pragma omp atomic read
#endif
    p.x = xyz[c*3+0];
#ifdef OPENMP
#pragma omp atomic read
#endif
    p.y = xyz[c*3+1];
#ifdef OPENMP
#pragma omp atomic read
#endif
    p.z = xyz[c*3+2];
    bary += p;
  }

  bary /= mt_real(npts);
  return bary;
}


mt_real tri_surf(vec3r & p0, vec3r & p1, vec3r & p2)
{
  vec3r e1 = p1 - p0, e2 = p2 - p0;
  vec3r n = e1.crossProd(e2);
  return n.length() * 0.5;
}

mt_real tri_surf(vec3r * points)
{
  return tri_surf(points[0], points[1], points[2]);
}

mt_real quad_surf(vec3r * points)
{
  vec3r ctr = (points[0] + points[1] + points[2] + points[4]) * 0.25;
  mt_real s = 0.0;

  s += tri_surf(points[0], points[1], ctr);
  s += tri_surf(points[1], points[2], ctr);
  s += tri_surf(points[2], points[3], ctr);
  s += tri_surf(points[3], points[0], ctr);

  return s;
}

mt_real tet_volume(vec3r * points) {
  return tet_volume(points[0], points[1], points[2], points[3]);
}

mt_real signed_tet_volume(const vec3r & p0,
                          const vec3r & p1,
                          const vec3r & p2,
                          const vec3r & p3)
{
#ifdef EXACT_PREDICATES
  mt_real p0arr[3] = {p0.x, p0.y, p0.z};
  mt_real p1arr[3] = {p1.x, p1.y, p1.z};
  mt_real p2arr[3] = {p2.x, p2.y, p2.z};
  mt_real p3arr[3] = {p3.x, p3.y, p3.z};
  const mt_real vol = orient3d(p1arr, p2arr, p3arr, p0arr) / 6.;
#else
  vec3r e1 = p1 - p0, e2 = p2 - p0, e3 = p3 - p0;
  vec3r n = e1.crossProd(e2);
  const mt_real vol = n.scaProd(e3) / 6.;
#endif
  return vol;
}

mt_real signed_tet_volume(vec3r * points) {
  return signed_tet_volume(points[0], points[1], points[2], points[3]);
}

mt_real signed_tri_surf(const vec3r & p0,
                        const vec3r & p1,
                        const vec3r & p2)
{
  mt_real p0arr[3] = {p0.x, p0.y, p0.z};
  mt_real p1arr[3] = {p1.x, p1.y, p1.z};
  mt_real p2arr[3] = {p2.x, p2.y, p2.z};

  mt_real vol = orient2d(p1arr, p2arr, p0arr) / 2.;

  return vol;
}

mt_real signed_tri_surf(vec3r * points) {
  return signed_tri_surf(points[0], points[1], points[2]);
}

mt_real tet_volume(const vec3r & p0,
                   const vec3r & p1,
                   const vec3r & p2,
                   const vec3r & p3)
{

  return fabs(signed_tet_volume(p0,p1,p2,p3));
}

mt_real pyr_volume(vec3r * points)
{
  // center point of the pyramid quad face
  vec3r fc = (points[0] + points[1] + points[2] + points[3]) * mt_real(0.25);

  mt_real pyr_vol = 0.;

  // we construct a tet for each edge. the third point of the tet base is the
  // center points of the pyramid's base.
  pyr_vol += tet_volume(points[0], points[1], fc, points[4]);
  pyr_vol += tet_volume(points[1], points[2], fc, points[4]);
  pyr_vol += tet_volume(points[2], points[3], fc, points[4]);
  pyr_vol += tet_volume(points[3], points[0], fc, points[4]);

  return pyr_vol;
}

mt_real hex_volume(vec3r * pts)
{
  vec3r ctr(0,0,0);
  for(mt_int i=0; i < 8; i++) ctr += pts[i];
  ctr *= (1. / 8.);

  mt_real hexvol = 0.;
  vec3r ppyr[5]; ppyr[4] = ctr;

  // we construct a pyramid for each face. the pyramid tip is the centerpoint.
  // face (0,1,2,3)
  ppyr[0] = pts[0]; ppyr[1] = pts[1];
  ppyr[2] = pts[2]; ppyr[3] = pts[3];
  hexvol += pyr_volume(ppyr);

  // face (0,4,7,1)
  ppyr[0] = pts[0]; ppyr[1] = pts[4];
  ppyr[2] = pts[7]; ppyr[3] = pts[1];
  hexvol += pyr_volume(ppyr);

  // face (1,7,6,2)
  ppyr[0] = pts[1]; ppyr[1] = pts[7];
  ppyr[2] = pts[6]; ppyr[3] = pts[2];
  hexvol += pyr_volume(ppyr);

  // face (2,6,5,3)
  ppyr[0] = pts[2]; ppyr[1] = pts[6];
  ppyr[2] = pts[5]; ppyr[3] = pts[3];
  hexvol += pyr_volume(ppyr);

  // face (3,5,4,0)
  ppyr[0] = pts[3]; ppyr[1] = pts[5];
  ppyr[2] = pts[4]; ppyr[3] = pts[0];
  hexvol += pyr_volume(ppyr);

  // face (4,5,6,7)
  ppyr[0] = pts[4]; ppyr[1] = pts[5];
  ppyr[2] = pts[6]; ppyr[3] = pts[7];
  hexvol += pyr_volume(ppyr);

  return hexvol;
}

mt_real prism_volume(vec3r * pts)
{
  vec3r ctr(0,0,0);
  for(mt_int i=0; i < 6; i++) ctr += pts[i];
  ctr *= (1. / 6.);

  mt_real prismvol = 0.;
  vec3r ppyr[5]; ppyr[4] = ctr;
  vec3r ptet[4]; ptet[3] = ctr;

  // we construct a pyramid for each quad face. the pyramid tip is the centerpoint.
  // face (0,3,5,1)
  ppyr[0] = pts[0]; ppyr[1] = pts[3];
  ppyr[2] = pts[5]; ppyr[3] = pts[1];
  prismvol += pyr_volume(ppyr);

  // face (0,2,4,3)
  ppyr[0] = pts[0]; ppyr[1] = pts[2];
  ppyr[2] = pts[4]; ppyr[3] = pts[3];
  prismvol += pyr_volume(ppyr);

  // face (1,5,4,2)
  ppyr[0] = pts[1]; ppyr[1] = pts[5];
  ppyr[2] = pts[4]; ppyr[3] = pts[2];
  prismvol += pyr_volume(ppyr);

  // we construct a tet for each tri face. the tet tip is the centerpoint.
  // fase (0,1,2)
  ptet[0] = pts[0], ptet[1] = pts[1], ptet[2] = pts[2];
  prismvol += tet_volume(ptet);
  // fase (3,4,5)
  ptet[0] = pts[3], ptet[1] = pts[4], ptet[2] = pts[5];
  prismvol += tet_volume(ptet);

  return prismvol;
}

vec3r vertices_centerpoint(const mt_vector<mt_real> & xyz)
{
  mt_vector<vec3r> pts;
  array_to_points(xyz, pts);

  return vertices_centerpoint(pts);
}

vec3r vertices_centerpoint(const mt_vector<vec3r> & pts)
{
  vec3r s = sum(pts);
  s /= mt_real(pts.size());

  return s;
}

void bounding_sphere(const mt_vector<vec3r> & pts, vec3r & ctr, mt_real & rad)
{
  ctr = vertices_centerpoint(pts);
  rad = 0.0;

  for(const vec3r & v : pts) {
    mt_real l = (v - ctr).length2();
    if(rad < l) rad = l;
  }

  rad = sqrt(rad);
}

mt_real vertices_seperation_distance(const mt_vector<vec3r> & xyz, vec3r & p1, vec3r & p2)
{
  mt_real min_dist = std::numeric_limits<mt_real>::max();

  mt_int i_min = -1, i_max = -1;
  const mt_int sz = xyz.size();
  const mt_int loop_max = sz*(sz-1)/2;

  #ifdef OPENMP
  #pragma omp parallel for schedule(guided) reduction(min:min_dist)
  #endif
  for(mt_int k=0; k < loop_max; k++)
  {
    // we use the symmetry to only have a single loop
    mt_int i = k / sz, j = k % sz;
    if(j<=i) i = sz - i - 2, j = sz - j - 1;

    // distance between i and j vertices
    mt_real sqr_dist = (xyz[i] - xyz[j]).length2();

    if(sqr_dist <= min_dist) {
      min_dist = sqr_dist;
      i_min = i;
      i_max = j;
    }
  }

  if (i_min == -1 || i_max == -1)
    return (std::numeric_limits<mt_real>::max());

  p1 = xyz[i_min];
  p2 = xyz[i_max];
  min_dist = 0.5 * std::sqrt(min_dist);
  return (min_dist);
}

mt_real vertices_seperation_distance(const mt_vector<mt_real> & xyz, vec3r & p1, vec3r & p2)
{
  mt_vector<vec3r> pts;
  array_to_points(xyz, pts);
  return vertices_seperation_distance(pts, p1, p2);
}

mt_real vertices_seperation_distance(const mt_vector<vec3r> & xyz, const mt_vector<mt_int> & subidx, vec3r & p1, vec3r & p2)
{
  assert(subidx.size() > 1);
  mt_vector<vec3r> subset(subidx.size());

  for(size_t i=0; i < subidx.size(); i++)
    subset[i] = xyz[subidx[i]];

  return vertices_seperation_distance(subset, p1, p2);
}

mt_real vertices_seperation_distance(const mt_vector<mt_real> & xyz, const mt_vector<mt_int> & subidx, vec3r & p1, vec3r & p2)
{
  mt_vector<vec3r> pts;
  array_to_points(xyz, pts);
  return vertices_seperation_distance(pts, subidx, p1, p2);
}

void rotate_points(mt_vector<vec3r> & pts,
                   const mt_real rotx,
                   const mt_real roty,
                   const mt_real rotz)
{
  dmat<mt_real> rot_mat, rotx_mat, roty_mat, rotz_mat, ident(3,3);
  ident.assign(0.0); ident[0][0] = 1.0, ident[1][1] = 1.0, ident[2][2] = 1.0;

  if(rotx) {
    rotx_mat.set_size(3,3); rotx_mat.assign(0.0); rotx_mat[0][0] = 1.0;
    rotx_mat[1][1] = cos(rotx), rotx_mat[1][2] = -sin(rotx);
    rotx_mat[2][1] = sin(rotx), rotx_mat[2][2] =  cos(rotx);
  }
  else
    rotx_mat = ident;

  if(roty) {
    roty_mat.set_size(3,3); roty_mat.assign(0.0); roty_mat[1][1] = 1.0;
    roty_mat[0][0] =  cos(roty), roty_mat[0][2] = sin(roty);
    roty_mat[2][0] = -sin(roty), roty_mat[2][2] = cos(roty);
  }
  else
    roty_mat = ident;

  if(rotz) {
    rotz_mat.set_size(3,3); rotz_mat.assign(0.0); rotz_mat[2][2] = 1.0;
    rotz_mat[0][0] = cos(rotz), rotz_mat[0][1] = -sin(rotz);
    rotz_mat[1][0] = sin(rotz), rotz_mat[1][1] =  cos(rotz);
  }
  else
    rotz_mat = ident;

  rot_mat = rotx_mat * roty_mat * rotz_mat;

  mt_real ww[3], vv[3];
  for(vec3r & p : pts) {
    ww[0] = p.x, ww[1] = p.y, ww[2] = p.z;
    rot_mat.mult(ww, vv);
    p.x = vv[0], p.y = vv[1], p.z = vv[2];
  }
}

mt_real tet_circumsphere(const vec3r & p0, const vec3r & p1, const vec3r & p2, const vec3r & p3, vec3r & ctr)
{
  mt_real radius = 0.;
  dmat<mt_real> A(3,3);
  mt_real rhs[3];

  A[0][0] = p1.x - p0.x; A[0][1] = p1.y - p0.y; A[0][2] = p1.z - p0.z;
  A[1][0] = p2.x - p0.x; A[1][1] = p2.y - p0.y; A[1][2] = p2.z - p0.z;
  A[2][0] = p3.x - p0.x; A[2][1] = p3.y - p0.y; A[2][2] = p3.z - p0.z;

  rhs[0] = 0.5 * (A[0][0] * A[0][0] + A[0][1] * A[0][1] + A[0][2] * A[0][2]);
  rhs[1] = 0.5 * (A[1][0] * A[1][0] + A[1][1] * A[1][1] + A[1][2] * A[1][2]);
  rhs[2] = 0.5 * (A[2][0] * A[2][0] + A[2][1] * A[2][1] + A[2][2] * A[2][2]);

  const bool is_inv = A.lu_decomp();
  if(!is_inv) {
    fprintf(stderr, "ERROR in %s: Degenerate tetrahedron! Can't compute circumsphere!\n", __func__);
    exit(EXIT_FAILURE);
  }

  A.lu_solve(rhs);
  ctr.x = p0.x + rhs[0];
  ctr.y = p0.y + rhs[1];
  ctr.z = p0.z + rhs[2];

  radius = sqrt(rhs[0] * rhs[0] + rhs[1] * rhs[1] + rhs[2] * rhs[2]);
  return radius;
}

mt_real tet_circumsphere(const vec3r* points, vec3r & ctr)
{
  return tet_circumsphere(points[0], points[1], points[2], points[3], ctr);
}

mt_real tri_circumcircle(const vec3r & p0, const vec3r & p1, const vec3r & p2, vec3r & ctr)
{
  mt_real radius = 0.;
  dmat<mt_real> A(3,3);
  mt_real rhs[3];

  vec3r e10 = p1 - p0;
  vec3r e20 = p2 - p0;
  vec3r n = e10.crossProd(e20);

  A[0][0] = e10.x; A[0][1] = e10.y; A[0][2] = e10.z;
  A[1][0] = e20.x; A[1][1] = e20.y; A[1][2] = e20.z;
  A[2][0] = n.x; A[2][1] = n.y; A[2][2] = n.z;

  rhs[0] = 0.5 * (A[0][0] * A[0][0] + A[0][1] * A[0][1] + A[0][2] * A[0][2]);
  rhs[1] = 0.5 * (A[1][0] * A[1][0] + A[1][1] * A[1][1] + A[1][2] * A[1][2]);
  rhs[2] = 0.0;

  const bool is_inv = A.lu_decomp();
  if(!is_inv) {
    fprintf(stderr, "ERROR in %s: Degenerate triangle! Can't compute circumcircle!\n", __func__);
    exit(EXIT_FAILURE);
  }

  A.lu_solve(rhs);
  ctr.x = p0.x + rhs[0];
  ctr.y = p0.y + rhs[1];
  ctr.z = p0.z + rhs[2];

  radius = sqrt(rhs[0] * rhs[0] + rhs[1] * rhs[1] + rhs[2] * rhs[2]);
  return radius;
}

mt_real tri_circumcircle(const vec3r* points, vec3r & ctr)
{
  return tri_circumcircle(points[0], points[1], points[2], ctr);
}

mt_real line_circumcircle(const vec3r & p0, const vec3r & p1, vec3r & ctr)
{
  ctr = (p0 + p1) * 0.5;
  mt_real rad = (p1-ctr).length();

  return rad;
}

mt_real line_circumcircle(const vec3r* points, vec3r & ctr)
{
  return line_circumcircle(points[0], points[1], ctr);
}

mt_real PointCloudInterpolator::eval_rbf(const mt_real rsq, const mt_real ck) const
{
  const double cksq = ck * ck;
  return 1. / std::sqrt(cksq + rsq);
}

void PointCloudInterpolator::compute_roi(const bool verbose)
{
  _roi.resize(_pts.size(), 0.0);
  PROGRESS<size_t>* progress = NULL;
  if(verbose) progress = new PROGRESS<size_t>(_pts.size(), "PointCloud::ComputeRadiusOfInfluence: ");

  #ifdef OPENMP
  #pragma omp parallel
  #endif
  {
    mt_vector<mixed_tuple<mt_real, int> > res;

    #ifdef OPENMP
    #pragma omp for schedule(guided)
    #endif
    for(size_t k=0; k < _pts.size(); k++)
    {
      const vec3r actp = _pts[k];
      _tree->k_closest_vertices(_NW+1, actp, res);
      // k_closest_vertices inserts squared distances!!!
      _roi[k] = std::sqrt(res.end()->v1);

      if(verbose) {
        #ifdef OPENMP
        #pragma omp critical
        #endif
        progress->next();
      }
    }
  }
  if(verbose) {
    progress->finish();
    delete progress;
  }
}

void PointCloudInterpolator::compute_rbf_coeffs(const mt_vector<mt_real> & data, const bool verbose)
{
  const size_t sz = _pts.size();
  const mt_int loopmax = _NQ*(_NQ+1)/2;

  _rbf_scal.resize(sz, 0.0);
  _rbf_coeffs.resize(sz);
  _lu_mats.resize(sz);
  _nbh_idx.resize(sz);

  for(size_t i=0; i < sz; i++) {
    _rbf_coeffs[i].resize(_dpn * _NQ);
    _nbh_idx[i]   .resize(_NQ);
    _lu_mats[i]   .set_size(_NQ, _NQ);
  }

  if(_roi.size() == 0)
    compute_roi();

  PROGRESS<size_t>* progress = NULL;
  if(verbose) progress = new PROGRESS<size_t>(_pts.size(), "PointCloud::ComputeRadialBasisFunctionCoefficients: ");

  #ifdef OPENMP
  #pragma omp parallel
  #endif
  {
    mt_vector<mixed_tuple<mt_real, int> > res;
    mt_vector<vec3r> nbh_points(_NQ);
    mt_vector<mt_vector<mt_real> > rbf_rhs(_dpn);
    for(mt_int i=0; i < (mt_int)_dpn; i++)
      rbf_rhs[i].assign(_NQ, 0.0);

    #ifdef OPENMP
    #pragma omp for schedule(guided)
    #endif
    for(size_t k=0; k < _pts.size(); k++)
    {
      dmat<mt_real> & rbf_mat = _lu_mats[k];

      //Get the NQ-neighborhood of xk excluding xk
      const vec3r actp = _pts[k];
      _tree->k_closest_vertices(_NQ, actp, res);
      size_t widx = 0;
      for(auto it=res.begin(); it != res.end(); ++it)
        _nbh_idx[k][widx++] = it->v2;
      std::sort(_nbh_idx[k].begin(), _nbh_idx[k].end());

      for(mt_int i=0; i < _NQ; i++)
        nbh_points[i] = _pts[_nbh_idx[k][i]];

      _rbf_scal[k] = _roi[k];
      //Fill the rbf_mat
      // symmetric so we collapse the loop
      for(mt_int t=0; t < loopmax; t++) {
        size_t i = t/(_NQ+1), j = t%(_NQ+1);
        if(j>i) i = _NQ - i -1, j = _NQ - j;
        rbf_mat[i][j] = rbf_mat[j][i] = eval_rbf((nbh_points[i] - nbh_points[j]).length2(), _rbf_scal[k]);
      }
      //LU Decomp
      const bool is_inv = rbf_mat.lu_decomp();
      if(!is_inv) {
        fprintf(stderr, "ERROR in %s: Encountered non invertible RBF system!\n", __func__);
        exit(EXIT_FAILURE);
      }
      //Fill rhs and solve
      for(mt_int i=0; i < _dpn; i++)
      {
        for(mt_int l=0; l < _NQ; l++)
          rbf_rhs[i][l] = data[_dpn * _nbh_idx[k][l] + i];
        rbf_mat.lu_solve(rbf_rhs[i].data());
      }

      for(mt_int i=0; i < _dpn; i++)
        for(mt_int l=0; l < _NQ; l++)
          _rbf_coeffs[k][_dpn * l + i] = rbf_rhs[i][l];

      #ifdef OPENMP
      #pragma omp critical
      #endif
      if(verbose) progress->next();
    }
  }
  if(verbose) {
    progress->finish();
    delete progress;
  }
}

void PointCloudInterpolator::recompute_rbf_coeffs(const mt_vector<mt_real> & data, const bool verbose)
{
  assert(_lu_mats.size() > 0);
  PROGRESS<size_t> *progress = NULL;
  if(verbose)
    progress = new PROGRESS<size_t>(_pts.size(), "PointCloud::RecomputeRadialBasisFunctionCoefficients: ");

  #ifdef OPENMP
  #pragma omp parallel for schedule(guided)
  #endif
  for(size_t k=0; k < _pts.size(); k++)
  {
    mt_vector<mt_vector<mt_real> > rbf_rhs(_dpn);
    for(mt_int i=0; i < (mt_int)_dpn; i++)
      rbf_rhs[i].resize(_NQ);

    //Fill rhs and solve
    for(mt_int i=0; i < _dpn; i++)
    {
      for(mt_int l=0; l < _NQ; l++)
        rbf_rhs[i][l] = data[_dpn * _nbh_idx[k][l] + i];
      _lu_mats[k].lu_solve(rbf_rhs[i].data());
    }

    for(mt_int i=0; i < _dpn; i++)
      for(mt_int l=0; l < _NQ; l++)
        _rbf_coeffs[k][_dpn * l + i] = rbf_rhs[i][l];

    if(verbose) {
      #ifdef OPENMP
      #pragma omp critical
      #endif
      progress->next();
    }
  }
  if(verbose) {
    progress->finish();
    delete progress;
  }
}


void PointCloudInterpolator::interpolate_rbf(const mt_vector<mt_real> & inter_xyz, mt_vector<mt_real> & odat, const bool verbose) const
{
  if(_roi.size() + _rbf_coeffs.size() == 0) {
    fprintf(stderr, "ERROR in %s: Neither Radius of Influence nor RBF_Coeffs have been calculated!", __func__);
    exit(1);
  }

  mt_int nnodes = inter_xyz.size() / 3;
  odat.assign(_dpn * nnodes, 0.0);

  mt_vector<vec3r> inter_pts;
  array_to_points(inter_xyz, inter_pts);

  const float max_roi = *std::max_element(_roi.begin(), _roi.end());

  PROGRESS<size_t>* prg = NULL;
  if(verbose) prg = new PROGRESS<size_t>((size_t)nnodes, "RBF-Interpolation progress: ");

  #ifdef OPENMP
  #pragma omp parallel
  #endif
  {
    mt_vector<mixed_tuple<mt_real, int> > int_vtx;

    #ifdef OPENMP
    #pragma omp for schedule(guided)
    #endif
    for(mt_int nidx = 0; nidx < nnodes; nidx++)
    {
      vec3r pt = inter_pts[nidx];
      //check whether the point is inside of bounding sphere
      const bool is_in_sphere = (pt - _bounding_sphere_center).length2() < _bounding_sphere_radius_sq;
      if(is_in_sphere) {
        mt_real wsum = 0.;
        _tree->vertices_in_sphere(pt, 1.5 * max_roi, int_vtx);
        for(auto it = int_vtx.begin(); it != int_vtx.end(); ++it)
        {
          int k = it->v2;
          const mt_real rk = (pt-_pts[k]).length();
          const mt_real wk = rk < _roi[k] ? POW2((_roi[k] - rk) / (rk * _roi[k])) : 0.;
          if(wk == 0.0)
            continue;
          wsum += wk;
          //Eval the local radial basis function interpolant
          mt_vector<mt_real> rbf_k(_dpn, 0.);
          for(mt_int j=0; j < _NQ; j++)
          {
            const mt_real r_xj = (pt - _pts[_nbh_idx[k][j]]).length2();
            const mt_real rbf_val = eval_rbf(r_xj, _rbf_scal[k]);
            for(mt_int s=0; s < _dpn; s++) {
              rbf_k[s] += rbf_val * _rbf_coeffs[k][_dpn*j+s];
            }
          }
          for(mt_int s=0; s < _dpn; s++)
            odat[_dpn*nidx+s] += wk * rbf_k[s];
        }
        if(wsum > 0.)
          for(mt_int s=0; s < _dpn; s++)
            odat[_dpn*nidx+s] /= wsum;
      }
      #ifdef OPENMP
      #pragma omp critical
      #endif
      if(verbose) prg->next();
    }
  }
  if(verbose) {
    prg->finish();
    delete prg;
  }
}

void PointCloudInterpolator::interpolate_shepard(const mt_vector<mt_real> & inter_xyz,
                                                 const mt_vector<mt_real> & idat,
                                                 const int dpn, const bool global,
                                                 mt_vector<mt_real>& odat, const bool verbose) const
{
  if(_roi.size() == 0) {
    fprintf(stderr, "ERROR in %s: Radius of Influence has not been calculated!", __func__);
    exit(1);
  }

  mt_int nnodes = inter_xyz.size() / 3;
  odat.assign(dpn * nnodes, 0.0);

  mt_vector<vec3r> inter_pts;
  array_to_points(inter_xyz, inter_pts);

  const mt_real max_roi = *std::max_element(_roi.begin(), _roi.end());

  PROGRESS<size_t>* prg = NULL;
  if(verbose)
    prg = new PROGRESS<size_t>((size_t)nnodes, "Shepard-Interpolation progress: ");


  #ifdef OPENMP
  #pragma omp parallel
  #endif
  {
    mt_vector<mixed_tuple<mt_real, int> > int_vtx;

    #ifdef OPENMP
    #pragma omp for schedule(guided)
    #endif
    for(mt_int nidx = 0; nidx < nnodes; nidx++)
    {
      vec3r pt = inter_pts[nidx];
      mt_real sisum = 0.;
      bool hitpoint = false;

      if(global) {
        size_t cnt = 0;
        while(!hitpoint && cnt < _pts.size())
        {
          mt_real rk = (pt-_pts[cnt]).length();
          if(rk < 1e-8 * 0.5 * (pt.length() + _pts[cnt].length()))
          {
            hitpoint = true;
            for(mt_int s=0; s < dpn; s++)
              odat[dpn*nidx + s] = idat[dpn * cnt + s];
            sisum = 0.;
          } else {
            mt_real sk = 1. / (rk * rk * rk * rk);
            sisum += sk;
            for(mt_int s=0; s < dpn; s++)
              odat[dpn*nidx+s] += sk*idat[dpn*cnt+s];
          }
          cnt++;
        }
        if(sisum > 0)
          for(mt_int s=0; s < dpn; s++)
            odat[dpn*nidx+s] /= sisum;
      }
      else {
        _tree->vertices_in_sphere(pt, max_roi * 1.5, int_vtx);

        const mt_real rp = (int_vtx.size() > 0) ? int_vtx[int_vtx.size()-1].v1 : 0.;
        for(auto it = int_vtx.begin(); it != int_vtx.end(); ++it)
        {
          int k = it->v2;

          mt_real rk = (pt-_pts[k]).length();
          if(rk < 1e-8 * 0.5 * (pt.length() + _pts[k].length()))
          {
            for(mt_int s=0; s < dpn; s++)
              odat[dpn*nidx + s] = idat[dpn * k + s];
            continue;
          }
          mt_real si;
          if(rk <= rp /3)
            si = 1. / rk;
          else if(rk > rp / 3 && rk <= rp)
            si = (27. / (4. * rp)) *(rk/rp-1.)*(rk/rp-1.);
          else
            si = 0.;
          sisum += si*si;

          for(mt_int s=0; s < dpn; s++)
            odat[dpn*nidx+s] += si*si*idat[dpn*k+s];
        }
        if(sisum > 0)
          for(mt_int s=0; s < dpn; s++)
            odat[dpn*nidx+s] /= sisum;
      }

      #ifdef OPENMP
      #pragma omp critical
      #endif
      if(verbose)
        prg->next();
    }
  }
  if(verbose) {
    prg->finish();
    delete prg;
  }
}
