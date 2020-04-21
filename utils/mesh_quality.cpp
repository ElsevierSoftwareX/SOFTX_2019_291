/**
* @file mesh_quality.cpp
* @brief Functions for computing mesh quality.
* @author Aurel Neic
* @version
* @date 2017-08-16
*/

#include <assert.h>
#include "mt_utils_base.h"
#include "mesh_quality.h"
#include "topology_utils.h"
#include "geometry_utils.h"
#include "io_utils.h"
#include "lookup_table.hpp"

mt_real tet_qmetric_volume(const mt_int* tet, const mt_real* xyz)
{
  mt_real Qw, avg_edge = 0, tetVol;
  const mt_int v1 = tet[0], v2 = tet[1], v3 = tet[2], v4 = tet[3];

  vec3r p0(xyz + 3 * v1);
  vec3r p1(xyz + 3 * v2);
  vec3r p2(xyz + 3 * v3);
  vec3r p3(xyz + 3 * v4);

  tetVol = signed_tet_volume(p0, p1, p2, p3);
  //Zero Volume means bad quality
  if(tetVol <= 0.0)
    return 1.0;

  tetVol *= 6.0; //To actually get the jacobian determinant

  avg_edge += (p1-p0).length2();
  avg_edge += (p2-p0).length2();
  avg_edge += (p3-p0).length2();
  avg_edge += (p2-p1).length2();
  avg_edge += (p3-p1).length2();
  avg_edge += (p3-p2).length2();
  avg_edge /= 6.0;
  avg_edge = sqrt(avg_edge);
  Qw = sqrt(2.0) * tetVol / (avg_edge * avg_edge * avg_edge);

  return (1. - Qw);
}

mt_real tri_qmetric_surface(const mt_int* tet, const mt_real* xyz)
{
  mt_real Qw, avg_edge = 0;
  mt_int v1 = tet[0], v2 = tet[1], v3 = tet[2];

  vec3r p1(xyz + v1*3), p2(xyz + v2*3), p3(xyz + v3*3);
  vec3r e12 = p2 - p1, e13 = p3 - p1, e23 = p3 - p2;

  const mt_real surf = tri_surf(p1, p2, p3);
  if(surf <= 0.0)
    return 1.0;

  avg_edge += e12.length2();
  avg_edge += e13.length2();
  avg_edge += e23.length2();
  avg_edge /= 3;

  // Area of an even-sided triangle: a*a*sqrt(3)/4
  const mt_real avg_surf = avg_edge * 0.4330127;

  Qw = surf / avg_surf;

  return (1. - Qw);
}

mt_real tet_qmetric_shapeDist(const mt_int* tet, const mt_real* xyz)
{
  mt_int v0 = tet[0], v1 = tet[1], v2 = tet[2], v3 = tet[3];

  #ifdef OPENMP
  double x00, x01, x02;
  double x10, x11, x12;
  double x20, x21, x22;
  double x30, x31, x32;

  #pragma omp atomic read
  x00 = xyz[v0*3+0];
  #pragma omp atomic read
  x01 = xyz[v0*3+1];
  #pragma omp atomic read
  x02 = xyz[v0*3+2];
  #pragma omp atomic read
  x10 = xyz[v1*3+0];
  #pragma omp atomic read
  x11 = xyz[v1*3+1];
  #pragma omp atomic read
  x12 = xyz[v1*3+2];
  #pragma omp atomic read
  x20 = xyz[v2*3+0];
  #pragma omp atomic read
  x21 = xyz[v2*3+1];
  #pragma omp atomic read
  x22 = xyz[v2*3+2];
  #pragma omp atomic read
  x30 = xyz[v3*3+0];
  #pragma omp atomic read
  x31 = xyz[v3*3+1];
  #pragma omp atomic read
  x32 = xyz[v3*3+2];
  #else
  double x00 = xyz[v0*3+0], x01 = xyz[v0*3+1], x02 = xyz[v0*3+2];
  double x10 = xyz[v1*3+0], x11 = xyz[v1*3+1], x12 = xyz[v1*3+2];
  double x20 = xyz[v2*3+0], x21 = xyz[v2*3+1], x22 = xyz[v2*3+2];
  double x30 = xyz[v3*3+0], x31 = xyz[v3*3+1], x32 = xyz[v3*3+2];
  #endif

  double Frob = (3*POW2(x00) + 3*POW2(x01) + 3*POW2(x02) + 3*POW2(x10) +
                 3*POW2(x11) + 3*POW2(x12) - 2*x10*x20 + 3*POW2(x20) -
                 2*x11*x21 + 3*POW2(x21) - 2*x12*x22 + 3*POW2(x22) -
                 2*(x10 + x20)*x30 + 3*POW2(x30) - 2*x00*(x10 + x20 + x30) -
                 2*(x11 + x21)*x31 + 3*POW2(x31) - 2*x01*(x11 + x21 + x31) -
                 2*(x12 + x22)*x32 + 3*POW2(x32) - 2*x02*(x12 + x22 + x32)) / 2.;
  double det = 1.414213562*(x00*x12*x21 - x00*x11*x22 - x12*x21*x30 + x11*x22*x30 -
               x00*x12*x31 + x12*x20*x31 + x00*x22*x31 - x10*x22*x31 +
               x02*(-(x10*x21) + x11*(x20 - x30) + x21*x30 + x10*x31 - x20*x31) +
               x00*x11*x32 - x11*x20*x32 - x00*x21*x32 + x10*x21*x32 +
               x01*(-(x12*x20) + x10*x22 + x12*x30 - x22*x30 - x10*x32 + x20*x32));

  double eta = Frob / ( 3. * pow(fabs(det), 0.666666667));
  return (1. - 1./eta);
}


mt_real tet_qmetric_minimalSinus(const mt_int* tet, const mt_real* xyz)
{
  mt_real point[4][3];
  mt_real edgelength[3][4];
  mt_real facenormal[4][3]; // normal of tetface
  mt_real facearea2[4]; // squared area of tetface
  mt_real tetVol, sine2, minsine2;
  mt_real dx, dy, dz, dot;
  mt_int i,j,k,l;
  // FILL in the points of the tet
  for(i = 0; i < 4; i++)
    for(j=0; j < 3; j++)
      point[i][j] = xyz[3*tet[i] + j];

  vec3r p0(point[0]);
  vec3r p1(point[1]);
  vec3r p2(point[2]);
  vec3r p3(point[3]);
  tetVol = tet_volume(p0,p1,p2,p3);
  if(tetVol <= 0.0)
    return 1.0;

  tetVol *= 6.0;

  // for each vertex/face of the tetrahedron
  for (i = 0; i < 4; i++) {
      j = (i + 1) & 3;
      if ((i & 1) == 0) {
          k = (i + 3) & 3;
          l = (i + 2) & 3;
      } else {
          k = (i + 2) & 3;
          l = (i + 3) & 3;
      }
      //Compute face normals
      facenormal[i][0] =
          (point[k][1] - point[j][1]) * (point[l][2] - point[j][2]) -
          (point[k][2] - point[j][2]) * (point[l][1] - point[j][1]);
      facenormal[i][1] =
          (point[k][2] - point[j][2]) * (point[l][0] - point[j][0]) -
          (point[k][0] - point[j][0]) * (point[l][2] - point[j][2]);
      facenormal[i][2] =
          (point[k][0] - point[j][0]) * (point[l][1] - point[j][1]) -
          (point[k][1] - point[j][1]) * (point[l][0] - point[j][0]);

      // compute (2 *area)^2 for this face
      facearea2[i] = facenormal[i][0] * facenormal[i][0] +
          facenormal[i][1] * facenormal[i][1] +
          facenormal[i][2] * facenormal[i][2];

      // compute edge lengths (squared)
      for (j = i + 1; j < 4; j++) {
          dx = point[i][0] - point[j][0];
          dy = point[i][1] - point[j][1];
          dz = point[i][2] - point[j][2];
          edgelength[i][j] = dx * dx + dy * dy + dz * dz;
      }
  }

  minsine2 = 1.0e100;     /* start with absurdly big value for sine */
  // for each edge in the tetrahedron
  for (i = 0; i < 3; i++) {
      for (j = i + 1; j < 4; j++) {
          k = (i > 0) ? 0 : (j > 1) ? 1 : 2;
          l = 6 - i - j - k;
          // compute the expression for minimum sine, squared, over 4
          // The reason it's over 4 is because the area values we have
          // are actually twice the area squared
          // if either face area is zero, the sine is zero
          if (facearea2[k] > 0 && facearea2[l] > 0)
          {
              sine2 = edgelength[i][j] / (facearea2[k] * facearea2[l]);
          }
          else
          {
              #ifdef DEBUG
              fprintf(stderr, "WARNING in %s: Encountered zero-area face.\n", __func__);
              fprintf(stderr, "Suspicious element is (%ld %ld %ld %ld):\n", (mt_int) tet[0], (mt_int) tet[1], (mt_int) tet[2], (mt_int) tet[3]);
              #endif
              sine2 = 0.0;
          }

          // check whether this angle is obtuse
          dot = facenormal[k][0] * facenormal[l][0]
                + facenormal[k][1] * facenormal[l][1]
                + facenormal[k][2] * facenormal[l][2];
          if (dot > 0)
          {
              /* if so, warp it down */
              sine2 = 0.75 * sqrt(sine2);
              sine2 *= sine2;
          }

          /* update minimum sine */
          if (sine2 < minsine2)
          {
              minsine2 = sine2;
          }
      }
  }

  mt_real minsine = sqrt(minsine2) * tetVol;
  //Best sinus for equilateral tet is Sin(70.53°) = 2 * sqrt(2) / 3
  mt_real best_value = 0.9428090415820634;

  const mt_real dist = minsine < best_value ? fabs(minsine - best_value) / best_value : fabs(minsine - best_value);

  return dist;
}


mt_real tet_qmetric_weightedShapeDist(const mt_int* tet, const mt_real* xyz)
{
  mt_int v0 = tet[0], v1 = tet[1], v2 = tet[2], v3 = tet[3];

  #ifdef OPENMP
  double x00, x01, x02;
  double x10, x11, x12;
  double x20, x21, x22;
  double x30, x31, x32;

  #pragma omp atomic read
  x00 = xyz[v0*3+0];
  #pragma omp atomic read
  x01 = xyz[v0*3+1];
  #pragma omp atomic read
  x02 = xyz[v0*3+2];
  #pragma omp atomic read
  x10 = xyz[v1*3+0];
  #pragma omp atomic read
  x11 = xyz[v1*3+1];
  #pragma omp atomic read
  x12 = xyz[v1*3+2];
  #pragma omp atomic read
  x20 = xyz[v2*3+0];
  #pragma omp atomic read
  x21 = xyz[v2*3+1];
  #pragma omp atomic read
  x22 = xyz[v2*3+2];
  #pragma omp atomic read
  x30 = xyz[v3*3+0];
  #pragma omp atomic read
  x31 = xyz[v3*3+1];
  #pragma omp atomic read
  x32 = xyz[v3*3+2];
  #else
  double x00 = xyz[v0*3+0], x01 = xyz[v0*3+1], x02 = xyz[v0*3+2];
  double x10 = xyz[v1*3+0], x11 = xyz[v1*3+1], x12 = xyz[v1*3+2];
  double x20 = xyz[v2*3+0], x21 = xyz[v2*3+1], x22 = xyz[v2*3+2];
  double x30 = xyz[v3*3+0], x31 = xyz[v3*3+1], x32 = xyz[v3*3+2];
  #endif

  double Frob1 = sqrt((6*POW2(x00 - x10) + 6*POW2(x01 - x11) + 6*POW2(x02 - x12) +
                       2*POW2(x00 + x10 - 2*x20) + 2*POW2(x01 + x11 - 2*x21) +
                       2*POW2(x02 + x12 - 2*x22) + POW2(x00 + x10 + x20 - 3*x30) +
                       POW2(x01 + x11 + x21 - 3*x31) + POW2(x02 + x12 + x22 - 3*x32))/6.);
  double Frob2sq = (POW2(x00)*POW2(x11) + POW2(x00)*POW2(x12) - x00*POW2(x11)*x20 -
     x00*POW2(x12)*x20 + POW2(x11)*POW2(x20) + POW2(x12)*POW2(x20) -
     POW2(x00)*x11*x21 + x00*x10*x11*x21 + x00*x11*x20*x21 - 2*x10*x11*x20*x21 +
     POW2(x00)*POW2(x21) - x00*x10*POW2(x21) + POW2(x10)*POW2(x21) +
     POW2(x12)*POW2(x21) - POW2(x00)*x12*x22 + x00*x10*x12*x22 +
     x00*x12*x20*x22 - 2*x10*x12*x20*x22 - 2*x11*x12*x21*x22 +
     POW2(x00)*POW2(x22) - x00*x10*POW2(x22) + POW2(x10)*POW2(x22) +
     POW2(x11)*POW2(x22) - x00*POW2(x11)*x30 - x00*POW2(x12)*x30 -
     POW2(x11)*x20*x30 - POW2(x12)*x20*x30 + x10*x11*x21*x30 +
     x11*x20*x21*x30 - x00*POW2(x21)*x30 - x10*POW2(x21)*x30 +
     x10*x12*x22*x30 + x12*x20*x22*x30 - x00*POW2(x22)*x30 -
     x10*POW2(x22)*x30 + POW2(x11)*POW2(x30) + POW2(x12)*POW2(x30) -
     x11*x21*POW2(x30) + POW2(x21)*POW2(x30) - x12*x22*POW2(x30) +
     POW2(x22)*POW2(x30) - POW2(x00)*x11*x31 + x00*x10*x11*x31 +
     x10*x11*x20*x31 - x11*POW2(x20)*x31 - POW2(x00)*x21*x31 -
     POW2(x10)*x21*x31 - POW2(x12)*x21*x31 + x00*x20*x21*x31 +
     x10*x20*x21*x31 + x11*x12*x22*x31 + x12*x21*x22*x31 - x11*POW2(x22)*x31 +
     x00*x11*x30*x31 - 2*x10*x11*x30*x31 + x11*x20*x30*x31 + x00*x21*x30*x31 +
     x10*x21*x30*x31 - 2*x20*x21*x30*x31 + POW2(x00)*POW2(x31) -
     x00*x10*POW2(x31) + POW2(x10)*POW2(x31) + POW2(x12)*POW2(x31) -
     x00*x20*POW2(x31) - x10*x20*POW2(x31) + POW2(x20)*POW2(x31) -
     x12*x22*POW2(x31) + POW2(x22)*POW2(x31) +
     POW2(x02)*(POW2(x10) + POW2(x11) + POW2(x20) + POW2(x21) -
        x20*x30 + POW2(x30) - x10*(x20 + x30) - x21*x31 + POW2(x31) -
        x11*(x21 + x31)) + (-(x12*POW2(x20)) + x11*x12*x21 - x12*POW2(x21) -
        POW2(x10)*x22 - POW2(x11)*x22 + x11*x21*x22 -
        POW2(x00)*(x12 + x22) + x10*x12*(x20 - 2*x30) + x12*x20*x30 -
        2*x20*x22*x30 + x10*x22*(x20 + x30) +
        x00*(x10*x12 + x20*x22 + (x12 + x22)*x30) - 2*x11*x12*x31 + x12*x21*x31 +
        x11*x22*x31 - 2*x21*x22*x31)*x32 +
     (POW2(x00) + POW2(x10) + POW2(x11) - x10*x20 + POW2(x20) -
        x00*(x10 + x20) - x11*x21 + POW2(x21))*POW2(x32) +
     POW2(x01)*(POW2(x10) + POW2(x12) + POW2(x20) + POW2(x22) -
        x20*x30 + POW2(x30) - x10*(x20 + x30) - x22*x32 + POW2(x32) -
        x12*(x22 + x32)) + x02*(-(x12*POW2(x20)) + x11*x12*x21 -
        x12*POW2(x21) - POW2(x11)*x22 + x11*x21*x22 + x20*x22*x30 -
        x12*POW2(x30) - x22*POW2(x30) + x11*x12*x31 + x21*x22*x31 -
        x12*POW2(x31) - x22*POW2(x31) -
        (POW2(x11) + POW2(x20) + POW2(x21) - x20*x30 - (x11 + x21)*x31)*
         x32 - POW2(x10)*(x22 + x32) +
        x10*(x20*x22 + x12*(x20 + x30) + x30*x32) +
        x00*(-2*x20*x22 + x22*x30 + x12*(x20 + x30) + x20*x32 - 2*x30*x32 +
           x10*(-2*x12 + x22 + x32))) +
     x01*(x10*x11*x20 - x11*POW2(x20) - POW2(x10)*x21 - POW2(x12)*x21 +
        x10*x20*x21 + x11*x12*x22 + x12*x21*x22 - x11*POW2(x22) + x10*x11*x30 +
        x20*x21*x30 - x11*POW2(x30) - x21*POW2(x30) - POW2(x10)*x31 -
        POW2(x12)*x31 - POW2(x20)*x31 - POW2(x22)*x31 + x10*x30*x31 +
        x20*x30*x31 + x00*(-2*x20*x21 + x21*x30 + x11*(x20 + x30) + x20*x31 -
           2*x30*x31 + x10*(-2*x11 + x21 + x31)) +
        (x11*x12 + x21*x22 + (x12 + x22)*x31)*x32 - (x11 + x21)*POW2(x32) +
        x02*(-2*x21*x22 + x22*x31 + x12*(x21 + x31) + x21*x32 - 2*x31*x32 +
           x11*(-2*x12 + x22 + x32))))/
   POW2(x00*x12*x21 - x00*x11*x22 - x12*x21*x30 + x11*x22*x30 - x00*x12*x31 +
     x12*x20*x31 + x00*x22*x31 - x10*x22*x31 +
     x02*(-(x10*x21) + x11*(x20 - x30) + x21*x30 + x10*x31 - x20*x31) +
     x00*x11*x32 - x11*x20*x32 - x00*x21*x32 + x10*x21*x32 +
     x01*(-(x12*x20) + x10*x22 + x12*x30 - x22*x30 - x10*x32 + x20*x32));
  double Frob2 = sqrt(Frob2sq);
  double kappa = Frob1 * Frob2;
  return fabs(1. - 3. / kappa);
}

void mesh_volumes(mt_meshdata & mesh,
                  mt_vector<mt_real> & xyz,
                  mt_vector<mt_real> & volumes)
{
  size_t numelem = mesh.e2n_cnt.size();
  volumes.resize(numelem);

  if(mesh.e2n_dsp.size() == 0)
    bucket_sort_offset(mesh.e2n_cnt, mesh.e2n_dsp);

  #ifdef OPENMP
  #pragma omp parallel for schedule(guided)
  #endif
  for(size_t i=0; i<numelem; i++) {
    const mt_int* con = mesh.e2n_con.data() + mesh.e2n_dsp[i];
    volumes[i] = volume(mesh.etype[i], con, xyz.data());
  }
}

void mesh_quality(const mt_meshdata & mesh,
                     mt_vector<mt_real> & qual)
{
  size_t numelem = mesh.e2n_cnt.size();
  qual.resize(numelem);
  const mt_real* xyz = mesh.xyz.data();

  #ifdef OPENMP
  #pragma omp parallel for schedule(dynamic, 10)
  #endif
  for(size_t i=0; i<numelem; i++) {
    const mt_int* con = mesh.e2n_con.data() + mesh.e2n_dsp[i];
    qual[i] = element_quality(mesh.etype[i], con, xyz);
  }
}

void mesh_quality(const mt_meshdata & mesh, mt_real & min, mt_real & max, mt_real & avrg)
{
  mt_vector<mt_real> qual;
  mesh_quality(mesh, qual);
  min = qual[0], max = qual[0], avrg = 0.0;

  for(const mt_real & q : qual)
  {
    avrg += q;
    if(min > q) min = q;
    if(max < q) max = q;
  }
  avrg /= mt_real(qual.size());
}

void mesh_quality(const mt_meshdata & mesh, mt_real & min, mt_real & max, mt_real & avrg,
                  const mt_real & thr, mt_real & percentage)
{
  mt_vector<mt_real> qual;
  mesh_quality(mesh, qual);
  min = qual[0], max = qual[0], avrg = 0.0;
  mt_int bad_elem_cnt = 0;
  for(const mt_real & q : qual)
  {
    avrg += q;
    if(min > q) min = q;
    if(max < q) max = q;
    if(q >= thr)
      bad_elem_cnt++;
  }
  avrg /= mt_real(qual.size());
  percentage = 100. * ((mt_real) bad_elem_cnt / (mt_real)(qual.size()));
}

mt_real nbhd_quality_average(const mt_meshdata & mesh,
                          const mt_int nidx)
{
  const mt_real* xyz = mesh.xyz.data();
  mt_real ret = 0;

  mt_int start = mesh.n2e_dsp[nidx], stop = start + mesh.n2e_cnt[nidx];

  for(mt_int i=start; i<stop; i++) {
    mt_int eidx = mesh.n2e_con[i];
    const mt_int* con = mesh.e2n_con.data() + mesh.e2n_dsp[eidx];
    ret += element_quality(mesh.etype[eidx], con, xyz);
  }

  ret /= (stop - start);
  return ret;
}

mt_real nbhd_quality_max(const mt_meshdata & mesh,
                         const mt_int nidx)
{
  const mt_real* xyz = mesh.xyz.data();
  mt_real ret = 0;

  mt_int start = mesh.n2e_dsp[nidx], stop = start + mesh.n2e_cnt[nidx];

  for(mt_int i=start; i<stop; i++) {
    mt_int eidx = mesh.n2e_con[i];
    const mt_int* con = mesh.e2n_con.data() + mesh.e2n_dsp[eidx];

    mt_real m = element_quality(mesh.etype[eidx], con, xyz);

    if(m != m) {
      fprintf(stderr, "%s error: Computed NaN element quality!\n", __func__);
      return m;
    }

    if(ret < m) ret = m;
  }
  return ret;
}

void nbhd_quality_max(const mt_meshdata & mesh,
                      const mt_int nidx,
                      mt_real & max_qual,
                      mt_int & max_eidx)
{
  const mt_real* xyz = mesh.xyz.data();
  max_qual = -1.0;
  max_eidx = -1;

  mt_int start = mesh.n2e_dsp[nidx], stop = start + mesh.n2e_cnt[nidx];
  for(mt_int i=start; i<stop; i++) {
    mt_int eidx = mesh.n2e_con[i];
    const mt_int* con = mesh.e2n_con.data() + mesh.e2n_dsp[eidx];

    mt_real m = element_quality(mesh.etype[eidx], con, xyz);

    if(m != m) {
      max_qual = -1.0;
      max_eidx = -1;
      fprintf(stderr, "%s error: Computed NaN element quality!\n", __func__);
      return;
    }

    if(max_qual < m) {
      max_qual = m;
      max_eidx = eidx;
    }
  }
}

void mesh_min_angles(const mt_meshdata & mesh,
                       mt_vector<mt_real> & minangles)
{
  size_t numelem = mesh.e2n_cnt.size();
  minangles.resize(numelem);
  const double normal_coeff = 180. * .3183098861837906715377675267450287;

  #ifdef OPENMP
  #pragma omp parallel for schedule(dynamic, 10)
  #endif
  for(size_t i=0; i<numelem; i++) {
    if(mesh.etype[i] == Tetra) {
      const mt_int*  con    = mesh.e2n_con.data() + mesh.e2n_dsp[i];
      const mt_real* xyz    = mesh.xyz.data();
      mt_int v0 = con[0], v1 = con[1], v2 = con[2], v3 = con[3];

      #ifdef OPENMP
      double x00, x01, x02;
      double x10, x11, x12;
      double x20, x21, x22;
      double x30, x31, x32;

      #pragma omp atomic read
      x00 = xyz[v0*3+0];
      #pragma omp atomic read
      x01 = xyz[v0*3+1];
      #pragma omp atomic read
      x02 = xyz[v0*3+2];
      #pragma omp atomic read
      x10 = xyz[v1*3+0];
      #pragma omp atomic read
      x11 = xyz[v1*3+1];
      #pragma omp atomic read
      x12 = xyz[v1*3+2];
      #pragma omp atomic read
      x20 = xyz[v2*3+0];
      #pragma omp atomic read
      x21 = xyz[v2*3+1];
      #pragma omp atomic read
      x22 = xyz[v2*3+2];
      #pragma omp atomic read
      x30 = xyz[v3*3+0];
      #pragma omp atomic read
      x31 = xyz[v3*3+1];
      #pragma omp atomic read
      x32 = xyz[v3*3+2];
      #else
      double x00 = xyz[v0*3+0], x01 = xyz[v0*3+1], x02 = xyz[v0*3+2];
      double x10 = xyz[v1*3+0], x11 = xyz[v1*3+1], x12 = xyz[v1*3+2];
      double x20 = xyz[v2*3+0], x21 = xyz[v2*3+1], x22 = xyz[v2*3+2];
      double x30 = xyz[v3*3+0], x31 = xyz[v3*3+1], x32 = xyz[v3*3+2];
      #endif
      double nabc = sqrt(POW2(-(x01*x10) + x00*x11 + x01*x20 - x11*x20 - x00*x21 + x10*x21) +
      POW2(x02*x10 - x00*x12 - x02*x20 + x12*x20 + x00*x22 - x10*x22) +
      POW2(-(x02*x11) + x01*x12 + x02*x21 - x12*x21 - x01*x22 + x11*x22));
      double nabd = sqrt(POW2(-(x01*x10) + x00*x11 + x01*x30 - x11*x30 - x00*x31 + x10*x31) +
      POW2(x02*x10 - x00*x12 - x02*x30 + x12*x30 + x00*x32 - x10*x32) +
      POW2(-(x02*x11) + x01*x12 + x02*x31 - x12*x31 - x01*x32 + x11*x32));
      double nacd = sqrt(POW2(-(x01*x20) + x00*x21 + x01*x30 - x21*x30 - x00*x31 + x20*x31) +
      POW2(x02*x20 - x00*x22 - x02*x30 + x22*x30 + x00*x32 - x20*x32) +
      POW2(-(x02*x21) + x01*x22 + x02*x31 - x22*x31 - x01*x32 + x21*x32));
      double nbcd = sqrt(POW2(-(x11*x20) + x10*x21 + x11*x30 - x21*x30 - x10*x31 + x20*x31) +
      POW2(x12*x20 - x10*x22 - x12*x30 + x22*x30 + x10*x32 - x20*x32) +
      POW2(-(x12*x21) + x11*x22 + x12*x31 - x22*x31 - x11*x32 + x21*x32));

      double sprod_abc_abd = (x01*(x10 - x20) + x11*x20 - x10*x21 + x00*(-x11 + x21))*
       (x01*(x10 - x30) + x11*x30 - x10*x31 + x00*(-x11 + x31)) +
      (x02*(x10 - x20) + x12*x20 - x10*x22 + x00*(-x12 + x22))*
       (x02*(x10 - x30) + x12*x30 - x10*x32 + x00*(-x12 + x32)) +
      (x02*(x11 - x21) + x12*x21 - x11*x22 + x01*(-x12 + x22))*
       (x02*(x11 - x31) + x12*x31 - x11*x32 + x01*(-x12 + x32));

      double sprod_abc_acd = (x01*(x10 - x20) + x11*x20 - x10*x21 + x00*(-x11 + x21))*
       (x01*(x20 - x30) + x21*x30 - x20*x31 + x00*(-x21 + x31)) +
      (x02*(x10 - x20) + x12*x20 - x10*x22 + x00*(-x12 + x22))*
       (x02*(x20 - x30) + x22*x30 - x20*x32 + x00*(-x22 + x32)) +
      (x02*(x11 - x21) + x12*x21 - x11*x22 + x01*(-x12 + x22))*
       (x02*(x21 - x31) + x22*x31 - x21*x32 + x01*(-x22 + x32));

      double sprod_abc_bcd = (x01*(x10 - x20) + x11*x20 - x10*x21 + x00*(-x11 + x21))*
       (x11*(x20 - x30) + x21*x30 - x20*x31 + x10*(-x21 + x31)) +
      (x02*(x10 - x20) + x12*x20 - x10*x22 + x00*(-x12 + x22))*
       (x12*(x20 - x30) + x22*x30 - x20*x32 + x10*(-x22 + x32)) +
      (x02*(x11 - x21) + x12*x21 - x11*x22 + x01*(-x12 + x22))*
       (x12*(x21 - x31) + x22*x31 - x21*x32 + x11*(-x22 + x32));

      double sprod_abd_acd = (-(x11*x30) + x01*(-x10 + x30) + x00*(x11 - x31) + x10*x31)*
       (-(x21*x30) + x01*(-x20 + x30) + x00*(x21 - x31) + x20*x31) +
      (-(x12*x31) + x02*(-x11 + x31) + x01*(x12 - x32) + x11*x32)*
       (-(x22*x31) + x02*(-x21 + x31) + x01*(x22 - x32) + x21*x32) +
      (x02*(x10 - x30) + x12*x30 - x10*x32 + x00*(-x12 + x32))*
       (x02*(x20 - x30) + x22*x30 - x20*x32 + x00*(-x22 + x32));

      double sprod_abd_bcd = (x01*(x10 - x30) + x11*x30 - x10*x31 + x00*(-x11 + x31))*
       (x11*(x20 - x30) + x21*x30 - x20*x31 + x10*(-x21 + x31)) +
      (x02*(x10 - x30) + x12*x30 - x10*x32 + x00*(-x12 + x32))*
       (x12*(x20 - x30) + x22*x30 - x20*x32 + x10*(-x22 + x32)) +
      (x02*(x11 - x31) + x12*x31 - x11*x32 + x01*(-x12 + x32))*
       (x12*(x21 - x31) + x22*x31 - x21*x32 + x11*(-x22 + x32));

      double sprod_acd_bcd = (x01*(x20 - x30) + x21*x30 - x20*x31 + x00*(-x21 + x31))*
       (x11*(x20 - x30) + x21*x30 - x20*x31 + x10*(-x21 + x31)) +
      (x02*(x20 - x30) + x22*x30 - x20*x32 + x00*(-x22 + x32))*
       (x12*(x20 - x30) + x22*x30 - x20*x32 + x10*(-x22 + x32)) +
      (x02*(x21 - x31) + x22*x31 - x21*x32 + x01*(-x22 + x32))*
       (x12*(x21 - x31) + x22*x31 - x21*x32 + x11*(-x22 + x32));


      std::vector<double> angles(6);
      angles[0] = acos(sprod_abc_abd / (nabc * nabd));
      angles[1] = acos(sprod_abc_acd / (nabc * nacd));
      angles[2] = acos(sprod_abc_bcd / (nabc * nbcd));
      angles[3] = acos(sprod_abd_acd / (nabd * nacd));
      angles[4] = acos(sprod_abd_bcd / (nabd * nbcd));
      angles[5] = acos(sprod_acd_bcd / (nacd * nbcd));

      auto result = std::min_element(angles.begin(), angles.end());
      const double optimal_angle = 70.52877936550931;
      double alphamin = normal_coeff * (*result);
      minangles[i] = std::min(1.0, fabs(optimal_angle - alphamin) / optimal_angle);
    }
    else {
      minangles[i] = INFINITY;
    }
  }
}

mt_real distance_to_centroid(const mt_meshdata & mesh,
                            const mt_meshdata & manifold,
                            const mt_vector<bool> & onManifold,
                            const mt_int vtx)
{
  assert(mesh.xyz.size());
  assert(mesh.n2n_cnt.size());
  assert(manifold.n2n_cnt.size());

  const mt_meshdata *m;
  if(onManifold[vtx]) m = &manifold;
  else                m = &mesh;

  mt_point<mt_real> rpt(mesh.xyz.data() + vtx*3);
  mt_point<mt_real> avrg;

  mt_int num_neigh = m->n2n_cnt[vtx];
  if(num_neigh > 1) {
    mt_int start = m->n2n_dsp[vtx], stop = start + num_neigh;
    for(mt_int i=start; i<stop; i++)
    {
      mt_int cvtx = m->n2n_con[i];
      avrg += mt_point<mt_real>(mesh.xyz.data() + cvtx*3);
    }
    avrg -= rpt;
    avrg /= mt_real(num_neigh - 1);
    avrg -= rpt;

    return avrg.length();
  }
  else
    return 0.0;
}

mt_real minimal_surf_normal_correlation(const mt_meshdata & surfmesh,
                                       const mt_vector<mt_real> & xyz,
                                       const mt_int vtx)
{
  mt_point<mt_real> nref;
  double min = 1.0;

  mt_int start = surfmesh.n2e_dsp[vtx], end = start + surfmesh.n2e_cnt[vtx];
  for(mt_int i=start; i<end; i++)
  {
    mt_int cele = surfmesh.n2e_con[i];
    mt_point<mt_real> p0(xyz.data() + surfmesh.e2n_con[cele*3+0]*3);
    mt_point<mt_real> p1(xyz.data() + surfmesh.e2n_con[cele*3+1]*3);
    mt_point<mt_real> p2(xyz.data() + surfmesh.e2n_con[cele*3+2]*3);
    mt_point<mt_real> e1 = p1 - p0, e2 = p2 - p0;
    mt_point<mt_real> n = e1.crossProd(e2);
    n.normalize();

    if(nref.length2() > 0) {
      // a refrerence normal has been already set. we can compute a current normal
      // and do the angle check
      mt_real pro = nref.scaProd(n);
      if(min > pro) min = pro;
    }
    else nref = n; // set reference
  }

  // we dont allow negative results, 0 is bad enough
  min = min < 0 ? 0 : min;

  return min;
}


mt_real nbhd_quality_functional(const mt_meshdata & mesh,
                                const mt_int nidx,
                                const mt_real thr)
{
  const mt_real* xyz = mesh.xyz.data();
  mt_real ret = 0;
  mt_int start = mesh.n2e_dsp[nidx], stop = start + mesh.n2e_cnt[nidx];

  for(mt_int i=start; i<stop; i++) {
    mt_int eidx = mesh.n2e_con[i];
    const mt_int* con = mesh.e2n_con.data() + mesh.e2n_dsp[eidx];
    const mt_real q = element_quality(mesh.etype[eidx], con, xyz);
    const mt_real m = q + (1. - thr);
    ret += POW10(m);
  }

  return ret;
}

mt_real nbhd_smoothness_functional(const mt_meshdata & surfmesh,
                                   const mt_vector<mt_real> & xyz,
                                   const mt_int nidx,
                                   const mt_real thr)
{
#if 0
  mt_real c = minimal_surf_normal_correlation(surfmesh, xyz, nidx);
  mt_real m =  (1.0 - c) + (1. - thr);
  mt_real ret = POW10(m);
#else
  mt_point<mt_real> nref;
  mt_real ret = 0.0;

  mt_int start = surfmesh.n2e_dsp[nidx], end = start + surfmesh.n2e_cnt[nidx];
  for(mt_int i=start; i<end; i++)
  {
    mt_int cele = surfmesh.n2e_con[i];
    mt_point<mt_real> p0(xyz.data() + surfmesh.e2n_con[cele*3+0]*3);
    mt_point<mt_real> p1(xyz.data() + surfmesh.e2n_con[cele*3+1]*3);
    mt_point<mt_real> p2(xyz.data() + surfmesh.e2n_con[cele*3+2]*3);
    mt_point<mt_real> e1 = p1 - p0, e2 = p2 - p0;
    mt_point<mt_real> n = e1.crossProd(e2);
    n.normalize();

    if(nref.length2() > 0) {
      mt_real s = nref.scaProd(n);
      s = s < 0.0 ? 0.0 : s;
      mt_real p = (1.0 - s) + (1.0 - thr);
      ret += POW10(p);
    }
    else nref = n;
  }
#endif

  return ret;
}

vec3r quality_gradient(mt_meshdata & mesh, const mt_int vidx, const mt_real thr, const mt_real delta)
{
  vec3r grad(0,0,0), dx(delta, 0, 0), dy(0, delta, 0), dz(0, 0, delta);
  mt_real* vp = mesh.xyz.data() + vidx*3;
  vec3r vert(vp);

  mt_real init_val = nbhd_quality_functional(mesh, vidx, thr);

  vert += dx; vert.set(vp);
  grad.x = nbhd_quality_functional(mesh, vidx, thr) - init_val;
  vert = vert - dx + dy; vert.set(vp);
  grad.y = nbhd_quality_functional(mesh, vidx, thr) - init_val;
  vert = vert - dy + dz; vert.set(vp);
  grad.z = nbhd_quality_functional(mesh, vidx, thr) - init_val;
  vert -= dz; vert.set(vp);

  return grad;
}

void improve_nodeset_gradient_method(mt_meshdata & mesh, const mt_meshdata & surf,
                                   const mt_mask & is_surf,
                                   const MT_USET<mt_int> & nodes, const mt_real thr,
                                   const mt_real relative_step)
{
  mt_real min_l, max_l, avg_l;
  mt_vector<mt_int> nod_vec; nod_vec.assign(nodes.begin(), nodes.end());
  size_t nnod = nod_vec.size();

  bool verbose = false;

  // #ifdef OPENMP
  // #pragma omp parallel for schedule(guided, 10)
  // #endif
  for(size_t i=0; i<nnod; i++)
  {
    mt_int nidx = nod_vec[i];
    mt_int eidx = mesh.n2e_con[mesh.n2e_dsp[nidx]];
    element_edges_stats(mesh.etype[eidx], mesh.e2n_con.data() + mesh.e2n_dsp[eidx],
                        mesh.xyz.data(), min_l, max_l, avg_l);

    mt_real  abs_step = min_l * relative_step;
    mt_real* vert_ptr = mesh.xyz.data() + nidx*3;
    vec3r    vert(vert_ptr);

    vec3r grad = unit_vector(quality_gradient(mesh, nidx, thr, 1e-4 * min_l));

    // if vertex is on surface, we remove the surface normal components from all connected
    // surface elements as to not damage the surface smoothness
    if(is_surf.count(nidx)) {
      mt_int start = surf.n2e_dsp[nidx], stop = start + surf.n2e_cnt[nidx];
      for(mt_int j=start; j<stop; j++) {
        mt_int seidx = surf.n2e_con[j];
        mt_int v0 = surf.e2n_con[surf.e2n_dsp[seidx]],
               v1 = surf.e2n_con[surf.e2n_dsp[seidx] + 1],
               v2 = surf.e2n_con[surf.e2n_dsp[seidx] + 2];

        vec3r e1 = vec3r(mesh.xyz.data()+v1*3) - vec3r(mesh.xyz.data()+v0*3);
        vec3r e2 = vec3r(mesh.xyz.data()+v2*3) - vec3r(mesh.xyz.data()+v0*3);
        vec3r n  = unit_vector(e1.crossProd(e2));
        grad -= n * n.scaProd(grad);
      }
    }

    // if the gradient is too small we skip the node
    if(grad.length() < 1e-4) continue;
    else                     grad.normalize();

    vec3r step = grad * (-abs_step);

    // improve the quality of the element neighbourhood by moving along the negative
    // element quality gradient direciton
    mt_real init_qual = nbhd_quality_functional(mesh, nidx, thr);

    if(verbose)
    printf("Vertex %d: initial quality %f, descent dir %f,%f,%f\n", int(nidx), float(init_qual),
           float(step.x), float(step.y), float(step.z));

    vec3r   new_pt    = vert + step; new_pt.set(vert_ptr);
    mt_real new_qual  = nbhd_quality_functional(mesh, nidx, thr);

    if(verbose)
    printf("Vertex %d: new quality %f, descent dir %f,%f,%f\n", int(nidx), float(new_qual),
           float(step.x), float(step.y), float(step.z));

    // do backtracking if we didnt improve the neighbourhood quality
    while(new_qual > init_qual + 1e-2) {
      step *= 0.8;
      new_pt = vert + step;
      new_pt.set(vert_ptr);
      new_qual = nbhd_quality_functional(mesh, nidx, thr);

      if(verbose)
      printf("Vertex %d: new quality %f, descent dir %f,%f,%f\n", int(nidx), float(new_qual),
             float(step.x), float(step.y), float(step.z));
    }
  }
}

void intersect_edges_with_faces(const mt_meshdata & mesh,
                                const mt_vector<mt_int> & loc2glob,
                                mt_vector<tuple<mt_int>> & intersec_edges,
                                mt_vector<mt_int>        & intersec_eidx)
{
  mt_mapping<mt_int> ele2edge, ele2face;
  MT_MAP<tuple<mt_int>, mt_int> edge_map;
  MT_MAP<triple<mt_int>, mt_int> face_map;

  compute_edges(mesh, ele2edge, edge_map);
  compute_faces(mesh, ele2face, face_map);

  size_t widx = 0;
  mt_vector<tuple<mt_int>> edges(edge_map.size());
  for(const auto & it : edge_map) {
    edges[widx++] = it.first;
  }

  ele2face.transpose();
  ele2face.setup_dsp();

  mt_vector<tri_elem> trivec(face_map.size());

  auto fit = face_map.begin();
  for(size_t i=0; i<trivec.size(); i++, ++fit) {
    const triple<mt_int> & tri = fit->first;
    trivec[i].v0 = vec3r(mesh.xyz.data() + tri.v1*3);
    trivec[i].v1 = vec3r(mesh.xyz.data() + tri.v2*3);
    trivec[i].v2 = vec3r(mesh.xyz.data() + tri.v3*3);
    mt_int face_idx = fit->second;
    mt_int elem_idx = ele2face.bwd_con[ele2face.bwd_dsp[face_idx]];
    trivec[i].eidx  = elem_idx;
  }

  kdtree tree(10);
  tree.build_tree(trivec);

  #ifdef OPENMP
  #pragma omp parallel for schedule(guided, 10)
  #endif
  for(size_t i=0; i<edges.size(); i++)
  {
    const tuple<mt_int> & edge = edges[i];
    tri_elem hit_ele;
    vec3r    hit_pos;

    vec3r p0(mesh.xyz.data() + edge.v1*3), p1(mesh.xyz.data() + edge.v2*3);
    ray r(p0, p1 - p0);

    bool did_hit = tree.closest_intersect(r, 0.01, 0.99, hit_ele, hit_pos);

    if(did_hit) {
#if 0
      fprintf(stderr, "Self intersection detected: edge (%ld, %ld) with elem %d at t = %f!\n",
              edge.v1, edge.v2, int(hit_ele.eidx), (hit_pos-p0).length() / (p1 - p0).length());
#endif

      #ifdef OPENMP
      #pragma omp critical
      #endif
      {
        intersec_edges.push_back(edge);
        intersec_eidx.push_back(loc2glob[hit_ele.eidx]);
      }
    }
  }
}

void check_self_intersection(const mt_meshdata & mesh,
                             mt_vector<tuple<mt_int>> & intersec_edges,
                             mt_vector<mt_int> & intersec_eidx)
{
  const size_t nele = mesh.e2n_cnt.size();
  mt_vector<mt_int> dsp(nele);
  bucket_sort_offset(mesh.e2n_cnt, dsp);

  bbox box;
  generate_bbox(mesh.xyz, box);
  bboxAxis longest_axis     = get_longest_axis(box);
  mt_real  longest_axis_len = get_longest_axis_length(box);

  const int     numblocks = 50;
  const mt_real blocksize = (longest_axis_len / numblocks) * 1.01;

  mt_meshdata cur_mesh_block;
  mt_vector<mt_int> sel_eidx;
  sel_eidx.reserve(nele / (numblocks * 0.5));

  PROGRESS<int> prg(numblocks, "Self-intersection check: ");


  for(int bidx = 0; bidx < numblocks; bidx++)
  {
    sel_eidx.resize(0);

    #ifdef OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for(size_t eidx = 0; eidx < nele; eidx++) {
      const mt_int * nod = mesh.e2n_con.data() + dsp[eidx];
      vec3r ctr = barycenter(mesh.e2n_cnt[eidx], nod, mesh.xyz.data());

      mt_real coord = -1.0;
      switch(longest_axis) {
        default: break;
        case X: coord = ctr.x - box.bounds[0].x; break;
        case Y: coord = ctr.y - box.bounds[0].y; break;
        case Z: coord = ctr.z - box.bounds[0].z; break;
      }

      if(int(coord / blocksize) == bidx) {
        #ifdef OPENMP
        #pragma omp critical
        #endif
        {
          sel_eidx.push_back(eidx);
        }
      }
    }

    extract_mesh(sel_eidx, mesh, cur_mesh_block);

    if(cur_mesh_block.e2n_cnt.size()) {
      compute_full_mesh_connectivity(cur_mesh_block, false);
      intersect_edges_with_faces(cur_mesh_block, sel_eidx, intersec_edges, intersec_eidx);
    }

    prg.next();
  }
  prg.finish();
}


