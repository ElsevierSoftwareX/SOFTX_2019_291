/**
* @file geometry_utils.h
* @brief Utility functions for geometric operations.
* @author Aurel Neic
* @version
* @date 2018-10-28
*/

#ifndef _GEOMETRY_UTILS
#define _GEOMETRY_UTILS

#include "kdtree.h"
#include "dense_mat.hpp"

/**
* @brief Compute the angle between two lines
*	 based on https://people.eecs.berkeley.edu/~wkahan/Triangle.pdf
* @param A  Base for angle
* @param B  First endpoint of line
* @param C  Second endpoint of line
*
* @return The angle in radians
*/
mt_real angle_between_line(const vec3r & A, const vec3r & B, const vec3r & C);

/**
* @brief Compute the normal of a triangle
*
* @param v1  First triangle vertex.
* @param v2  Second triangle vertex.
* @param v3  Third triangle vertex.
* @param xyz Coords vector.
*
* @return The normal unit vector.
*/
vec3r triangle_normal(const mt_int v1, const mt_int v2, const mt_int v3,
                      const mt_vector<mt_real> & xyz);
vec3r triangle_normal(const vec3r & p0, const vec3r & p1, const vec3r & p2);

/**
* @brief Compute the centerpoint of a triangle
*
* @param v1  First triangle vertex.
* @param v2  Second triangle vertex.
* @param v3  Third triangle vertex.
* @param xyz Coords vector.
*
* @return The centerpoint coordinates.
*/
vec3r triangle_centerpoint(const mt_int v1, const mt_int v2, const mt_int v3,
                           const mt_vector<mt_real> & xyz);

/**
* @brief Compute best L2 fit of a sphere w.r.t a set of vertices
*
* @param [in]  vert   The set of vertices.
* @param [out] R      The sphere radius.
* @param [out] pos    The sphere center point
*/
void fit_sphere(const mt_vector<vec3r> & vert, mt_real & R, vec3r & pos);

/**
* @brief Compute the volume of a single element.
*
* @param type  Element type.
* @param con   Pointer to connectivity of the element.
* @param xyz   Pointer to start of all coords.
*
* @return Volume of element.
*/
mt_real volume(const elem_t type, const mt_int* con, const mt_real* xyz);

/// overload to compute volume on an entire mesh
inline mt_real volume(const mt_meshdata & mesh)
{
  mt_real ret = 0.0;
  const mt_int* con = mesh.e2n_con.data();

  for(size_t eidx = 0; eidx < mesh.e2n_cnt.size(); eidx++)
  {
    ret += volume(mesh.etype[eidx], con, mesh.xyz.data());
    con += mesh.e2n_cnt[eidx];
  }

  return ret;
}

/// Return triangle surface
mt_real tri_surf(vec3r & p0, vec3r & p1, vec3r & p2);
/// Wrapper of tri_surf for vec3r arrays
mt_real tri_surf(vec3r * points);
/// Return Quad surface
mt_real quad_surf(vec3r * points);

/**
* @brief Compute volume of a tetrahedra.
*
* @param points Array pointer holding the four tetra vertices.
*
* @return The tetra volume.
*/
mt_real tet_volume(mt_point<mt_real> * points);

/**
* @brief Compute volume of a tetrahedra.
*
* The four tetra vertices are passed individually.
*
* @return The tetra volume.
*/
mt_real tet_volume(const mt_point<mt_real> & p0,
                   const mt_point<mt_real> & p1,
                   const mt_point<mt_real> & p2,
                   const mt_point<mt_real> & p3);

/**
* @brief Compute signed volume of a tetrahedra.
*
* @param points Array pointer holding the four tetra vertices.
*
* @return The tetra volume. negative if tet is wrongly oriented
*/
mt_real signed_tet_volume(mt_point<mt_real> * points);

/**
* @brief Compute signed volume of a tetrahedra.
*
* The four tetra vertices are passed individually.
*
* @return The tetra volume. negative if tet is wrongly oriented
*/
mt_real signed_tet_volume(const mt_point<mt_real> & p0,
                   const mt_point<mt_real> & p1,
                   const mt_point<mt_real> & p2,
                   const mt_point<mt_real> & p3);


mt_real signed_tri_surf(const vec3r & p0, const vec3r & p1, const vec3r & p2);
mt_real signed_tri_surf(vec3r * points);

/**
* @brief Compute volume of a pyramid.
*
* @param points Array pointer holding the five pyramid vertices.
*
* @return The pyramid volume.
*/
mt_real pyr_volume(mt_point<mt_real> * points);

/**
* @brief Compute volume of a hexahedra.
*
* @param points Array pointer holding the eight hex vertices.
*
* @return The hex volume.
*/
mt_real hex_volume(mt_point<mt_real> * points);

/**
* @brief Compute volume of a prism.
*
* @param points Array pointer holding the six prism vertices.
*
* @return The prism volume.
*/
mt_real prism_volume(mt_point<mt_real> * pts);

/**
* @brief Rotate a points array w.r.t. x-, y- and z-axes.
*
* @param pts   The points.
* @param rotx  The angle (in degrees) we rotate around the x axis.
* @param roty  The angle (in degrees) we rotate around the y axis.
* @param rotz  The angle (in degrees) we rotate around the z axis.
*/
void rotate_points(mt_vector<vec3r> & pts,
                   const mt_real rotx,
                   const mt_real roty,
                   const mt_real rotz);


/**
* @brief Compute tetrahedral circumsphere.
*
* @param[in] p0  First tet point.
* @param[in] p1  Second tet point.
* @param[in] p2  Third tet point.
* @param[in] p3  Fourth tet point.
* @param[out] ctr Midpoint of tetrahedral circumsphere.
*
* @return Radius of circumsphere.
*/

mt_real tet_circumsphere(const vec3r & p0, const vec3r & p1, const vec3r & p2,
                         const vec3r & p3, vec3r & ctr);

/**
* @brief Compute tetrahedral circumsphere.
*
* @param[in] points  Four points spanning the tetrahedron.
* @param[out] ctr Midpoint of tetrahedral circumsphere.
*
* @return Radius of circumsphere.
*/

mt_real tet_circumsphere(const vec3r * points, vec3r & ctr);


/**
* @brief Compute triangle circumcircle.
*
* @param[in] p0  First tri point.
* @param[in] p1  Second tri point.
* @param[in] p2  Third tri point.
* @param[out] ctr Midpoint of triangular circumcircle.
*
* @return Radius of circumcircle.
*/
mt_real tri_circumcircle(const vec3r & p0, const vec3r & p1, const vec3r & p2,
                         vec3r & ctr);

/**
* @brief Compute triangular circumcircle.
*
* @param[in] points  Three points spanning the triangle.
* @param[out] ctr Midpoint of triangular circumcircle.
*
* @return Radius of circumcircle.
*/
mt_real tri_circumcircle(const vec3r * points, vec3r & ctr);

mt_real line_circumcircle(const vec3r & p0, const vec3r & p1, vec3r & ctr);
mt_real line_circumcircle(const vec3r* points, vec3r & ctr);

/**
* @brief Compute the circumsphere of an element.
*
* @param[in] type  Element type.
* @param[in] con   Pointer to connectivity of the element.
* @param[in] xyz   Pointer to start of all coords.
* @param[out] ctr  Midpoint of circumsphere.
*
* @return Radius of circumsphere.
*/
mt_real circumsphere(const elem_t type, const mt_int* con, const mt_real* xyz, vec3r & ctr);


/**
* @brief Compute the centroid (barycenter) of an element.
*
* @param[in] type  Element type.
* @param[in] con   Pointer to connectivity of the element.
* @param[in] xyz   Pointer to start of all coords.
*
* @return Barycenter.
*/
vec3r barycenter(const mt_int npts, const mt_int* con, const mt_real* xyz);

inline vec3r
barycenter(const mt_meshdata & mesh, const mt_int eidx)
{
  return barycenter(mesh.e2n_cnt[eidx], mesh.e2n_con.data() + mesh.e2n_dsp[eidx],
                    mesh.xyz.data());
}

/// compute center-point of an array of coordinates
vec3r vertices_centerpoint(const mt_vector<mt_real> & xyz);
/// compute center-point of an array of coordinates
vec3r vertices_centerpoint(const mt_vector<vec3r> & pts);

/**
 * @brief Compute the seperation distance of a set of points defined by
 *        s = 1/2 * min_{x_i, x_j \in Set} || x_i - x_j ||_2
 * @param [in] xyz The points in xyz-format
 * @param [out] p1 The first point of the set
 * @param [out] p2 The second point of the set
 *
 * @return seperation distance
 */
mt_real vertices_seperation_distance(const mt_vector<mt_real> & xyz, vec3r & p1, vec3r & p2);

/// callback and comfort routine
mt_real vertices_seperation_distance(const mt_vector<vec3r> & xyz, vec3r & p1, vec3r & p2);

/// compute seperation distance on a subset of xyz
mt_real vertices_seperation_distance(const mt_vector<mt_real> & xyz,
                                     const mt_vector<mt_int> & subidx,
                                     vec3r & p1, vec3r & p2);
/// callback and comfort routine
mt_real vertices_seperation_distance(const mt_vector<vec3r> & xyz,
                                     const mt_vector<mt_int> & subidx,
                                     vec3r & p1, vec3r & p2);

/**
* @brief Get bounding sphere for a vector of points.
*
* @param [in]  pts  The points.
* @param [out] ctr  Bounding sphere center.
* @param [out] rad  Bounding sphere radius.
*/
void bounding_sphere(const mt_vector<vec3r> & pts, vec3r & ctr, mt_real & rad);

/**
 * Implementation of radial basis function interpolation of point cloud data
 * see https://en.wikipedia.org/wiki/Radial_basis_function for definition of RBFs
 * see https://doi.org/10.1016/S0377-0427(01)00485-X for paper describing the interpolator
 */
class PointCloudInterpolator
{
  private:
  mt_vector<vec3r> _pts; ///< The pointcloud where data is defined on
  mt_vector<mt_real> _roi; ///< The radius of influence for each point in the input pointcloud
  mt_vector<mt_vector<mt_real> > _rbf_coeffs; ///< The coefficients of the local RBF interpolator
  mt_vector<mt_real> _rbf_scal; ///< The nodewise scaling parameter of of the radial basis function
  mt_vector<mt_vector<mt_int> > _nbh_idx; ///< Array holding the indices of neighbors inside sphere
                                          ///< of given radius
  mt_vector<dmat<mt_real> > _lu_mats; ///< Array holding the LU-Decomps of the local rbf
                                      ///< interpolationn-problem matrices
  int _NW; ///< Amount of neighbors included for calculating radius of influence
  int _NQ; ///< Amount of neighbors included for local RBF interpolator
  int _dpn; ///< Degrees of freedom per input node

  kdtree* _tree; ///< KD-Tree for speeding up lookup of data

  vec3r _bounding_sphere_center;
  mt_real _bounding_sphere_radius_sq;

  /**
   * @brief Compute the radius of influence
   *
   * @param verbose verbosity option
   */
  void compute_roi(const bool verbose = false);
  /**
   * @brief Compute the coefficients of the local RBF interpolator
   *
   * @param data Input data
   * @param verbose verbosity option
   */
  void compute_rbf_coeffs(const mt_vector<mt_real> & data, const bool verbose = false);
  /**
   * @brief Interpolate data using the constructed local RBF interpolator
   *
   * @param inter_pts Points to interpolate to
   * @param odat Output data
   * @param verbose verbosity option
   */
  void interpolate_rbf(const mt_vector<mt_real> & inter_pts, mt_vector<mt_real> & odat, const bool verbose = false) const;


  /**
   * @brief Evaluate the radial basis function
   *
   * @param r Squared radius
   * @param ck scaling coefficient
   * @return mt_real Value of RBF
   */
  mt_real eval_rbf(const mt_real r, const mt_real ck) const;

  public:
  /**
   * @brief Emtpy constructor
   *
   */
  PointCloudInterpolator() :
  _pts(), _roi(), _rbf_coeffs(), _rbf_scal(), _nbh_idx(),
  _NW(0), _NQ(0), _dpn(0), _tree(nullptr)
  {}
  /**
   * @brief Construct a new Point Cloud Interpolator object with input points
   *
   * @param pts Input points
   */
  PointCloudInterpolator(const mt_vector<vec3r> & pts) :
                         _pts(), _roi(), _rbf_coeffs(), _rbf_scal(),
                         _nbh_idx(),
                         _NW(0), _NQ(0), _dpn(0), _tree(nullptr)
  {
    set_pts(pts);
  }
  /**
   * @brief Construct a new Point Cloud Interpolator object with input points in xyz format
   *
   * @param xyz Input points in xyz format
   */
  PointCloudInterpolator(const mt_vector<mt_real> & xyz) :
                         _roi(), _rbf_coeffs(), _rbf_scal(),
                         _nbh_idx(),
                         _NW(0), _NQ(0), _dpn(0), _tree(nullptr)
  {
    mt_vector<vec3r> pts;
    array_to_points(xyz, pts);

    set_pts(pts);
  }
  /**
   * @brief Destroy the Point Cloud Interpolator object
   *
   */
  virtual ~PointCloudInterpolator()
  {
    if(_tree) delete _tree;
  }
  /**
   * @brief Setter for _pts. Also deletes the kdtree and rebuilds it in case the pts have been to            bad
   *
   * @param pts The data points
   */
  void set_pts(const mt_vector<vec3r> & pts)
  {
    const int min_leafs = 50;

    // it might be that we call this->set_pts(this->_pts). In the case we cant do a deep-copy
    if(_pts.data() != pts.data())
      _pts = pts;

    bounding_sphere(_pts, _bounding_sphere_center, _bounding_sphere_radius_sq);
    _bounding_sphere_radius_sq = POW2(_bounding_sphere_radius_sq);

#if 1
    if(_tree) {
      int min, max;
      mt_real avrg;
      _tree->balance_stats(min, max, avrg);
      //fprintf(stderr, "DEBUG in %s: max = %d, min = %d, avrg = %g\n", __func__, max, min, avrg);
      if( max > int(100 * avrg) ) {
        delete _tree;
        _tree = new kdtree(min_leafs);
        _tree->build_vertex_tree(_pts, new cartesian_csys(), 100.0);
      }
      else {
        _tree->clear_vertices();
        mt_vector<mt_real> cpy_pts;
        points_to_array(_pts, cpy_pts);
        _tree->insert_vertices(cpy_pts);
      }
    }
    else {
      _tree = new kdtree(min_leafs);
      _tree->build_vertex_tree(_pts, new cartesian_csys(), 100.0);
    }
#else
    if(_tree) delete _tree;
    _tree = new kdtree(min_leafs);
    _tree->build_vertex_tree(_pts, new cartesian_csys(), 100.0);
#endif
  }

  void set_pts(const mt_vector<mt_real> & pts)
  {
    array_to_points(pts, _pts);
    set_pts(_pts);
  }

  /**
   * @brief Public setup routine for radius of influence
   *
   * @param NW Neighborhood size for calculating radius of influence
   * @param verbose verbosity option
   */
  void setup_roi(const int NW = 7, const bool verbose = false) {
    _NW = NW;
    compute_roi(verbose);
  }
  /**
   * @brief Public setup routine for calculating the local RBF coefficients
   *
   * @param data Input data
   * @param dpn Degrees of freedom per node
   * @param NQ Neighborhood size for the local RBF interpolator
   * @param verbose verbosity option
   */
  void setup_rbf_coeffs(const mt_vector<mt_real> & data, const int dpn, const int NQ = 18, const bool verbose = false) {
    assert(data.size() > 0);
    _dpn = dpn;
    _NQ = NQ+1;
    compute_rbf_coeffs(data, verbose);
  }
  /**
   * @brief Public routine for interpolating data
   *
   * @param inter_pts Points to interpolate to
   * @param odat Output data
   * @param verbose verbosity option
   */
  void interpolate_data(const mt_vector<mt_real> & inter_pts, mt_vector<mt_real> & odat, const bool verbose = false) const
  {
    interpolate_rbf(inter_pts, odat, verbose);
  }
  /**
   * @brief Public routine for interpolating data
   *
   * @param inter_pts Points to interpolate to in vec3r form
   * @param odat Output data
   */
  void interpolate_data(const mt_vector<vec3r> & inter_pts, mt_vector<mt_real> & odat, const bool verbose = false) const
  {
    mt_vector<mt_real> xyz;
    points_to_array(inter_pts, xyz);
    interpolate_data(xyz, odat, verbose);
  }

  /**
   * @brief Interpolate with a modified shepard interpolator
   *
   * @param inter_pts Points to interpolate to
   * @param idat Input data
   * @param dpn degrees of freedom per node
   * @param global use a global shepard interpolator
   * @param odat Output data
   */
  void interpolate_shepard(const mt_vector<mt_real> & inter_pts,
                           const mt_vector<mt_real> &idat,
                           const int dpn, const bool global,
                           mt_vector<mt_real> & odat, const bool verbose = false) const;

  void interpolate_shepard(const mt_vector<vec3r> & inter_pts,
                           const mt_vector<mt_real> &idat,
                           const int dpn, const bool global,
                           mt_vector<mt_real> & odat, const bool verbose = false) const
  {
    mt_vector<mt_real> xyz;
    points_to_array(inter_pts, xyz);
    interpolate_shepard(xyz, idat, dpn, global, odat, verbose);
  }

  /**
   * @brief Recompute the RBF coeffs for new data
   *
   * @param data New data
   */
  void recompute_rbf_coeffs(const mt_vector<mt_real> & data, const bool verbose = false);
};

#endif

