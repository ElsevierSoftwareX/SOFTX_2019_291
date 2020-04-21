/**
* @file kdtree.h
* @brief K-dimensional search tree implementation and related data structs and algorithms
* @author Aurel Neic
* @version
* @date 2018-07-28
*/
#ifndef _KDTREE_H
#define _KDTREE_H

#include <cstdio>
#include <cmath>
#include <cassert>
#include <limits>

#include "mt_utils_base.h"

static const mt_real EPS_MB = 0.000001;
static const mt_real EPS_KDTREEE = 1.0e-4;
static const mt_real FLT_MAX = std::numeric_limits<mt_real>::max();

/**
* @brief Bounding box struct
*/
struct bbox {
  /// the bounds. bounds[0] = lower left (min), bounds[1] = upper right (max)
  vec3r bounds[2];
};

/**
* @brief Triangle element struct used inside the kdtree.
*/
struct tri_elem {
  vec3r   v0, v1, v2;
  mt_int  eidx;        ///< triangle element index
};

/**
* @brief Class holding ray information.
*/
class ray
{
  public:
  vec3r A;   ///< Ray origin.
  vec3r B;   ///< Ray direction.
  vec3r iB;  ///< Ray direction with inversed components. (needed for bbox intersection)
  unsigned char sign[3]; ///< sign of inversed ray direction components (needed for bbox intersection)

  ray(){}
  ray(vec3r a, vec3r b) : A(a), B(b)
  {
    iB.x = 1.0f / B.x;
    iB.y = 1.0f / B.y;
    iB.z = 1.0f / B.z;

    sign[0] = iB.x < 0;
    sign[1] = iB.y < 0;
    sign[2] = iB.z < 0;
  }

  vec3r & origin() {return A;}
  vec3r & direction() {return B;}
  const vec3r & origin() const {return A;}
  const vec3r & direction() const {return B;}

  vec3r point_at_param(mt_real t) const {return A + B*t;}
};
/**
 * @brief Calculate barycentric coordinates for a given point inside a triangle spanned by (v0, v1, v2). Algorithm taken from Numerical Recipes in C++
 *
 * @param v0 First triangle vertex
 * @param v1 Second triangle vertex
 * @param v2 Third triangle vertex
 * @param q Point to calculate barycentric coords for
 * @param u 1st barycentric coord
 * @param v 2nd --"--
 * @param w 3rf --"--
 */
void barycentric_coordinates_triangle(const vec3r & v0, const vec3r & v1, const vec3r & v2, const vec3r & q, mt_real & u, mt_real & v, mt_real & w);

/**
* @brief Intersect a ray with a triangle. Algorithm from Numerical Recipes in C++
*
* @param r     Ray datastruct
* @param v0    First triangle vertex
* @param v1    Second triangle vertex
* @param v2    Third triangle vertex
* @param t     Ray direction scale
* @param u     First triangle barycentric coord
* @param v     Second triangle barycentric coord
*
* @return Whether an intersection was found
*/
bool ray_triangle_intersect(const ray & r,
                            const vec3r &v0, const vec3r &v1, const vec3r &v2,
                            mt_real &t, mt_real &u, mt_real &v);

/**
* @brief Intersect a ray with a triangle.
         optimized Mueller-Trombone algorithm
         https://github.com/erich666/jgt-code/blob/master/Volume_02/Number_1/Moller1997a/raytri.c
         intersect_triangle_1
*
* @param r     Ray datastruct
* @param v0    First triangle vertex
* @param v1    Second triangle vertex
* @param v2    Third triangle vertex
* @param t     Ray direction scale
* @param u     First triangle barycentric coord
* @param v     Second triangle barycentric coord
*
* @return Whether an intersection was found
*/
bool ray_triangle_intersect_mueller_trombone (const ray & r,
                                              const tri_elem & tri,
                                              mt_real &t, mt_real &u, mt_real &v);
#define KDMIN(A, B) (A < B ? A : B)
#define KDMAX(A, B) (A < B ? B : A)

template<typename S>
void generate_bbox(const mt_vector<mt_point<S> > & xyz, bbox & box)
{
  vec3r & min = box.bounds[0];
  vec3r & max = box.bounds[1];

  size_t nnod = xyz.size();
  min.x = xyz[0].x, min.y = xyz[0].y, min.z = xyz[0].z;
  max.x = xyz[0].x, max.y = xyz[0].y, max.z = xyz[0].z;

  for(size_t i=1; i<nnod; i++) {
    S x = xyz[i].x, y = xyz[i].y, z = xyz[i].z;
    if(min.x > x) min.x = x;
    if(min.y > y) min.y = y;
    if(min.z > z) min.z = z;
    if(max.x < x) max.x = x;
    if(max.y < y) max.y = y;
    if(max.z < z) max.z = z;
  }
}

template<typename S>
void generate_bbox(const mt_vector<S> & xyz, bbox & box)
{
  vec3r & min = box.bounds[0];
  vec3r & max = box.bounds[1];

  size_t nnod = xyz.size() / 3;
  min.x = xyz[0*3+0], min.y = xyz[0*3+1], min.z = xyz[0*3+2];
  max.x = xyz[0*3+0], max.y = xyz[0*3+1], max.z = xyz[0*3+2];

  for(size_t i=1; i<nnod; i++) {
    S x = xyz[i*3+0], y = xyz[i*3+1], z = xyz[i*3+2];
    if(min.x > x) min.x = x;
    if(min.y > y) min.y = y;
    if(min.z > z) min.z = z;
    if(max.x < x) max.x = x;
    if(max.y < y) max.y = y;
    if(max.z < z) max.z = z;
  }
}

template<typename T, typename S>
void generate_bbox(const mt_vector<mt_point<S> > & xyz, const mt_vector<T> & nod, bbox & box)
{
  vec3r & min = box.bounds[0];
  vec3r & max = box.bounds[1];

  size_t nnod = nod.size();
  min.x = xyz[nod[0]].x, min.y = xyz[nod[0]].y, min.z = xyz[nod[0]].z;
  max.x = xyz[nod[0]].x, max.y = xyz[nod[0]].y, max.z = xyz[nod[0]].z;

  for(size_t i=1; i<nnod; i++) {
    T v = nod[i];
    S x = xyz[v].x, y = xyz[v].y, z = xyz[v].z;
    if(min.x > x) min.x = x;
    if(min.y > y) min.y = y;
    if(min.z > z) min.z = z;
    if(max.x < x) max.x = x;
    if(max.y < y) max.y = y;
    if(max.z < z) max.z = z;
  }
}

template<typename T, typename S>
void generate_bbox(const mt_vector<S> & xyz, const mt_vector<T> & nod, bbox & box)
{
  vec3r & min = box.bounds[0];
  vec3r & max = box.bounds[1];

  size_t nnod = nod.size();
  min.x = xyz[nod[0]*3+0], min.y = xyz[nod[0]*3+1], min.z = xyz[nod[0]*3+2];
  max.x = xyz[nod[0]*3+0], max.y = xyz[nod[0]*3+1], max.z = xyz[nod[0]*3+2];

  for(size_t i=1; i<nnod; i++) {
    T v = nod[i];
    S x = xyz[v*3+0], y = xyz[v*3+1], z = xyz[v*3+2];
    if(min.x > x) min.x = x;
    if(min.y > y) min.y = y;
    if(min.z > z) min.z = z;
    if(max.x < x) max.x = x;
    if(max.y < y) max.y = y;
    if(max.z < z) max.z = z;
  }
}

/// check if bbox is hit by a ray
bool bbox_is_hit_1(const ray & r, const bbox & box, const mt_real scale_min, const mt_real scale_max);
/// check if bbox is hit by a ray
bool bbox_is_hit_2(const ray & r, const bbox & box, const mt_real scale_min, const mt_real scale_max);
/// check if bbox is hit by a ray
bool bbox_is_hit_3(const ray & r, const bbox & box, const mt_real scale_min, const mt_real scale_max);


/// whether a is subinterval of b
template<typename T>
bool is_subinterval(T a0, T a1, T b0, T b1)
{
  return a0 >= b0 && a1 <= b1;
}

template<typename T>
bool intervals_intersect(T a0, T a1, T b0, T b1)
{
  return a0 < b1 && b0 < a1;
}

vec3r random_point_in_unit_sphere_fast();


/**
* @brief bbox axis enum
*/
enum bboxAxis {
  X = 0,
  Y,
  Z,
  UNSET
};

/**
* @brief Get the longest axis of a bbox
*/
inline bboxAxis
get_longest_axis(const bbox & box)
{
  const vec3r & min = box.bounds[0];
  const vec3r & max = box.bounds[1];

  mt_real x = max.x - min.x;
  mt_real y = max.y - min.y;
  mt_real z = max.z - min.z;

  return x > y && x > z ? X : y > z ? Y : Z;
}


inline mt_real
get_longest_axis_length(const bbox & box)
{
  const vec3r & min = box.bounds[0];
  const vec3r & max = box.bounds[1];

  mt_real x = max.x - min.x;
  mt_real y = max.y - min.y;
  mt_real z = max.z - min.z;

  return x > y && x > z ? (x) : y > z ? (y) : (z);
}

class basic_csys
{
  public:
  /// scale a bounding box around its center by a given factor s
  virtual void bbox_scale_centered(bbox & box, const mt_real s) const = 0;
  /**
  * @brief Check whether box a is subbox of b
  */
  virtual bool is_subbox(const bbox & a, const bbox & b) const = 0;
  /**
  * @brief Check whether two bounding boxes have a non-empty intersection.
  */
  virtual bool bboxes_intersect(const bbox & a, const bbox & b) const = 0;

  virtual bbox bbox_from_sphere(const vec3r ctr, const mt_real rad) const = 0;

  virtual void sphere_from_bbox(const bbox & in, mt_real & rad, vec3r & ctr) const = 0;

  virtual mt_real distance(const vec3r v1, const vec3r v2) const = 0;
  virtual mt_real distance2(const vec3r v1, const vec3r v2) const = 0;

  /// Whether a vert is in a given bounding box
  virtual bool vert_in_bbox(const vec3r & p, const bbox & box) const = 0;

  virtual ~basic_csys() {}

  virtual bboxAxis get_longest_axis(const bbox & box) const = 0;

};

class cartesian_csys : public basic_csys
{
  public:

  inline void bbox_scale_centered(bbox & box, const mt_real s) const
  {
    vec3r ctr = (box.bounds[0] + box.bounds[1]) * 0.5;
    box.bounds[0] = ((box.bounds[0] - ctr) * s) + ctr;
    box.bounds[1] = ((box.bounds[1] - ctr) * s) + ctr;
  }

  inline bool is_subbox(const bbox & a, const bbox & b) const
  {
    bool ret = is_subinterval(a.bounds[0].x, a.bounds[1].x, b.bounds[0].x, b.bounds[1].x) &&
               is_subinterval(a.bounds[0].y, a.bounds[1].y, b.bounds[0].y, b.bounds[1].y) &&
               is_subinterval(a.bounds[0].z, a.bounds[1].z, b.bounds[0].z, b.bounds[1].z);

    return ret;
  }

  inline bool bboxes_intersect(const bbox & a, const bbox & b) const
  {
    bool ret = intervals_intersect(a.bounds[0].x, a.bounds[1].x, b.bounds[0].x, b.bounds[1].x) &&
               intervals_intersect(a.bounds[0].y, a.bounds[1].y, b.bounds[0].y, b.bounds[1].y) &&
               intervals_intersect(a.bounds[0].z, a.bounds[1].z, b.bounds[0].z, b.bounds[1].z);

    return ret;
  }

  inline bbox bbox_from_sphere(const vec3r ctr, const mt_real rad) const
  {
    bbox ret;

    ret.bounds[0] = {ctr.x - rad, ctr.y - rad, ctr.z - rad};
    ret.bounds[1] = {ctr.x + rad, ctr.y + rad, ctr.z + rad};

    return ret;
  }

  inline void sphere_from_bbox(const bbox & in, mt_real & rad, vec3r & ctr) const
  {
    ctr = (in.bounds[0] + in.bounds[1]) * 0.5;

    //Radius of bounding sphere is half of the diagonal length
    vec3r diam = in.bounds[1] - in.bounds[0];
    rad  = 0.5 * diam.length();
  }

  inline bool vert_in_bbox(const vec3r & p, const bbox & box) const
  {
    bool in_x = p.x >= box.bounds[0].x && p.x <= box.bounds[1].x;
    bool in_y = p.y >= box.bounds[0].y && p.y <= box.bounds[1].y;
    bool in_z = p.z >= box.bounds[0].z && p.z <= box.bounds[1].z;

    return in_x && in_y && in_z;
  }

  inline mt_real distance(const vec3r v1, const vec3r v2) const
  {
    return (v1 - v2).length();
  }

  inline mt_real distance2(const vec3r v1, const vec3r v2) const
  {
    return (v1 - v2).length2();
  }

  inline bboxAxis get_longest_axis(const bbox & box) const
  {
    const vec3r & min = box.bounds[0];
    const vec3r & max = box.bounds[1];

    mt_real x = max.x - min.x;
    mt_real y = max.y - min.y;
    mt_real z = max.z - min.z;

    return x > y && x > z ? X : y > z ? Y : Z;
  }

  ~cartesian_csys() {}
};

class cylindrical_csys : public basic_csys
{
  // we use (r, phi, z) coordinate convention

  public:

  #define ANG_LB -MT_PI
  #define ANG_RB  MT_PI

  inline void bbox_scale_centered(bbox & box, const mt_real s) const
  {
    vec3r ctr = (box.bounds[0] + box.bounds[1]) * 0.5;
    box.bounds[0] = ((box.bounds[0] - ctr) * s) + ctr;
    box.bounds[1] = ((box.bounds[1] - ctr) * s) + ctr;

    if(box.bounds[0].y < ANG_LB) box.bounds[0].y = ANG_LB;
    if(box.bounds[1].y > ANG_RB) box.bounds[1].y = ANG_RB;
  }

  inline bool is_subbox(const bbox & a, const bbox & b) const
  {
    bool chk_r   = is_subinterval(a.bounds[0].x, a.bounds[1].x, b.bounds[0].x, b.bounds[1].x);
    bool chk_phi = is_subinterval(a.bounds[0].y, a.bounds[1].y, b.bounds[0].y, b.bounds[1].y);
    bool chk_z   = is_subinterval(a.bounds[0].z, a.bounds[1].z, b.bounds[0].z, b.bounds[1].z);

    return chk_r && chk_phi && chk_z;
  }

  inline bool bboxes_intersect(const bbox & a, const bbox & b) const
  {
    // with cylindrical coordinates, we have the special sitiuation with the angular coordinate:
    // if bbox a crosses ANG_LB or ANG_RB, we have to do two checks instead of one.
    // we expect that bbox b does not cross ANG_LB / ANG_RB.
    bool chk_phi = false;

    if(a.bounds[0].y < ANG_LB) {
      chk_phi = intervals_intersect<mt_real>(a.bounds[0].y+(ANG_RB-ANG_LB), ANG_RB, b.bounds[0].y, b.bounds[1].y) ||
                intervals_intersect<mt_real>(ANG_LB, a.bounds[1].y, b.bounds[0].y, b.bounds[1].y);
    }
    else if (a.bounds[1].y > ANG_RB) {
      chk_phi = intervals_intersect<mt_real>(ANG_LB, a.bounds[0].y-(ANG_RB-ANG_LB), b.bounds[0].y, b.bounds[1].y) ||
                intervals_intersect<mt_real>(a.bounds[0].y, ANG_RB, b.bounds[0].y, b.bounds[1].y);

    }
    else
      chk_phi = intervals_intersect<mt_real>(a.bounds[0].y, a.bounds[1].y, b.bounds[0].y, b.bounds[1].y);

    bool chk_r   = intervals_intersect<mt_real>(a.bounds[0].x, a.bounds[1].x, b.bounds[0].x, b.bounds[1].x);
    bool chk_z   = intervals_intersect<mt_real>(a.bounds[0].z, a.bounds[1].z, b.bounds[0].z, b.bounds[1].z);

    return chk_r && chk_phi && chk_z;
  }

  inline bbox bbox_from_sphere(const vec3r ctr, const mt_real rad) const
  {
    bbox ret;

    ret.bounds[0] = {ctr.x - rad, ctr.y - rad/ctr.x, ctr.z - rad};
    ret.bounds[1] = {ctr.x + rad, ctr.y + rad/ctr.x, ctr.z + rad};

    return ret;
  }

  inline void sphere_from_bbox(const bbox & in, mt_real & rad, vec3r & ctr) const
  {
    ctr = (in.bounds[0] + in.bounds[1]) * 0.5;

    //Radius of bounding sphere is half of the diagonal length
    rad  = 0.5f * distance(in.bounds[1], in.bounds[0]);
  }

  inline bool vert_in_bbox(const vec3r & p, const bbox & box) const
  {
    bool in_x = p.x >= box.bounds[0].x && p.x <= box.bounds[1].x;
    bool in_y = p.y >= box.bounds[0].y && p.y <= box.bounds[1].y;
    bool in_z = p.z >= box.bounds[0].z && p.z <= box.bounds[1].z;

    return in_x && in_y && in_z;
  }

  inline mt_real distance(const vec3r v1, const vec3r v2) const
  {
    mt_real dr = (v1.x - v2.x), dphi = (v1.y - v2.y), dz = (v1.z - v2.z);
    return std::sqrt(POW2(dr) + POW2(v1.x)*POW2(dphi) + POW2(dz));
  }

  inline mt_real distance2(const vec3r v1, const vec3r v2) const
  {
    mt_real dr = (v1.x - v2.x), dphi = (v1.y - v2.y), dz = (v1.z - v2.z);
    return POW2(dr) + POW2(v1.x)*POW2(dphi) + POW2(dz);
  }

  inline bboxAxis get_longest_axis(const bbox & box) const
  {
    const vec3r & min = box.bounds[0];
    const vec3r & max = box.bounds[1];

    mt_real x = max.x - min.x;
    mt_real y = ((min.x + max.x)*0.5f) * (max.y - min.y);
    mt_real z = max.z - min.z;

    return x > y && x > z ? X : y > z ? Y : Z;
  }

  ~cylindrical_csys() {}
};

template<class V>
inline void cartesian_to_cylindrical(const mt_point<V> & crt_pt, mt_point<V> & cyl_pt)
{
  cyl_pt.x = std::sqrt(POW2(crt_pt.x) + POW2(crt_pt.y));
  cyl_pt.y = atan2f(crt_pt.y, crt_pt.x);
  cyl_pt.z = crt_pt.z;
}

template<class V>
inline void cylindrical_to_cartesian(const mt_point<V> & cyl_pt, mt_point<V> & crt_pt)
{
  crt_pt.x = cyl_pt.x * cos(cyl_pt.y);
  crt_pt.y = cyl_pt.x * sin(cyl_pt.y);
  crt_pt.z = cyl_pt.z;
}

// #define KDTREE_VERB
//
/**
* @brief KDTREE class
*/
class kdtree {

  public:

  /**
  * @brief Node definition of a kdtree.
  */
  class node {
    public:
    int up;     ///< index of the upper node
    int left;   ///< index of the node to the left
    int right;  ///< index of the node to the right
    bboxAxis split_axis;  ///< The axis used to split left and right at

    node(): up(-1), left(-1), right(-1), split_axis(UNSET)
    {}
  };

  mt_vector<node> nodes;                    ///< nodes in the kdtree.
  mt_vector<bbox> boxes;                    ///< bboxes in the kdtree.
  mt_vector<mt_vector<tri_elem>*> tris;     ///< triangles in the kdtree.
  mt_vector<mt_vector<int>*>      verts;    ///< vertices in the kdtree.
  int init_size;
  int last_nod_idx;
  int items_per_leaf;
  /** the coordinate system used in vertex search algorithms
   *
   * This allows to abstract away differences between kartesian and cylindrical
   * coordinate systems.
   */
  const basic_csys* csys;

  mt_vector<vec3r> _xyz;                    ///< coordinates

  /**
  * @brief Constructor
  *
  * @param t  Target amount of triangles per kdtree leaf.
  */
  kdtree(int t);
  kdtree() : kdtree(10) {}

  ~kdtree();

  /**
  * @brief Set up the kdtree structure
  *
  * @param triangles  ///< the triangles
  */
  void build_tree(const mt_vector<tri_elem> & triangles);

  /**
  * @brief Set up the kdtree structure
  *
  * @param trimesh The triangle mesh we set the kdtree up for
  */
  void build_tree(const mt_meshdata & trimesh);

  /**
  * @brief Build a kdtree for vertices.
  *
  * @param xyz         The vertices we want to build the kdtree for.
  * @param bbox_scale  The factor to scale the starting bounding box with. This allows the kdtree to span a bigger volume, thus
  *                    preventing bbox misses on queries.
  */
  void build_vertex_tree(const mt_vector<mt_real> & xyz, const basic_csys* used_csys = new cartesian_csys(),
                         const mt_real bbox_scale = 2.0f);

  /// build_vertex_tree wrapper for arbitary vertex datatypes
  template<class V> inline
  void build_vertex_tree(const mt_vector<mt_point<V> > & xyz,
                         const basic_csys* used_csys = new cartesian_csys(),
                         const mt_real bbox_scale = 2.0f)
  {
    mt_vector<mt_real> xyzbuff;
    points_to_array(xyz, xyzbuff);
    build_vertex_tree(xyzbuff, used_csys, bbox_scale);
  }

  /// get kdtree balance statistics
  void balance_stats(int & min_per_leaf, int & max_per_leaf, mt_real & avrg_per_leaf);

  /**
  * @brief Get the closest intersection between a ray and the triangles in the kdtree.
  *
  * @param r               The ray we check the intersection for.
  * @param scale_min       Minimum valid scale.
  * @param scale_max       Maximum valid scale.
  * @param hit_ele         The triangle element with the current closest hit.
  * @param hit_pos         The position of the current closest hit.
  *
  * @return
  */
  bool closest_intersect(const ray & r, const mt_real scale_min, const mt_real scale_max,
                         tri_elem & hit_ele, vec3r & hit_pos) const;

  /**
  * @brief Get the number of intersections between a ray and a kdtree of triangles.
  *
  * @param r               The ray we check the intersection for.
  * @param scale_min       Minimum valid scale.
  * @param scale_max       Maximum valid scale.
  *
  * @return The number of intersections.
  */
  int count_intersects(const ray & r, const mt_real scale_min, const mt_real scale_max) const;

  /**
  * @brief Find the vertex closest to a reference point.
  *
  * @param [in]  ref     Reference point.
  * @param [out] idx     Index of closest vertex.
  * @param [out] closest Coords of closest vertex.
  * @param [out] len2    Squared length to closest vertex.
  *
  * @return Whether a closest vertex was found.
  */
  bool closest_vertex(const vec3r ref, int & idx, vec3r & closest, mt_real & len2) const;

  /**
  * @brief Find the k-closest vertices to a reference point.
  *
  * @param [in]  k       Amount of vertices
  * @param [in]  ref     Reference point.
  * @param [out] vtx     set holding the squared distances and the indices
  */
  void k_closest_vertices(const int k, const vec3r ref, mt_vector<mixed_tuple<mt_real, int> > & vtx) const;

  /**
  * @brief Find all vertices in a given sphere definition.
  *
  * @param [in]  ref   The center of the sphere.
  * @param [in]  rad   The radius of the sphere.
  * @param [out] vtx   The vertices inside the sphere.
  */
  void vertices_in_sphere(const vec3r ref, const mt_real rad, mt_vector<mixed_tuple<mt_real,int> > & vtx) const;

  /// clear the current vertices, but keep the volume hirarchy
  void clear_vertices();

  /**
  * @brief Insert a set of vertices into the kdtree.
  *
  * Note that the old vertices are still present in the kdtree. As such, the new indexing
  * of the new vertices is implicitely after the indexing of the old ones. Use clear_vertices()
  * berforehand if you just want to replace the old vertex coordinates with new ones.
  *
  * @param xyz  A vector holding the new set of coordinates to insert.
  */
  void insert_vertices(const mt_vector<mt_real> & xyz);

  /**
   * @brief Check whether a reference point is in the root bbox
   * @param ref A test vector
   *
   * @return True if point is inside the bbox False else
   */

  inline bool in_root_bbox(const vec3r ref) const
  {
    const bbox & root_bbox = boxes[0];

    if(csys)
      return csys->vert_in_bbox(ref, root_bbox);
    else {
      fprintf(stderr, "%s error: No coordinate system set up! Aborting!\n", __func__);
      exit(1);
    }
  }

  private:
  /**
  * @brief Constructs a kdtree node and recurses to the left and right
  *        nodes if possible.
  *
  * @param xyz           Global vertex coords.
  * @param current_tris  Triangels of current node.
  * @param depth         kdtree level
  *
  * @return The generated node index.
  */
  int build(mt_vector<tri_elem>* current_tris, int depth, int up);

  /**
  * @brief Build kdtree for vertices
  *
  * @param current_vert  The pool of vertices on the current level.
  * @param depth         The level index.
  *
  * @return The generated node index.
  */
  int build(mt_vector<int>* current_vert, const bbox & box, int up);

  /**
  * @brief Add a nod to the kdtree data arrays.
  *
  * @return Index of the newly created node.
  */
  int add_node();

  /**
  * @brief Recursive function. Searches for the minimum-distance hit between a ray and the
  *        triangles in the kdtree.
  *
  *        To start at root, call with nod_idx = 0.
  *
  * @param nod_idx         Current node index in the kdtree.
  * @param r               The ray we check the intersection for.
  * @param scale_min       Minimum valid scale.
  * @param closest_so_far  The scale of the current closesest hit.
  * @param hit_ele         The triangle element with the current closest hit.
  * @param hit_pos         The position of the current closest hit.
  *
  * @return Whether any hit was detected.
  */
  bool intersect_ray(const int nod_idx, const ray & r, const mt_real scale_min,
                     mt_real & closest_so_far, tri_elem & hit_ele, vec3r & hit_pos) const;
  /**
  * @brief Recursive function. Counts the number of triangles that intersect with a ray.
  *
  *        To start at root, call with nod_idx = 0.
  *
  * @param nod_idx         Current node index in the kdtree.
  * @param r               The ray we check the intersection for.
  * @param scale_min       Minimum valid scale.
  * @param scale_max       Maximum valid scale.
  * @param count           The current number of intersections.
  *
  */
  void count_ray_intersects(const int nod_idx, const ray & r, const mt_real scale_min,
                           const mt_real scale_max, int & count) const;

  /**
  * @brief Get the index of the leaf node that has the reference vertex in its bounding box.
  *
  * @param nod_idx  The node we start searching.
  * @param ref      The reference coordinate we search the closest vertex for.
  * @param depth    The search depth counter. Not really necessary (may always start at 0).
  * @param leaf_idx The leaf node index.
  *
  * @return Whether a the reference coordinate was located insied a leaf node.
  */
  bool find_enclosing_leaf(const int nod_idx, const vec3r ref, int & leaf_idx) const;

  /**
  * @brief Recurse over parts of a kdtree, search EVERY leaf for the closest vertex.
  *
  * @param [in]  nod_idx   Starting node index.
  * @param [in]  ref       Reference vertex.
  * @param [out] idx       Index of current closest vert.
  * @param [out] len2      Distance to reference of the current closest vertex.
  */
  void closest_vertex(const int nod_idx, const vec3r ref, int & idx, mt_real & len2) const;
  /**
  * @brief Get the closest vertex from the leafs that overlap with a given region
  *
  * @param [in]  nod_idx  Starting node index.
  * @param [in]  ref      Reference vertex.
  * @param [in]  reg      Search region.
  * @param [out] idx      Index of the current closest vertex.
  * @param [out] len2     Distance to reference of the current closest vertex.
  */
  void closest_vertex_in_region(const int nod_idx, const vec3r ref, const bbox & reg,
                                int & idx, mt_real & len2) const;

  void vertices_in_sphere(const int nidx, const vec3r ref, const mt_real rad,
                          const bbox & search_region, mt_vector<mixed_tuple<mt_real,int> > & vtx) const;

  /**
  * @brief Get the median value for a set of vertices.
  *
  * @param vert         The vertices.
  * @param split_axis   The axis we use to compute the median on.
  *
  * @return The median value.
  */
  mt_real get_median(const mt_vector<int> & vert, const bboxAxis split_axis) const;

  /**
  * @brief Get the index of the first upward node that encloses a given reference
  *        bounding box. Return root index if given bbox is not enclosed by any node.
  *
  * @param start_idx  Starting node index.
  * @param ref_box    The reference bounding box we compare for.
  *
  * @return The enclosing node index.
  */
  int enclosing_node_upwards(const int start_idx, const bbox & ref_box) const;

  /**
  * @brief Get the index of the last downward node that encloses a given reference
  *        bounding box. Return start_idx index if given bbox is not enclosed by any downward node.
  *
  * @param start_idx  Starting node index.
  * @param ref_box    The reference bounding box we compare for.
  *
  * @return The enclosing node index.
  */
  int enclosing_node_downwards(const int start_idx, const bbox & ref_box) const;

};

/**
* @brief Check whether a point is inside a closed surface.
*
* @param tree   Kdtree holding the surface triangles
* @param pos    The point position we test
*
* @return       Whether the point is inside the surface.
*/
bool inside_closed_surface(const kdtree & tree, const vec3r pos);

#endif

