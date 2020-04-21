/**
* @file kdtree.cpp
* @brief K-dimensional search tree implementation and related data structs and algorithms
* @author Aurel Neic
* @version
* @date 2018-07-28
*/
#include "io_utils.h"
#include "kdtree.h"

void generate_bbox(const mt_vector<tri_elem> & tri, bbox & box)
{
  vec3r & min = box.bounds[0];
  vec3r & max = box.bounds[1];

  size_t nele = tri.size();
  min.x = tri[0].v0.x, min.y = tri[0].v0.y, min.z = tri[0].v0.z;
  max.x = tri[0].v0.x, max.y = tri[0].v0.y, max.z = tri[0].v0.z;

  for(size_t i=0; i<nele; i++) {
    mt_real x = tri[i].v0.x, y = tri[i].v0.y, z = tri[i].v0.z;
    if(min.x > x) min.x = x;
    if(min.y > y) min.y = y;
    if(min.z > z) min.z = z;
    if(max.x < x) max.x = x;
    if(max.y < y) max.y = y;
    if(max.z < z) max.z = z;

    x = tri[i].v1.x, y = tri[i].v1.y, z = tri[i].v1.z;
    if(min.x > x) min.x = x;
    if(min.y > y) min.y = y;
    if(min.z > z) min.z = z;
    if(max.x < x) max.x = x;
    if(max.y < y) max.y = y;
    if(max.z < z) max.z = z;

    x = tri[i].v2.x, y = tri[i].v2.y, z = tri[i].v2.z;
    if(min.x > x) min.x = x;
    if(min.y > y) min.y = y;
    if(min.z > z) min.z = z;
    if(max.x < x) max.x = x;
    if(max.y < y) max.y = y;
    if(max.z < z) max.z = z;
  }
}

bool bbox_is_hit_1(const ray & r, const bbox & box, const mt_real scale_min, const mt_real scale_max)
{
  const vec3r ray_orig = r.origin();
  const vec3r ray_dir  = r.direction();

  const vec3r & min = box.bounds[0];
  const vec3r & max = box.bounds[1];

  mt_real tmin, tmax, tymin, tymax, tzmin, tzmax;
  mt_real divx = 1.0f / ray_dir.x,
        divy = 1.0f / ray_dir.y,
        divz = 1.0f / ray_dir.z;

  if (r.direction().x >= 0) {
    tmin = (min.x - ray_orig.x) * divx;
    tmax = (max.x - ray_orig.x) * divx;
  }
  else {
    tmin = (max.x - ray_orig.x) * divx;
    tmax = (min.x - ray_orig.x) * divx;
  }

  if (r.direction().y >= 0) {
    tymin = (min.y - ray_orig.y) * divy;
    tymax = (max.y - ray_orig.y) * divy;
  }
  else {
    tymin = (max.y - ray_orig.y) * divy;
    tymax = (min.y - ray_orig.y) * divy;
  }

  if ( (tmin > tymax) || (tymin > tmax) )
    return false;

  if (tymin > tmin) tmin = tymin;
  if (tymax < tmax) tmax = tymax;

  if (r.direction().z >= 0) {
    tzmin = (min.z - ray_orig.z) * divz;
    tzmax = (max.z - ray_orig.z) * divz;
  }
  else {
    tzmin = (max.z - ray_orig.z) * divz;
    tzmax = (min.z - ray_orig.z) * divz;
  }

  if ( (tmin > tzmax) || (tzmin > tmax) )
    return false;

  if (tzmin > tmin) tmin = tzmin;
  if (tzmax < tmax) tmax = tzmax;

  return ( (tmin < scale_max) && (tmax > scale_min) );
}

bool bbox_is_hit_2(const ray & r, const bbox & box, const mt_real scale_min, const mt_real scale_max)
{
  const vec3r ray_orig = r.origin();
  const vec3r ray_dir  = r.direction();

  mt_real divx = 1.0 / ray_dir.x;
  mt_real divy = 1.0 / ray_dir.y;
  mt_real divz = 1.0 / ray_dir.z;

  const vec3r & min = box.bounds[0];
  const vec3r & max = box.bounds[1];

  mt_real t1 = (min.x - ray_orig.x)*divx;
  mt_real t2 = (max.x - ray_orig.x)*divx;
  mt_real t3 = (min.y - ray_orig.y)*divy;
  mt_real t4 = (max.y - ray_orig.y)*divy;
  mt_real t5 = (min.z - ray_orig.z)*divz;
  mt_real t6 = (max.z - ray_orig.z)*divz;

  mt_real tmin = KDMAX(KDMAX(KDMIN(t1, t2), KDMIN(t3, t4)), KDMIN(t5, t6));
  mt_real tmax = KDMIN(KDMIN(KDMAX(t1, t2), KDMAX(t3, t4)), KDMAX(t5, t6));

  // if tmax < 0, ray is intersecting AABB, but the whole AABB is behind us
  if (tmax < 0) {
    // *t = tmax;
    return false;
  }

  // if tmin > tmax, ray doesn't intersect AABB
  if (tmin > tmax) {
    // *t = tmax;
    return false;
  }

  // *t = tmin;
  // return true;
  return ( (tmin < scale_max) && (tmax > scale_min) );
}

bool bbox_is_hit_3(const ray & r, const bbox & box, const mt_real scale_min, const mt_real scale_max)
{
  mt_real tmin, tmax, tymin, tymax, tzmin, tzmax;
  const vec3r & r_orig = r.origin(), & r_idir = r.iB;

  tmin  = (box.bounds[r.sign[0]    ].x - r_orig.x) * r_idir.x;
  tmax  = (box.bounds[1 - r.sign[0]].x - r_orig.x) * r_idir.x;
  tymin = (box.bounds[r.sign[1]    ].y - r_orig.y) * r_idir.y;
  tymax = (box.bounds[1 - r.sign[1]].y - r_orig.y) * r_idir.y;

  if ( (tmin > tymax) || (tymin > tmax) )
    return false;

  if (tymin > tmin) tmin = tymin;
  if (tymax < tmax) tmax = tymax;

  tzmin = (box.bounds[r.sign[2]    ].z - r_orig.z) * r_idir.z;
  tzmax = (box.bounds[1 - r.sign[2]].z - r_orig.z) * r_idir.z;

  if ( (tmin > tzmax) || (tzmin > tmax) )
    return false;

  if (tzmin > tmin) tmin = tzmin;
  if (tzmax < tmax) tmax = tzmax;

  return ( (tmin < scale_max) && (tmax > scale_min) );
}

void barycentric_coordinates_triangle(const vec3r & v0, const vec3r & v1, const vec3r & v2,
                                      const vec3r & q, mt_real & u, mt_real & v, mt_real & w)
{
  const vec3r ap = v0 - v2;
  const vec3r bp = v1 - v2;
  const vec3r qp = q - v2;

  mt_real apsq = ap.length2();
  mt_real bpsq = bp.length2();
  mt_real apdotbp = ap.scaProd(bp);
  mt_real denom = apsq * bpsq - apdotbp*apdotbp; //only zero for degenerate triangles

  u = (bpsq * ap.scaProd(qp) - ap.scaProd(bp) * bp.scaProd(qp)) / denom;
  v = (apsq * bp.scaProd(qp) - ap.scaProd(bp) * ap.scaProd(qp)) / denom;
  w = 1 - u - v;
}

bool ray_triangle_intersect(const ray & r,
                            const tri_elem & tri,
                            mt_real &t, mt_real &u, mt_real &v)
{
  const vec3r & orig = r.origin();
  const vec3r & dir  = r.direction();
  vec3r qvec;
  // find vectors for two edges sharing v0
  const vec3r edge1 = tri.v1 - tri.v0;
  const vec3r edge2 = tri.v2 - tri.v0;

  // begin calculating determinant - also used to calculate u parameter
  const vec3r pvec = dir.crossProd(edge2);
  // if determinant is near zero, ray lies in plane of triangle
  const mt_real det = edge1.scaProd(pvec);

  if (det > EPS_MB) {
    // calculate distance from v0 to ray origin
    const vec3r tvec = orig - tri.v0;
    // calculate u parameter and test bounds
    u = tvec.scaProd(pvec);
    if (u < 0.0 || u > det)
      return false;

    // prepare to test v parameter
    qvec = tvec.crossProd(edge1);
    // calculate v parameter and test bounds
    v = dir.scaProd(qvec);
    if (v < 0.0 || u + v > det)
      return false;
  }
  else if(det < -EPS_MB) {
    // calculate distance from vert0 to ray origin
    const vec3r tvec = orig - tri.v0;
    // calculate u parameter and test bounds
    u = tvec.scaProd(pvec);
    if (u > 0.0 || u < det)
      return false;
    // prepare to test V parameter
    qvec = tvec.crossProd(edge1);
    // calculate v parameter and test bounds
    v = dir.scaProd(qvec);
    if (v > 0.0 || u + v < det)
      return false;
  }
  else {
      return false;  // ray is parallel to the plane of the triangle
  }

  const mt_real inv_det = 1.0 / det;
  // calculate t, ray intersects triangle
  t = edge2.scaProd(qvec) * inv_det;
  u *= inv_det;
  v *= inv_det;
  return true;
}

kdtree::kdtree(int t) : init_size(1000), last_nod_idx(0), items_per_leaf(t), csys(nullptr)
{
  nodes.resize(init_size);
  boxes.resize(init_size);
}

kdtree::~kdtree()
{
  if(tris.size())
    for(int i=0; i<last_nod_idx; i++)
      if(tris[i]) delete tris[i];

  if(verts.size())
    for(int i=0; i<last_nod_idx; i++)
      if(verts[i]) delete verts[i];

  if(csys) delete csys;
}

void kdtree::build_tree(const mt_vector<tri_elem> & triangles)
{
  tris.resize(init_size, NULL);

  if(!csys)
    csys = new cartesian_csys();

  mt_vector<tri_elem>* tri = new mt_vector<tri_elem>(triangles);
  build(tri, -1, 0);
}

void kdtree::build_tree(const mt_meshdata & trimesh)
{
  tris.resize(init_size, NULL);

  if(!csys)
    csys = new cartesian_csys();

  mt_vector<tri_elem>* triangles = new mt_vector<tri_elem>(trimesh.e2n_cnt.size());
  const mt_real* xyz = trimesh.xyz.data();
  for(size_t i=0; i<triangles->size(); i++) {
    tri_elem & t = (*triangles)[i];
    t.v0.get(xyz + trimesh.e2n_con[i*3+0]*3);
    t.v1.get(xyz + trimesh.e2n_con[i*3+1]*3);
    t.v2.get(xyz + trimesh.e2n_con[i*3+2]*3);
    t.eidx = i;
  }

  build(triangles, -1, 0);
}

void kdtree::build_vertex_tree(const mt_vector<mt_real> & xyz, const basic_csys* used_csys,
                               const mt_real bbox_scale)
{
  csys = used_csys;

  array_to_points(xyz, _xyz);
  verts.resize(init_size, NULL);

  mt_vector<int>* vertices = new mt_vector<int>(xyz.size() / 3);

  for(size_t i=0; i<vertices->size(); i++)
    (*vertices)[i] = i;

  bbox startbox;
  generate_bbox(_xyz, *vertices, startbox);
  csys->bbox_scale_centered(startbox, bbox_scale);

  build(vertices, startbox, -1);
}

int kdtree::build(mt_vector<tri_elem>* current_tris, int up, int depth)
{
  int n        = add_node();
  nodes[n].up  = up;
  int num_poly = current_tris->size();

  if (num_poly == 0) {
    assert(0);
#ifdef KDTREE_VERB
    printf("%d: stopping empty. does this even happen? \n", depth);
#endif
    return n;
  }

  generate_bbox(*current_tris, boxes[n]);
  bbox & box = boxes[n];

  if (num_poly <= items_per_leaf) {
    assert(num_poly > 0);

    tris[n] = current_tris;

#ifdef KDTREE_VERB
    printf("%d: stopping with %d poly. success.\n", depth, (int) num_poly);
#endif
    return n;
  }

  // printf("%d: bbox: min (%lf, %lf, %lf), max (%lf, %lf, %lf)\n", depth,
  // node->bbox->start.x, node->bbox->start.y, node->bbox->start.z,
  // node->bbox->end.x, node->bbox->end.y, node->bbox->end.z);

  vec3r midPoint = (box.bounds[0] + box.bounds[1]) * 0.5;
  mt_vector<tri_elem>* tris_left  = new mt_vector<tri_elem>();
  mt_vector<tri_elem>* tris_right = new mt_vector<tri_elem>();
  tris_left->reserve(num_poly / 2);
  tris_right->reserve(num_poly / 2);

  bboxAxis axis       = get_longest_axis(box);
  nodes[n].split_axis = axis;

  for (const tri_elem t : *(current_tris)) {
    const vec3r polyMidPoint = (t.v0 + t.v1 + t.v2) / mt_real(3.0);

    if (((axis == X) && (midPoint.x < polyMidPoint.x)) ||
        ((axis == Y) && (midPoint.y < polyMidPoint.y)) ||
        ((axis == Z) && (midPoint.z < polyMidPoint.z)))
    {
      tris_right->push_back(t);
    }
    else {
      tris_left->push_back(t);
    }
  }

  bool empty_side = false;
  if ( (tris_left ->size() == 0 && tris_right->size() > 0) ||
      (tris_right->size() == 0 && tris_left ->size() > 0) )
    empty_side = true;

  if (empty_side == false) {
#ifdef KDTREE_VERB
    printf("%d: splitting %d -> ( %d , %d )\n", depth, num_poly, int(tris_left->size()),
        int(tris_right->size()));
    fflush(stdout);
#endif
    //Recurse down both left and right sides
    int lidx = build(tris_left, n, depth + 1);
    int ridx = build(tris_right, n, depth + 1);
    nodes[n].left  = lidx;
    nodes[n].right = ridx;
    // this is not a leaf node -> we can delete the current_tris
    delete current_tris;
  }
  else {
    //Stop here
    tris[n] = current_tris;

#ifdef KDTREE_VERB
    printf("%d: stopping with %d polys because of bad splitting (one side empty). \n", depth, num_poly);
    fflush(stdout);
#endif
  }

  return n;
}

mt_real kdtree::get_median(const mt_vector<int> & vert, const bboxAxis split_axis) const
{
  const size_t sz = vert.size();
  mt_vector<mt_real> med_coords(sz);

  switch(split_axis) {
    case X:
      for(size_t i=0; i < sz; i++)
        med_coords[i] = _xyz[vert[i]].x;
      break;
    case Y:
      for(size_t i=0; i < sz; i++)
        med_coords[i] = _xyz[vert[i]].y;
      break;
    case Z:
      for(size_t i=0; i < sz; i++)
        med_coords[i] = _xyz[vert[i]].z;
      break;
    default: break;
  }

  size_t med = med_coords.size() / 2;
  // we use std::nth_element to get the sorted value at median index without sorting
  // the whole sequence
  std::nth_element(med_coords.begin(), med_coords.begin() + med, med_coords.end());

  return med_coords[med];
}


int kdtree::build(mt_vector<int>* current_vert, const bbox & box, const int up)
{
  int n         = add_node();
  int num_verts = current_vert->size();
  nodes[n].up   = up;
  boxes[n]      = box;

  if (num_verts == 0) {
#ifdef KDTREE_VERB
    printf("kdtree::vertex_build: stopping empty. does this even happen? \n");
#endif
    assert(0);
    return n;
  }

  if (num_verts <= items_per_leaf) {
#ifdef KDTREE_VERB
    printf("kdtree::vertex_build: stopping with %d verts. success.\n", (int) num_verts);
#endif
    verts[n] = current_vert;
    return n;
  }

  mt_vector<int>* verts_left  = new mt_vector<int>();
  mt_vector<int>* verts_right = new mt_vector<int>();
  verts_left->reserve(num_verts / 2);
  verts_right->reserve(num_verts / 2);

  bbox auxbox = box;
  bool unbalanced = false;
  short iter = 0;
  mt_real median = 0.0;

  do {
    iter++;
    bboxAxis axis = csys->get_longest_axis(auxbox);
    nodes[n].split_axis = axis;

    //vec3r midPoint = (auxbox.bounds[0] + auxbox.bounds[1]) * 0.5f;
    median = get_median(*current_vert, axis);
    verts_left->resize(0);
    verts_right->resize(0);

    const mt_vector<int> & v = *(current_vert);
    for (int i = 0; i < num_verts; i++)
    {
      const vec3r & p = _xyz[v[i]];

      if (((axis == X) && (median < p.x)) ||
          ((axis == Y) && (median < p.y)) ||
          ((axis == Z) && (median < p.z)))
      {
        verts_right->push_back(v[i]);
      }
      else {
        verts_left->push_back(v[i]);
      }
    }

    #if 1
    mt_real ls = verts_left ->size(), rs = verts_right->size();
    unbalanced = (fabs(ls - rs) / num_verts) > 0.6f;
    #else
    unbalanced = (verts_left ->size() == 0 && verts_right->size() > 0) ||
                 (verts_right->size() == 0 && verts_left ->size() > 0);
    #endif

    // if we have an unbalanced kdtree node, we want to redo the partitioning for the
    // next, biggest axis. For that, we zero out the current axis bounds, so that the
    // next call to get_longest_axis will give us the next longest
    if(unbalanced) {
      switch(axis) {
        case X: auxbox.bounds[0].x = auxbox.bounds[1].x = 0.0f; break;
        case Y: auxbox.bounds[0].y = auxbox.bounds[1].y = 0.0f; break;
        case Z: auxbox.bounds[0].z = auxbox.bounds[1].z = 0.0f; break;
        default: break;
      }
    }
  } while(unbalanced == true && iter < 3);

  if (unbalanced == false) {
    // compute bboxes for children
    bbox left_split = box, right_split = box;
    switch(nodes[n].split_axis) {
      case X:
        left_split.bounds[1].x = median; right_split.bounds[0].x = median;
        break;
      case Y:
        left_split.bounds[1].y = median; right_split.bounds[0].y = median;
        break;
      case Z:
        left_split.bounds[1].z = median; right_split.bounds[0].z = median;
        break;
      case UNSET: break;
    }

    //Recurse down both left and right sides
    int lidx = build(verts_left, left_split, n);
    int ridx = build(verts_right, right_split, n);
    nodes[n].left  = lidx;
    nodes[n].right = ridx;
    // this is not a leaf node -> we can delete the current_vert
    delete current_vert;
  }
  else {
    //Stop here
    verts[n] = current_vert;
    nodes[n].split_axis = UNSET;

    if(num_verts > items_per_leaf * 2)
      printf("Warning! Node %d: stopping with %d vertices because of unbalanced splitting.\n", n, num_verts);
#ifdef KDTREE_VERB
    printf("Warning! Node %d: stopping with %d vertices because of unbalanced splitting.\n", n, num_verts);
    fflush(stdout);
#endif
  }

  return n;
}


int kdtree::add_node()
{
  int ret = last_nod_idx;
  last_nod_idx++;

  if(last_nod_idx == int(nodes.size()))
  {
    nodes.resize(last_nod_idx*2);
    boxes.resize(last_nod_idx*2);
    if(tris.size())
      tris.resize (last_nod_idx*2, NULL);
    if(verts.size())
      verts.resize (last_nod_idx*2, NULL);
  }

  return ret;
}

bool kdtree::intersect_ray(const int nod_idx, const ray & r, const mt_real scale_min,
                           mt_real & closest_so_far, tri_elem & hit_ele, vec3r & hit_pos) const
{
  const bbox & box       = boxes[nod_idx];
  const kdtree::node & n = nodes[nod_idx];

  if (bbox_is_hit_3(r, box, scale_min, closest_so_far)) {
    bool hasHit = false;
    if (tris[nod_idx] == NULL) {
      // current node is not a leaf yet -> recurse down both sides
      bool hitLeft  = intersect_ray(n.left,  r, scale_min, closest_so_far, hit_ele, hit_pos);
      bool hitRight = intersect_ray(n.right, r, scale_min, closest_so_far, hit_ele, hit_pos);

      return hitLeft || hitRight;
    }
    else {
      // current node is a leaf -> check all polys
      const mt_vector<tri_elem> & tr = *(tris[nod_idx]);
      int num_tri = tr.size();

      for (int i = 0; i < num_tri; i++)
      {
        double u, v, scale;
        bool did_hit = ray_triangle_intersect(r, tr[i], scale, u, v);
        if (did_hit && scale > scale_min && scale < closest_so_far) {
          hasHit = true;
          hit_ele = tr[i];
          hit_pos = r.point_at_param(scale);
          closest_so_far = scale;
        }
      }
      return hasHit;
    }
  }
  else
    return false;
}

bool kdtree::closest_intersect(const ray & r, const mt_real scale_min, const mt_real scale_max,
                               tri_elem & hit_ele, vec3r & hit_pos) const
{
  mt_real closest_so_far = scale_max;
  bool ret = intersect_ray(0, r, scale_min, closest_so_far, hit_ele, hit_pos);
  return ret;
}

bool kdtree::find_enclosing_leaf(const int nod_idx, const vec3r ref, int & leaf_idx) const
{
  const bbox & box       = boxes[nod_idx];
  const kdtree::node & n = nodes[nod_idx];

  if (csys->vert_in_bbox(ref, box)) {
    if (verts[nod_idx] == NULL) {
      // current node is not a leaf yet -> recurse down both sides
      bool in_left  = find_enclosing_leaf(n.left,  ref, leaf_idx);
      bool in_right = find_enclosing_leaf(n.right, ref, leaf_idx);

      if(!(in_left || in_right)) {
        fprintf(stderr, "%s error: KDtree node %d: Vertex is inside me but not inside my sons!\n",
                __func__, nod_idx);
      }
      return in_left || in_right;
    }
    else {
      leaf_idx = nod_idx;
      return true;
    }
  }
  else
    return false;
}

void kdtree::closest_vertex(const int nod_idx, const vec3r ref, int & idx, mt_real & len2) const
{
  const kdtree::node & n = nodes[nod_idx];

  if (verts[nod_idx] == NULL) {
    // current node is not a leaf yet -> recurse down both sides
    closest_vertex(n.left,  ref, idx, len2);
    closest_vertex(n.right, ref, idx, len2);
  }
  else {
    // current node is a leaf -> check all verts
    const mt_vector<int> & v = *verts[nod_idx];
    int num_verts = v.size();

    for (int i = 0; i < num_verts; i++) {
      const vec3r p = _xyz[v[i]];
      mt_real dist2 = csys->distance2(ref, p);

      if (dist2 < len2) {
        len2 = dist2;
        idx = v[i];
      }
    }
  }
}

void kdtree::clear_vertices()
{
  _xyz.resize(0);

  for(size_t i=0; i<verts.size(); i++)
    if(verts[i] != NULL) verts[i]->resize(0);
}

void kdtree::insert_vertices(const mt_vector<mt_real> & xyz)
{
  mt_vector<vec3r> newcoords;
  array_to_points(xyz, newcoords);

  mt_vector<int> newinds(newcoords.size());

  size_t old_numverts = _xyz.size();
  _xyz.append(newcoords.begin(), newcoords.end());

  for(size_t i=0; i<newinds.size(); i++)
    newinds[i] = old_numverts + i;

  for(size_t i = 0; i < newcoords.size(); i++) {
    int idx   = newinds[i];
    vec3r & v = newcoords[i];

    int leaf_idx;
    bool is_inside = find_enclosing_leaf(0, v, leaf_idx);

    if(is_inside) {
      verts[leaf_idx]->push_back(idx);
    }
    else {
      fprintf(stderr, "%s warning: Could not insert vertex since it is outside of kdtree bounding box!\n",
              __func__);
    }
  }
}

void kdtree::balance_stats(int & min_per_leaf, int & max_per_leaf, mt_real & avrg_per_leaf)
{
  min_per_leaf = 100000;
  max_per_leaf = 0;
  avrg_per_leaf = 0.0;
  mt_real numleafs = 0.0;

  if(verts.size() > 0) {
    for(size_t i=0; i<verts.size(); i++) {
      if(verts[i] != NULL) {
        int sz = verts[i]->size();
        if(min_per_leaf > sz) min_per_leaf = sz;
        if(max_per_leaf < sz) max_per_leaf = sz;
        avrg_per_leaf += sz;
        numleafs += 1.0f;
      }
    }
    avrg_per_leaf /= numleafs;
  }
  else {
    for(size_t i=0; i<tris.size(); i++) {
      if(tris[i] != NULL) {
        int sz = tris[i]->size();
        if(min_per_leaf > sz) min_per_leaf = sz;
        if(max_per_leaf < sz) max_per_leaf = sz;
        avrg_per_leaf += sz;
        numleafs += 1.0f;
      }
    }
    avrg_per_leaf /= numleafs;
  }
}

void kdtree::closest_vertex_in_region(const int nod_idx, const vec3r ref, const bbox & reg,
                                      int & idx, mt_real & len2) const
{
  const bbox & box       = boxes[nod_idx];
  const kdtree::node & n = nodes[nod_idx];

  if(csys->bboxes_intersect(reg, box)) {
    if (verts[nod_idx] == NULL) {
      // current node is not a leaf yet -> recurse down both sides
      closest_vertex_in_region(n.left,  ref, reg, idx, len2);
      closest_vertex_in_region(n.right, ref, reg, idx, len2);
    }
    else {
      // current node is a leaf -> check all verts
      const mt_vector<int> & v = *verts[nod_idx];
      int num_verts = v.size();

      for (int i = 0; i < num_verts; i++) {
        const vec3r p = _xyz[v[i]];
        mt_real dist2 = csys->distance2(ref, p);

        if (dist2 < len2) {
          len2 = dist2;
          idx = v[i];
        }
      }
    }
  }
}

bool kdtree::closest_vertex(const vec3r ref, int & idx, vec3r & closest, mt_real & len2) const
{
  if(verts.size() == 0) {
    fprintf(stderr, "%s error: This function needs the kdtree to be filled with\n"
                    "vertices (not triangles)! Aborting!\n", __func__);
    exit(1);
  }

  #if 0
  struct timeval t1, t2;
  gettimeofday(&t1, NULL);
  #endif

  int leaf_idx = -1;
  len2 = FLT_MAX;

  bool is_inside = find_enclosing_leaf(0, ref, leaf_idx);

  if(is_inside) {
    closest_vertex(leaf_idx, ref, idx, len2);

    vec3r cur = _xyz[idx];
    mt_real rad = (ref - cur).length();
    bbox search_region = csys->bbox_from_sphere(ref, rad);

    // we search upwards for a node that encloses the search region,
    // if we find one other than our current leaf node, we have to redo
    // the search from there
    int enc_idx = enclosing_node_upwards(leaf_idx, search_region);
    if(enc_idx != leaf_idx)
      closest_vertex_in_region(0, ref, search_region, idx, len2);
  }
  else {
    closest_vertex(0, ref, idx, len2);
  }

  closest = _xyz[idx];

  #if 0
  gettimeofday(&t2, NULL);
  printf("KDtree vertex lookup in %g sec\n", timediff_sec(t1, t2));
  #endif

  return is_inside;
}

void kdtree::k_closest_vertices(const int k, const vec3r ref, mt_vector<mixed_tuple<mt_real, int> > & vtx) const
{
  if(verts.size() == 0) {
    fprintf(stderr, "%s error: This function needs the kdtree to be filled with\n"
        "vertices (not triangles)! Aborting!\n", __func__);
    exit(EXIT_FAILURE);
  }

  if(size_t(k) > _xyz.size()) {
    fprintf(stderr, "%s error: k is larger than number of vertices in kdtree: k = %d ! Aborting!\n",
            __func__, k);
    exit(EXIT_FAILURE);
  }

  #if 0
  struct timeval t1, t2;
  gettimeofday(&t1, NULL);
  #endif

    // we set the current size to 0, since we want to start populating vtx anew in each call
  vtx.resize(0);
  // we reserve k*2 many items to have good performance when using vtx.push_back()
  vtx.reserve(k*2);

  // we first look at the vertices of the leaf that includes the reference point.
  int leaf_idx = -1;
  bool is_inside = find_enclosing_leaf(0, ref, leaf_idx);
  if(is_inside == false) {
    // should the reference point be outside the bounding box, we first search the
    // closest vertex and then pick the leaf of that vertex
    int closest_vert_idx;
    mt_real closest_vert_len2;
    closest_vertex(0, ref, closest_vert_idx, closest_vert_len2);
    find_enclosing_leaf(0, _xyz[closest_vert_idx], leaf_idx);
  }

  // we take the k closest vertices in this leaf. this gives us a good starting point for
  // our search volume
  for(const int v : *verts[leaf_idx])
    vtx.push_back({csys->distance2(_xyz[v], ref), v});

  std::sort(vtx.begin(), vtx.end());
  if(vtx.size() > size_t(k)) vtx.resize(k);

  // the starting radius is the distance of the furthest away point from vtx.
  // In case we don't have enough vertices we take the upper node bbox
  mt_real current_rad;

  if(vtx.size() > size_t(k * 0.5))
   current_rad = std::sqrt(vtx[vtx.size()-1].v1);
  else {
   vec3r ctr;
   csys->sphere_from_bbox(boxes[nodes[leaf_idx].up], current_rad, ctr);
  }

  bbox search_region = csys->bbox_from_sphere(ref, current_rad);

  // if we havent found enough vertices in leaf node, or the search region is not a
  // subbox of our leaf (i.e. it also intersects other leafs), we need to do a volume search
  if(vtx.size() < size_t(k) || csys->is_subbox(search_region, boxes[leaf_idx]) == false)
  {
    int start_idx = enclosing_node_upwards(leaf_idx, search_region);
    do {
      // we reset vtx so that we get unique values
      vtx.resize(0);

      vertices_in_sphere(start_idx, ref, current_rad, search_region, vtx);

      if(vtx.size() < size_t(k)) {
        current_rad *= 2.0;
        search_region = csys->bbox_from_sphere(ref, current_rad);
        start_idx = enclosing_node_upwards(start_idx, search_region);
      }
    } while(vtx.size() < size_t(k));

    std::sort(vtx.begin(), vtx.end());
    vtx.resize(k);
  }

  if(vtx.size() != (size_t)(k)) {
   fprintf(stderr, "ERROR in %s: couldn't find %d neighbors! The final %zd neighbors are: \n",
           __func__, k, vtx.size());

   for(size_t i=0; i < vtx.size(); i++)
    fprintf(stderr, "(%g, %d) ", vtx[i].v1, vtx[i].v2);

   fprintf(stderr, "\n");
   exit(EXIT_FAILURE);
  }

  #if 0
  gettimeofday(&t2, NULL);
  printf("KDtree vertex lookup in %g sec\n", timediff_sec(t1, t2));
  #endif
}

void kdtree::vertices_in_sphere(const int nidx, const vec3r ref, const mt_real rad,
                                const bbox & search_region, mt_vector<mixed_tuple<mt_real,int> > & vtx) const
{
  const bbox & box       = boxes[nidx];
  const kdtree::node & n = nodes[nidx];

  if(csys->bboxes_intersect(search_region, box)) {
    if (verts[nidx] == NULL) {
      // current node is not a leaf yet -> recurse down both sides
      vertices_in_sphere(n.left,  ref, rad, search_region, vtx);
      vertices_in_sphere(n.right, ref, rad, search_region, vtx);
    }
    else {
      // current node is a leaf -> check all verts
      const mt_vector<int> & v = *(verts[nidx]);
      int num_verts = v.size();

      for (int i = 0; i < num_verts; i++) {
        const vec3r p = _xyz[v[i]];
        mt_real dist2 = csys->distance2(ref, p);

        if (dist2 <= rad * rad)
          vtx.push_back({dist2, v[i]});
      }
    }
  }
}

void kdtree::vertices_in_sphere(const vec3r ref, const mt_real rad, mt_vector<mixed_tuple<mt_real,int> > & vtx) const
{
  if(verts.size() == 0) {
    fprintf(stderr, "%s error: This function needs the kdtree to be filled with\n"
                    "vertices (not triangles)! Aborting!\n", __func__);
    exit(1);
  }

  // we set the current size to 0, since we want to start populating vtx anew in each call
  vtx.resize(0);
  // for performance reasons, we dont want to start at size 1
  vtx.reserve(128);

  bbox search_region = csys->bbox_from_sphere(ref, rad);
  vertices_in_sphere(0, ref, rad, search_region, vtx);

  std::sort(vtx.begin(), vtx.end());
}

void kdtree::count_ray_intersects(const int nod_idx, const ray & r, const mt_real scale_min,
                                  const mt_real scale_max, int & count) const
{
  const bbox & box       = boxes[nod_idx];
  const kdtree::node & n = nodes[nod_idx];

  if (bbox_is_hit_3(r, box, scale_min, scale_max)) {
    if (tris[nod_idx] == NULL) {
      // current node is not a leaf yet -> recurse down both sides
      count_ray_intersects(n.left,  r, scale_min, scale_max, count);
      count_ray_intersects(n.right, r, scale_min, scale_max, count);
    }
    else
    {
      // current node is a leaf -> check all polys
      const mt_vector<tri_elem> & tr = *(tris[nod_idx]);
      int num_tri = tr.size();

      for (int i = 0; i < num_tri; i++)
      {
        double u, v, scale;
        bool did_hit = ray_triangle_intersect(r, tr[i], scale, u, v);
        if (did_hit && scale > scale_min && scale < scale_max)
          count++;
      }
    }
  }
}

int kdtree::count_intersects(const ray & r, const mt_real scale_min, const mt_real scale_max) const
{
  int cnt = 0;
  count_ray_intersects(0, r, scale_min, scale_max, cnt);
  return cnt;
}

int kdtree::enclosing_node_upwards(const int start_idx, const bbox & ref_box) const
{
  int curidx = start_idx;
  if(curidx < 0) curidx = 0;

  const bbox* testbox = &boxes[curidx];

  while(curidx > 0 && csys->is_subbox(ref_box, *testbox) == false) {
    curidx = nodes[curidx].up;
    if(curidx < 0) curidx = 0;

    testbox = &boxes[curidx];
  }

  return curidx;
}


int kdtree::enclosing_node_downwards(const int start_idx, const bbox & ref_box) const
{
  int curidx = start_idx, nextidx = start_idx;

  do {
    curidx = nextidx;
    // check if we have children
    if(nodes[curidx].left > 0) {
      if(csys->is_subbox(ref_box, boxes[nodes[curidx].left]))        // check left
        nextidx = nodes[curidx].left;
      else if(csys->is_subbox(ref_box, boxes[nodes[curidx].right]))  // check right
        nextidx = nodes[curidx].right;
      else
        nextidx = curidx;
    }
    else
      nextidx = curidx;
  } while(nextidx != curidx);

  return curidx;
}

vec3r random_point_in_unit_sphere_fast()
{

  vec3r v;
  do {
    v.x = drand48() * 2. - 1.;
    v.y = drand48() * 2. - 1.;
    v.z = drand48() * 2. - 1.;
  } while(v.length2() > 1.0);

  return v;
}

bool inside_closed_surface(const kdtree & tree, const vec3r pos)
{
  mt_real testcnt = 0.0;
  const int init_tests = 3, full_tests = 8;

  for(int i=0; i < init_tests; i++) {
    int num_intersect = tree.count_intersects(ray(pos, random_point_in_unit_sphere_fast()), 0.0, FLT_MAX);
    if(num_intersect % 2) testcnt += 1.0;
  }

  if((testcnt != mt_real(init_tests)) && (testcnt != 0)) {
    // the decision is not clear, we sample more and decide statistically
    const int num_tests = full_tests - init_tests;

    for(int i=0; i<num_tests; i++) {
      int num_intersect = tree.count_intersects(ray(pos, random_point_in_unit_sphere_fast()), 0.0, FLT_MAX);
      if(num_intersect % 2) testcnt += 1.0;
    }

    testcnt /= mt_real(full_tests);
    return testcnt > 0.5f;
  }
  else {
    // the decision is clear
    return testcnt == mt_real(init_tests);
  }
}


