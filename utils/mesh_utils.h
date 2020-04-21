/**
* @file mesh_utils.h
* @brief Main mesh utilities header.
*
* This file includes the IO utils.
*
* @authors Aurel Neic
* @version
* @date 2016-12-13
*/



#ifndef _MESH_UTILS
#define _MESH_UTILS

#include "geometry_utils.h"
#include "kdtree.h"

/**
* @brief Reduce a vector to all values not matching a certain flag.
*
* The vector is reduced in-place.
*
* @tparam T            Any type with a meaningful "!=" operator.
* @param vec [in,out]  Vector to reduce.
* @param fl  [in]      Flag value.
*/
template<class T>
void reduce_by_flag(mt_vector<T> &vec, T fl)
{
    size_t size = vec.size();
    size_t idx_read, idx_write;

    for(idx_read=0, idx_write=0; idx_read < size; idx_read++)
    {
        if(vec[idx_read] != fl)
        {
            vec[idx_write] = vec[idx_read];
            idx_write++;
        }
    }
    vec.resize(idx_write);
}


/**
* @brief Re-index the vertices of a mesh.
*
* @param mesh     The mesh.
* @param nod_out  The new-to-old index map.
*/
void reindex_nodes(mt_meshdata & mesh, mt_vector<mt_int> & nod_out, bool verbose = true);

void reindex_nodes(mt_meshdata & mesh, bool verbose = true);

/**
* @brief Restrict the elements of a mesh.
*
* @param keep     [in]     Vector defining which element to keep.
* @param mesh     [in,out] The mesh.
* @param eidx_out [out]    The mapping vector from the current element indexing to the former one.
* @param verbose  [in]     Whether to print out progress.
*/
void restrict_elements(const mt_vector<bool> & keep,
                       mt_meshdata & mesh,
                       mt_vector<mt_int> & eidx_out,
                       bool verbose = true);

/**
* @brief Restrict a mesh to a certain subset of elements.
*
* We keep the element with the index "eidx" if and only if
* keep[eidx] == true.
*
* The mesh data is restricted in-place. The nodes are re-indexed.
*
* @param keep     [in]     Vector defining which element to keep.
* @param mesh     [in,out] The mesh.
* @param nod_out  [out]    The mapping vector from the current node indexing to the former one.
* @param eidx_out [out]    The mapping vector from the current element indexing to the former one.
* @param verbose  [in]     Whether to print out progress.
*/
inline
void restrict_meshdata(const mt_vector<bool> & keep,
                       mt_meshdata & mesh,
                       mt_vector<mt_int> & nod_out,
                       mt_vector<mt_int> & eidx_out,
                       bool verbose = true)
{
  restrict_elements(keep, mesh, eidx_out, verbose);
  reindex_nodes(mesh, nod_out, verbose);
}

/**
* @brief Restrict a mesh to a certain subset of elements.
*
* We keep the element with the index "eidx" if and only if
* keep[eidx] == true.
*
* The mesh data is restricted in-place.
*
* @param keep     [in]     Vector defining which element to keep.
* @param mesh     [in,out] The mesh.
* @param eidx_out [out]    The mapping vector from the current element indexing to the former one.
* @param verbose  [in]     Whether to print out progress.
*/
inline
void restrict_meshdata(const mt_vector<bool> & keep,
                       mt_meshdata & mesh,
                       mt_vector<mt_int> & eidx_out,
                       bool verbose = true)
{
  restrict_elements(keep, mesh, eidx_out, verbose);
}

/**
* @brief Extract a meshgraph from a mesh based on tags.
*
* @param tags  [in]  Tags defining what to extract.
* @param mesh  [in]  The mesh.
* @param graph [out] The extracted meshgraph.
*/
void extract_tagged_meshgraph(const MT_USET<mt_int> & tags,
                              const mt_meshdata & mesh,
                              mt_meshgraph & graph);

/**
* @brief Extract a meshgraph from a mesh based on element indices.
*
* @param eidx  [in]  Element indices defining what to extract.
* @param mesh  [in]  The mesh.
* @param graph [out] The extracted meshgraph.
*/
void extract_meshgraph(const mt_vector<mt_int> & eidx,
                       const mt_meshdata & mesh,
                       mt_meshgraph & graph);

/**
* @brief Extract a mesh from based on element indices.
*
* @param eidx  [in]  Element indices defining what to extract.
* @param mesh  [in]  The mesh.
* @param graph [out] The extracted mesh.
*/
void extract_mesh(const mt_vector<mt_int> & eidx,
                       const mt_meshdata & mesh,
                       mt_meshdata & outmesh);


/**
* @brief Insert new element tags into a mesh
*
* @param [out] mesh   The mesh we insert into.
* @param [in]  etags  The tags to insert.
* @param [in]  eidx   The indices where to insert the tags.
*/
void insert_etags(struct mt_meshdata & mesh, const mt_vector<mt_int> & etags, const mt_vector<mt_int> & eidx);

/**
* @brief Insert new element fibers into a mesh
*
* @param [out] mesh   The mesh we insert into.
* @param [in]  lon    The fibers to insert.
* @param [in]  eidx   The indices where to insert the fibers.
*/
void insert_fibers(struct mt_meshdata & mesh, const mt_vector<mt_real> & lon, const mt_vector<mt_int> & eidx);

/**
* @brief Insert new vertex coordinates into a mesh
*
* @param [out] mesh The mesh we insert into.
* @param [in]  xyz  The new coordinates.
* @param [in]  nod  The indices where to insert the coordinates.
*/
void insert_points(struct mt_meshdata & mesh, const mt_vector<mt_real> & xyz, const mt_vector<mt_int> & nod);

/// sort the "in" tuple into the "out" tuple
template<class T>
void sortTuple(const T in1, const T in2, T & out1, T & out2)
{
   if(in1 < in2) out1 = in1, out2 = in2;
   else          out1 = in2, out2 = in1;
}

/// sort the "in" triple into the "out" triple
template<class T>
void sortTriple(const T in1, const T in2, const T in3, T & out1, T & out2, T & out3)
{
    bool t12 = in1 < in2;
  bool t13 = in1 < in3;
  bool t23 = in2 < in3;

  if(t12) {
    //(123),(312),(132)
    if(t23) {
      //123
      out1=in1; out2=in2; out3=in3;
    }
    else {
      //(312),(132)
      if(t13) {
        //132
        out1=in1; out2=in3; out3=in2;
      }
      else {
        //312
        out1=in3; out2=in1; out3=in2;
      }
    }
  }
  else {
    //(213),(231),(321)
    if(t23) {
      //(213),(231)
      if(t13) {
        //213
        out1=in2; out2=in1; out3=in3;
      }
      else {
        //231
        out1=in2; out2=in3; out3=in1;
      }
    }
    else {
      //321
      out1=in3; out2=in2; out3=in1;
    }
  }
}



/**
* @brief Insert the triangle formed by indices (n1, n2, n3)
* into a map.
*
* @param n1       First triangle index.
* @param n2       Second triangle index.
* @param n3       Third triangle index.
* @param eidx     Index of the associated element.
* @param surf     Key buffer.
* @param sele     Value buffer.
* @param surfmap  The map we insert into.
*/
void insert_surf_tri(mt_int n1, mt_int n2, mt_int n3, size_t eidx,
                     triple<mt_int> & surf, tri_sele & sele,
                     MT_MAP<triple<mt_int>, tri_sele> & surfmap);


/**
* @brief Insert the quad formed by indices (n1, n2, n3, n4)
* into a map.
*
* @param n1       First quad index.
* @param n2       Second quad index.
* @param n3       Third quad index.
* @param n4       Forth quad index.
* @param eidx     Index of the associated element.
* @param buff     Sorting buffer.
* @param surf     Key buffer.
* @param sele     Value buffer.
* @param surfmap  The map we insert into.
*/
void insert_surf_quad(mt_int n1, mt_int n2, mt_int n3, mt_int n4, size_t eidx,
                      mt_vector<mt_int> & buff,
                      quadruple<mt_int> & surf,
                      quad_sele & sele,
                      MT_MAP<quadruple<mt_int>, quad_sele> & surfmap);

/// Insert all surfaces of a tetrahedral element into a map.
void insert_surf_tet(const mt_int* nod,
                     const size_t eidx,
                     MT_MAP<triple<mt_int>, tri_sele> & surfmap);

/// Insert all surfaces of a pyramid element into a map.
void insert_surf_pyr(const mt_int* nod, const size_t eidx, mt_vector<mt_int> & buff,
                     MT_MAP<triple<mt_int>, tri_sele> & surfmap,
                     MT_MAP<quadruple<mt_int>, quad_sele> & qsurfmap);

/// Insert all surfaces of a prism element into a map.
void insert_surf_pri(const mt_int* nod, const size_t eidx, mt_vector<mt_int> & buff,
                     MT_MAP<triple<mt_int>, tri_sele> & surfmap,
                     MT_MAP<quadruple<mt_int>, quad_sele> & qsurfmap);

/// Insert all surfaces of a hexahedral element into a map.
void insert_surf_hex(const mt_int* nod, const size_t eidx, mt_vector<mt_int> & buff,
                     MT_MAP<quadruple<mt_int>, quad_sele> & surfmap);


/**
* @brief Compute the surface of a mesh graph.
*
* @param [in]  etype     Element type identifier.
* @param [in]  e2n_cnt   Counts of the mesh graph.
* @param [in]  e2n_con   Connectivities of the mesh graph.
* @param [out] surfmesh  The surface mesh
* @param [in]  full_mesh Whether to generate the full mesh data (cnt, type)
*/
void compute_surface(const mt_vector<elem_t> & etype,
                     const mt_vector<mt_int> & e2n_cnt,
                     const mt_vector<mt_int> & e2n_con,
                     mt_meshdata & surfmesh,
                     const bool full_mesh = true);

/**
* @brief Compute the surface of a mesh graph.
*
* @param [in]  etype     Element type identifier.
* @param [in]  e2n_cnt   Counts of the mesh graph.
* @param [in]  e2n_con   Connectivities of the mesh graph.
* @param [out] tri_surf  Triangle surface elements.
* @param [out] quad_surf Quat surface elements.
*/
void compute_surface(const mt_vector<elem_t> & etype,
                     const mt_vector<mt_int> & e2n_cnt,
                     const mt_vector<mt_int> & e2n_con,
                     MT_MAP<triple<mt_int>, tri_sele>     & tri_surf,
                     MT_MAP<quadruple<mt_int>, quad_sele> & quad_surf);

/**
* @brief Compute the surface of a mesh graph.
*
* @param [in]  etype     Element type identifier.
* @param [in]  e2n_cnt   Counts of the mesh graph.
* @param [in]  e2n_con   Connectivities of the mesh graph.
* @param [in]  ref_eidx The element indices of the original mesh.
* @param [out] tri_surf  Triangle surface elements.
* @param [out] quad_surf Quat surface elements.
*/
void compute_surface(const mt_vector<elem_t> & etype,
                     const mt_vector<mt_int> & e2n_cnt,
                     const mt_vector<mt_int> & e2n_con,
                     const mt_vector<mt_int> & ref_eidx,
                     MT_MAP<triple<mt_int>, tri_sele>     & tri_surf,
                     MT_MAP<quadruple<mt_int>, quad_sele> & quad_surf);


/**
* @brief Compute surface of tag regions
*
* This version does not require the user to compute a meshgraph first when restricting
* the mesh to a subregion.
*
* @param [in]  mesh       The mesh.
* @param [in]  tags       The selected tag regions
* @param [out] tri_surf   Tri surface elements
* @param [out] quad_surf  Quad surface elements
*/
void compute_surface(const mt_meshdata & mesh,
                     const MT_USET<mt_int> & tags,
                     MT_MAP<triple<mt_int>, tri_sele>     & tri_surf,
                     MT_MAP<quadruple<mt_int>, quad_sele> & quad_surf);

/**
* @brief Compute surface of tag regions
*
* This version does not require the user to compute a meshgraph first when restricting
* the mesh to a subregion.
*
* @param [in]  mesh      The mesh.
* @param [in]  tags      The selected tag regions
* @param [out] surfmesh  The surface mesh
* @param [in]  full_mesh Whether to generate the full mesh data (cnt, type)
*/
void compute_surface(const mt_meshdata & mesh,
                     const MT_USET<mt_int> & tags,
                     mt_meshdata & surfmesh,
                     const bool full_mesh = true);

/**
* @brief Convert a surface of triangle and quad maps to a triangle surface.
*
* @param tri         [in]  Triangle map.
* @param quad        [in]  Quad map.
* @param surf_con    [out] Vector of triangle indices.
* @param elem_origin [out] Vector with index of the associated volumetric element.
*/
void surfmap_to_vector(const MT_MAP<struct triple<mt_int>, struct tri_sele> & tri,
                       const MT_MAP<struct quadruple<mt_int>, struct quad_sele> & quad,
                       mt_vector<mt_int> & surf_con,
                       mt_vector<mt_int> & elem_origin);

/**
* @brief Convert a surface from tri / quad maps to mt_meshdata
*
* @param tri         [in]   Triangle map.
* @param quad        [in]   Quad map.
* @param surfmesh    [out]  The surface mesh.
* @param elem_origin [out]  The volumetric element index the surface element belongs to.
*/
void surfmap_to_surfmesh(const MT_MAP<struct triple<mt_int>, struct tri_sele> & tri,
                         const MT_MAP<struct quadruple<mt_int>, struct quad_sele> & quad,
                         mt_meshdata & surf_mesh,
                         mt_vector<mt_int> & elem_origin);

/**
* @brief Add triangle surface elements from vector into a surface map.
*
* @param surf_con  Connectivity vector.
* @param trimap    Surface map.
*/
void vector_to_surfmap(const mt_vector<mt_int> & surf_con,
                       MT_MAP<struct triple<mt_int>, struct tri_sele> & trimap);

/**
* @brief Add surface mesh elements into surface elements maps.
*
* @param surfmesh  The surface mesh.
* @param eidx      Vector with mesh element indices associated to each surface.
* @param tri       The trianlge elements map.
* @param quad      The quad elements map.
*/
void surfmesh_to_surfmap(const mt_meshdata & surfmesh,
                         const mt_vector<mt_int> & eidx,
                         MT_MAP<triple<mt_int>, tri_sele> & tri,
                         MT_MAP<quadruple<mt_int>, quad_sele> & quad);

/**
* @brief Compute a set difference for two surfaces.
*
* The result will be stored in the lhs surface.
*
* @param lhs_tri    Triangle surface elements of lhs.
* @param lhs_quad   Quad surface elements of lhs.
* @param rhs_tri    Triangle surface elements of rhs.
* @param rhs_quad   Quad surface elements of rhs.
*/
void surface_difference(MT_MAP<triple<mt_int>, tri_sele> & lhs_tri,
                        MT_MAP<quadruple<mt_int>, quad_sele> & lhs_quad,
                        const MT_MAP<triple<mt_int>, tri_sele> & rhs_tri,
                        const MT_MAP<quadruple<mt_int>, quad_sele> & rhs_quad);

/**
* @brief Compute a set intersection for two surfaces.
*
* The result will be stored in the lhs surface.
*
* @param lhs_tri    Triangle surface elements of lhs.
* @param lhs_quad   Quad surface elements of lhs.
* @param rhs_tri    Triangle surface elements of rhs.
* @param rhs_quad   Quad surface elements of rhs.
*/
void surface_intersection(MT_MAP<triple<mt_int>, tri_sele> & lhs_tri,
                          MT_MAP<quadruple<mt_int>, quad_sele> & lhs_quad,
                          const MT_MAP<triple<mt_int>, tri_sele> & rhs_tri,
                          const MT_MAP<quadruple<mt_int>, quad_sele> & rhs_quad);

/**
* @brief Compute a set union for two surfaces.
*
* The result will be stored in the lhs surface.
*
* @param lhs_tri    Triangle surface elements of lhs.
* @param lhs_quad   Quad surface elements of lhs.
* @param rhs_tri    Triangle surface elements of rhs.
* @param rhs_quad   Quad surface elements of rhs.
*/
void surface_union(MT_MAP<triple<mt_int>, tri_sele> & lhs_tri,
                   MT_MAP<quadruple<mt_int>, quad_sele> & lhs_quad,
                   const MT_MAP<triple<mt_int>, tri_sele> & rhs_tri,
                   const MT_MAP<quadruple<mt_int>, quad_sele> & rhs_quad);

/**
* @brief Generate the data necessary for writing the *.nbc file of a surface.
*
* @param [in]  mesh         The mesh.
* @param [in]  surf         The surface mesh.
* @param [in]  elem_orig    The mesh elements associated with the surface elements.
* @param [out] nbc          The nbc data.
*/
void generate_nbc_data(const mt_meshdata & mesh,
                       const mt_meshdata & surf,
                       const mt_vector<mt_int> & elem_orig,
                       struct nbc_data & nbc);

/**
* @brief Map data from global indexing to local indexing.
*
*  We map data, which consists of indices in the range (glob[0], ..., glob[N])
*  to the range (0, ..., N). Note that glob needs to be sorted in ascending order.
*
* @param glob   The global indices, indexing the local domain data.
* @param data   The data we need to map from global to local indexing.
*/
template<class T>
void map_glob2loc(const mt_vector<T> & glob, mt_vector<T> & data)
{
  size_t dsize = data.size(), gsize = glob.size();

  mt_vector<T> perm(dsize);
  for(size_t i=0; i<dsize; i++) perm[i] = i;

  binary_sort_copy(data, perm);
  for (size_t i=0, idx=0; i<dsize; i++)
  {
    while (glob[idx] < data[i] && idx < (gsize-1)) idx++;
    data[i] = idx;
  }
  binary_sort_copy(perm, data);
}

/**
* @brief Map data from global indexing to local indexing.
*
*  We map data, which consists of indices in the range (glob[0], ..., glob[N])
*  to the range (0, ..., N). Non-existent indices will be removed.
*
* @param glob   A map between the global indices and the local ones.
* @param data   The data we need to map from global to local indexing.
*/
template<class T>
void map_glob2loc(const MT_MAP<T,T> & g2l, mt_vector<T> & data)
{
  size_t widx=0;
  for(size_t ridx=0; ridx<data.size(); ridx++)
  {
    auto it = g2l.find(data[ridx]);
    if(it != g2l.end())
      data[widx++] = it->second;
  }
  data.resize(widx);
}

/**
* @brief Map data from global indexing to local indexing.
*
*  We map data, which consists of indices in the range (glob[0], ..., glob[N])
*  to the range (0, ..., N). Note that glob needs to be sorted in ascending order.
*
* @param glob   The global indices, indexing the local domain data.
* @param data   The data we need to map from global to local indexing.
*/
template<class T>
void map_connectivity_glob2loc(const MT_MAP<T,T> & g2l,
                               mt_vector<T> & e2n_cnt,
                               mt_vector<T> & e2n_con,
                               mt_vector<bool> & removed)
{
  T widx=0, widx_con=0, ridx_con=0;
  removed.assign(e2n_cnt.size(), false);

  for(size_t ridx=0; ridx < e2n_cnt.size(); ridx++)
  {
    T esize = e2n_cnt[ridx];
    bool good = true;

    for(int i=0; i<esize; i++)
    {
      auto it = g2l.find(e2n_con[ridx_con + i]);
      if(it != g2l.end())
        e2n_con[widx_con + i] = it->second;
      else {
        good = false;
        break;
      }
    }
    ridx_con += esize;

    if(good)
    {
      e2n_cnt[widx] = esize;
      widx++;
      widx_con += esize;
    }
    else
      removed[ridx] = true;
  }
  e2n_cnt.resize(widx);
  e2n_con.resize(widx_con);
}

/**
* @brief Compute union of the surfaces defined by the given list of tag sets.
*
* Surface "i" is computed for the submesh formed by all tags in set tags[i].
*
* @param [in]  mesh       The mesh.
* @param [in]  tags       A sequence of tag sets.
* @param [out] surfmesh   The unified surface mesh.
* @param [out] vtx_set    A set of all vertices in the provided submeshes.
*/
void unified_surface_from_tags(const mt_meshdata & mesh,
                               const mt_vector<MT_USET<mt_int> > & tags,
                               mt_meshdata & surfmesh,
                               MT_USET<mt_int> * vtx_set);

/**
* @brief Compute a unified surface definition from a list of .surf files.
*
* @param list        The list of surface files.
* @param delimiter   The list entry delimiter.
* @param surfmesh    The unified surface.
*/
void unified_surface_from_list(const std::string & list,
                               const char delimiter,
                               mt_meshdata & surfmesh);

/**
* @brief Resize the element data vectors of a mesh.
*
* @param [out] mesh       The mesh to resize.
* @param [in]  size_elem  The new number of elements.
* @param [in]  size_con   The new number of connectivity entries.
*/
void mesh_resize_elemdata(mt_meshdata & mesh, size_t size_elem, size_t size_con);

/**
* @brief Convert surface definition into a full mesh.
*
* @param [in]  mesh       The full mesh.
* @param [in]  surf_con   The surface definition.
* @param [in]  nbc        Neumann BC data.
* @param [out] surfmesh   The surface mesh.
*/
void surface_to_mesh(const mt_meshdata & mesh,
                     const mt_vector<mt_int> & surf_con,
                     const nbc_data & nbc,
                     mt_meshdata & surfmesh);


/**
* @brief Remove faces from a surface connectivity if they contain certain nodes.
*
*/
void remove_nodes_from_surf(const MT_USET<mt_int> & nodes,
                            mt_vector<mt_int> & surf_cnt,
                            mt_vector<mt_int> & surf_con,
                            struct nbc_data & nbc);

void remove_nodes_from_surf(const MT_USET<mt_int> & nodes,
                            mt_vector<mt_int> & surf_cnt,
                            mt_vector<mt_int> & surf_con);

/// remove faces. similar to restrict meshdata, but beter suited for lightweight surfs
void remove_elems_from_surf(const mt_vector<bool> & keep,
                            mt_vector<mt_int> & surf_cnt,
                            mt_vector<mt_int> & surf_con,
                            struct nbc_data & nbc);

/**
* @brief Insert nbc data and mesh vertices into a surface mesh
*/
void nbc_data_into_surface(const mt_meshdata & mesh, const nbc_data & nbc,
                           mt_meshdata & surfmesh);
/**
* @brief Convert element based data to node based data
*
* @tparam V    Data value type. Should support arithmetic operations.
* @param [in]  mesh  The mesh.
* @param [in]  edat  Element data.
* @param [out] ndat  Node data.
*/
template<class V>
void elemData_to_nodeData(const mt_meshdata & mesh,
                          const mt_vector<V> & edat,
                          mt_vector<V> & ndat)
{
  assert(edat.size() == mesh.e2n_cnt.size());
  assert(mesh.n2n_cnt.size() > 0);

  bool vol_weight = mesh.xyz.size() > 0;
  size_t nnodes = mesh.n2n_cnt.size();
  ndat.resize(nnodes, V());

  #ifdef OPENMP
  #pragma omp parallel for schedule(guided)
  #endif
  for(size_t nidx = 0; nidx < nnodes; nidx++)
  {
    if(mesh.n2e_cnt[nidx])
    {
      V avrg = V();
      mt_real total_weight = 0.0;

      mt_int start = mesh.n2e_dsp[nidx], stop = start + mesh.n2e_cnt[nidx];
      for(mt_int i = start; i < stop; i++)
      {
        mt_int eidx = mesh.n2e_con[i];
        const mt_int* con = mesh.e2n_con.data() + mesh.e2n_dsp[eidx];

        mt_real cur_weight = vol_weight ? volume(mesh.etype[eidx], con, mesh.xyz.data()) : 1.0;
        total_weight += cur_weight;
        avrg += edat[eidx] * cur_weight;
      }

      if(total_weight == 0.0)
        printf("%s warning: zero weight!\n", __func__);
      else
        ndat[nidx] = avrg / total_weight;
    }
  }
}

/**
* @brief Convert node based data to element based data
*
* @tparam V    Data value type. Should support arithmetic operations.
* @param [in]  mesh  The mesh.
* @param [in]  ndat  Node data.
* @param [out] edat  Element data.
*/
template<class V>
void nodeData_to_elemData(const mt_meshdata & mesh,
                          const mt_vector<V> & ndat,
                          mt_vector<V> & edat)
{
  assert(ndat.size() == mesh.n2n_cnt.size());
  size_t nelem = mesh.e2n_cnt.size();
  edat.resize(nelem);

  #ifdef OPENMP
  #pragma omp parallel for schedule(guided)
  #endif
  for(size_t eidx = 0; eidx < nelem; eidx++)
  {
    V avrg = V();
    mt_int start = mesh.e2n_dsp[eidx], stop = start + mesh.e2n_cnt[eidx];
    for(mt_int i = start; i < stop; i++)
    {
      mt_int nidx = mesh.e2n_con[i];
      avrg += ndat[nidx];
    }
    edat[eidx] = avrg / mt_real(mesh.e2n_cnt[eidx]);
  }
}

/// compute surface normals in each element
void compute_element_surface_normals(const mt_meshdata & surfmesh,
                                     const mt_vector<mt_real> & xyz,
                                     mt_vector<mt_real> & snrml);

/// compute surface normals in each node
void compute_nodal_surface_normals(const mt_meshdata & surfmesh,
                                   const mt_vector<mt_real> & xyz,
                                   mt_vector<mt_real> & snrml);

/**
* @brief Make sure the surface normals have a desired orientation
*
* @param surf     the surface connectivity
* @param inward   true fore inward facing normals, false for outware facing
*
* @return whether a normals flip was needed
*/
bool apply_normal_orientation(mt_meshdata & surf, const bool inward);

bool apply_normal_orientation(mt_meshdata & surf, kdtree & surf_tree, mt_vector<mt_int> & eidx, const bool inward);

/// remove n2n edges connecting nodal surfaces with a negative scalar product
void remove_bad_surface_edges(const mt_meshdata & mesh,
                              mt_meshdata & surfmesh);

/**
* @brief Convert the data of a PS struct into a mesh.
*
* @param [in]  ps     PS data struct
* @param [out] psmesh PS mesh.
*/
void psdata_to_mesh(const mt_psdata & ps, mt_meshdata & psmesh);


/**
* @brief Generate a unordered map from a coordinate triple to a coordinate index
*
* @param [in]  mesh           The mesh holding the coordinates.
* @param [out] cmap           The unordered map.
* @param [out] scale_factor   The scale factor used when converting the floating point coords to integers.
*/
void generate_coord_map(mt_meshdata & mesh,
                        MT_MAP<triple<mt_int>,mt_int> & cmap,
                        int & scale_factor);

/**
* @brief Unify two meshes into one (mesh1).
*
* @param [in, out] mesh1 First mesh. Holds mesh union after completition.
* @param [in] mesh2 Second mesh.
* @param [in] error_on_empty_intf  Whether an empty interface should yield an error abort.
*/
void mesh_union(mt_meshdata & mesh1, const mt_meshdata & mesh2,
                const bool error_on_empty_intf);


/**
* @brief Function evaluating if nodes are inside a closed surface or not.
*
* @param mesh     The mesh with the nodes we evaluate.
* @param tree     The kdtree of the surface.
* @param in_surf  The vector wheather a node is inside.
*/
void nodes_in_surface(const mt_meshdata & mesh,
                      const kdtree & tree,
                      mt_vector<bool> & in_surf);

/**
* @brief For a given surface kdtree, retag the elements who are determined "inside".
*        There are two criteria for when an element is inside.
*
* @param mesh             The mesh we sample the elements in.
* @param tree             The kdtree of the surface.
* @param newtag           The new tag we set the inside elements to.
* @param sampling_type    The sampling type index. 0 = node based sampling, 1 = element center based sampling.
*/
void sample_elem_tags(mt_meshdata & mesh, mt_mask & tag_found, const kdtree & tree, const mt_int newtag,
                     const int sampling_type);

/**
 * @brief Add an element into a meshdata struct
 * 
 * @param mesh The meshdata struct we add into.
 * @param type The element type to add.
 * @param con  Pointer to the connectivity to add from.
 * @param tag  The element tag to add.
 */
void mesh_add_elem(mt_meshdata & mesh, const elem_t type, const mt_int* con, const mt_int tag);

/**
* @brief Check wheter the elements in a dataset are in a given interval.
*
* Each element must be equal or greater than the interval start and smaller
* than the interval end.
*
* @tparam InputIterator  Iterator type.
* @tparam T              Type that the iterator dereference to.
* @param s               Start of dataset.
* @param e               End of dataset.
* @param i_start         Start of interval.
* @param i_end           End of interval.
*
* @return Wether the elements are in the given interal.
*/
template<typename InputIterator, typename T>
bool set_in_interval(InputIterator s, InputIterator e, const T i_start, const T i_end)
{
  bool in_interv = true;

  while(s != e) {
    T val = *s;
    if(val < i_start || val >= i_end) {
      in_interv = false;
      break;
    }

    ++s;
  }

  return in_interv;
}

void tet_gradient(const mt_int* con, const mt_vector<mt_real> & xyz,
                  const mt_vector<mt_real> & data, mt_real* grad);

void compute_gradient(const mt_meshdata & mesh, const mt_vector<mt_real> & data,
                      const bool nodal_input, const bool nodal_output,
                      mt_vector<vec3r> & grad, mt_vector<mt_real> & mag);

/**
* @brief Read as much surface information as possible from disk. Try to recover the rest.
*
* @param mesh       The mesh associated to the surface. Needed for recovery.
* @param surf       The surface meshdata struct.
* @param eidx       The indices of the volumetric elements associated to the surface elements.
* @param basefile   The path to the surface data on disk.
*/
void read_surf_info(const mt_meshdata & mesh, mt_meshdata & surf, mt_vector<mt_int> & eidx,
                    std::string basefile);
/**
* @brief Write .surf, .vtx and .neubc files for a surface
*
* @param surface  Surface meshdata
* @param nbc      NBC data struct pointer, NULL if .neubc file should not be written
* @param numele   Number of elements in overall mesh
* @param npts     Number of points in overall mesh
* @param writeNBC Whether to write .neubc file.
* @param bname    Base name for all output files.
*/
void write_surf_info(const mt_meshdata & surface, const nbc_data* nbc,
                     const size_t numele, const size_t npts, std::string bname);

#endif

