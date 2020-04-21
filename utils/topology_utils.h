/**
* @file topology_utils.h
* @brief Topology manipulation utility functions.
* @author Aurel Neic
* @version
* @date 2017-02-13
*/


#ifndef _TOPOLOGY_UTILS_H
#define _TOPOLOGY_UTILS_H


#include "mesh_utils.h"
#include "io_utils.h"

#define CORR_THR_DFLT 1000


/**
* @brief Convert a elements set to a nodes set.
*
* This is done by traversing the e2n (elements-to-nodes)
* mesh topology.
*
* @param mesh     Mesh.
* @param nodeSet  Node set.
* @param elemSet  Element set.
*/
void elemSet_to_nodeSet(const mt_meshdata & mesh,
                        const MT_USET<mt_int> & elemSet,
                        MT_USET<mt_int> & nodeSet);

/**
* @brief Convert a node set to a elements set.
*
* This is done by traversing the n2e (node-to-element)
* mesh topology. Therefore, this topology must be set up
* prior to the call.
*
* @param mesh     Mesh.
* @param nodeSet  Node set.
* @param elemSet  Element set.
*/
void nodeSet_to_elemSet(const mt_meshdata & mesh,
                        const MT_USET<mt_int> & nodeSet,
                        MT_USET<mt_int> & elemSet);

/**
* @brief  Functor for the sparse matrix mutliply-transpose operation.
*
* C := (A*B)^T
*
* @tparam T Integer type.
* @tparam S Floating point type.
*/
template<class T, class S>
class sparse_multiply_transpose
{
public:

     /**
     * @brief Execute the matrix mutliply-transpose operation
     *
     */
    void operator()(const mt_vector<T> &_acnt,
                    const mt_vector<T> &_acol,
                    const mt_vector<S> &_aele,
                    const mt_vector<T> &_bcnt,
                    const mt_vector<T> &_bcol,
                    const mt_vector<S> &_bele,
                    const size_t _csize,
                    mt_vector<T> &_ccnt,
                    mt_vector<T> &_ccol,
                    mt_vector<S> &_cele)
    {
        size_t arows = _acnt.size();
        const T *acnt = &_acnt[0];
        const T *acol = &_acol[0];
        const S *aele = &_aele[0];
        const T *bcnt = &_bcnt[0];
        const T *bcol = &_bcol[0];
        const S *bele = &_bele[0];

        size_t brows = _bcnt.size();
        mt_vector<T> _bdsp(brows);
        bucket_sort_offset(_bcnt, _bdsp);

        _ccnt.resize(_csize);
        _ccnt.zero();

        mt_vector<T> _clst(_csize, -1);
        T *ccnt = _ccnt.data();
        T *clst = _clst.data();
        T *bdsp = _bdsp.data();
        size_t csize = 0;
        for(size_t m = 0, i = 0; i < arows; i++)
        {
            for(T j = 0; j < acnt[i]; j++)
            {
                T c = acol[m];
                m++;
                for(T k = bdsp[c]; k < bdsp[c] + bcnt[c]; k++)
                {
                    T d = bcol[k];
                    if(clst[d] != T(i))
                    {
                        clst[d] = i;
                        ccnt[d]++;
                        csize++;
                    }
                }
            }
        }

        _clst.assign(_csize, -1);
        clst = _clst.data();

        mt_vector<T> _cdsp(_csize);
        bucket_sort_offset(_ccnt, _cdsp);

        _ccol.resize(csize);
        _cele.resize(csize);
        T *ccol = _ccol.data();
        S *cele = _cele.data();
        T *cdsp = _cdsp.data();
        for(size_t m = 0, i = 0; i < arows; i++)
        {
            for(T j = 0; j < acnt[i]; j++)
            {
                T c = acol[m];
                S s = aele[m];
                m++;
                for(T k = bdsp[c]; k < bdsp[c] + bcnt[c]; k++)
                {
                    T d = bcol[k];
                    T e = cdsp[d];
                    S t = s * bele[k];
                    if(clst[d] != T(i))
                    {
                        clst[d] = i;
                        ccol[e] = i;
                        cele[e] = t;
                        cdsp[d] = e + 1;
                    }
                    else
                    {
                        cele[e - 1] += t;
                    }
                }
            }
        }
    }
};

/*
 * This a CRS matrix-matrix multiplication routine that
 * does not ouptut a transposed result. hasnt been
 * integrated and tested yet.
 *
 */
template<class T, class S>
void mat_mult_mat_crs(const mt_vector<T> & acnt,
                      const mt_vector<T> & acol,
                      const mt_vector<S> & aele,
                      const mt_vector<T> & bcnt,
                      const mt_vector<T> & bdsp,
                      const mt_vector<T> & bcol,
                      const mt_vector<S> & bele,
                      mt_vector<T> & ccnt,
                      mt_vector<T> & ccol,
                      mt_vector<S> & cele)
{
  T tsize = 0;
  for(size_t i=0, aidx=0; i<acnt.size(); i++)
    for(T j=0; j<acnt[i]; j++, aidx++)
      tsize += bcnt[acol[aidx]];

  mt_vector<T> trow(tsize), tcol(tsize);
  mt_vector<S> tele(tsize);

  T tidx = 0;
  for(size_t i=0, aidx=0; i<acnt.size(); i++) {
    for(T j=0; j<acnt[i]; j++, aidx++) {
      T c = acol[aidx];
      for(T bidx=bdsp[c]; bidx<bdsp[c]+bcnt[c]; bidx++) {
        trow[tidx] = i;
        tcol[tidx] = bcol[bidx];
        tele[tidx] = bele[bidx]*aele[aidx];
        tidx++;
      }
    }
  }

  binary_sort_sort_copy(trow, tcol, tele);
  unique_accumulate(trow, tcol, tele);

  ccnt.resize(acnt.size()); ccnt.zero();
  bucket_sort_count(trow, ccnt);

  ccol.assign(tcol.begin(), tcol.end());
  cele.assign(tele.begin(), tele.end());
}

/**
* @brief Transpose CRS matrix graph A into B.
*
* @tparam T Integer type
* @param a_cnt [in]  A matrix counts
* @param a_con [in]  A matrix column indices
* @param b_cnt [out] B matrix counts
* @param b_con [out] B matrix column indices
*/
template<class T>
void transpose_connectivity(const mt_vector<T> & a_cnt,
                            const mt_vector<T> & a_con,
                            mt_vector<T> & b_cnt,
                            mt_vector<T> & b_con,
                            bool sort_con = true)
{
  if(a_con.size() == 0) return;

  // get the largest node index to determine number of nodes
  T numnodes = *(std::max_element(a_con.begin(), a_con.end())) + 1;
  mt_vector<T> b_row;

  // we compute bcol := arow
  b_con.resize(a_con.size());
  for(size_t i=0, k=0; i < a_cnt.size(); i++)
    for(int j=0; j < a_cnt[i]; j++, k++) b_con[k] = i;

  b_row.assign(a_con.begin(), a_con.end());

  if(sort_con)
    binary_sort_sort(b_row, b_con);
  else
    binary_sort_copy(b_row, b_con);

  b_cnt.resize(numnodes); b_cnt.zero();
  bucket_sort_count(b_row, b_cnt);
}

/**
* @brief Transpose CRS matrix graph A into B.
*
* @tparam T Integer type
* @param a_cnt [in]  A matrix counts
* @param a_con [in]  A matrix column indices
* @param a_ele [in]  A matrix entry values
* @param b_cnt [out] B matrix counts
* @param b_con [out] B matrix column indices
* @param b_ele [out] B matrix entry values
*/
template<class T>
void transpose_connectivity(const mt_vector<T> & a_cnt,
                            const mt_vector<T> & a_con,
                            const mt_vector<T> & a_ele,
                            mt_vector<T> & b_cnt,
                            mt_vector<T> & b_con,
                            mt_vector<T> & b_ele,
                            bool sort = true)
{
  if(a_con.size() == 0) return;

  // get the largest node index to determine number of nodes
  T numnodes = *(std::max_element(a_con.begin(), a_con.end())) + 1;
  mt_vector<T> b_row;

  // we compute bcol := arow
  b_con.resize(a_con.size());
  for(size_t i=0, k=0; i < a_cnt.size(); i++)
    for(int j=0; j < a_cnt[i]; j++, k++) b_con[k] = i;

  b_row = a_con;
  b_ele = a_ele;

  if(sort)
    binary_sort_sort_copy(b_row, b_con, b_ele);
  else
    binary_sort_copy_copy(b_row, b_con, b_ele);

  b_cnt.resize(numnodes); b_cnt.zero();
  bucket_sort_count(b_row, b_cnt);
}

/**
 *  Matrix graph (aka connectivity) multiplication function.
 *
 *  We compute the matrix multiplications C := (A * B)'.
 *  Also the values of the multiplications are returned.
 *
 *  @tparam T Integer type
 *  @param a_cnt [in]  A matrix counts
 *  @param a_con [in]  A matrix column indices
 *  @param b_cnt [in]  B matrix counts
 *  @param b_con [in]  B matrix column indices
 *  @param c_cnt [out] C matrix counts
 *  @param c_con [out] C matrix column indices
 *  @param c_ele [out] C matrix values
 */
template<class T>
void multiply_connectivities(const mt_vector<T> & a_cnt,
                             const mt_vector<T> & a_con,
                             const mt_vector<T> & b_cnt,
                             const mt_vector<T> & b_con,
                             mt_vector<T> & c_cnt,
                             mt_vector<T> & c_con,
                             mt_vector<T> & c_ele)
{
  size_t numnodes = a_cnt.size();

  mt_vector<T> a_ele, b_ele;
  a_ele.assign(a_con.size(), 1);
  b_ele.assign(b_con.size(), 1);

  sparse_multiply_transpose<T, T> multiply;
  multiply(a_cnt, a_con, a_ele, b_cnt, b_con, b_ele, numnodes, c_cnt, c_con, c_ele);
}

/**
 *  Matrix graph (aka connectivity) multiplication function.
 *
 *  We compute the matrix multiplications C := (A * B)'.
 *  Only the matrix graphs are exposed to the outside.
 *
 *  @tparam T Integer type
 *  @param a_cnt [in]  A matrix counts
 *  @param a_con [in]  A matrix column indices
 *  @param b_cnt [in]  B matrix counts
 *  @param b_con [in]  B matrix column indices
 *  @param b_cnt [out] C matrix counts
 *  @param b_con [out] C matrix column indices
 */
template<class T>
void multiply_connectivities(const mt_vector<T> & a_cnt,
                             const mt_vector<T> & a_con,
                             const mt_vector<T> & b_cnt,
                             const mt_vector<T> & b_con,
                             mt_vector<T> & c_cnt,
                             mt_vector<T> & c_con)
{
  mt_vector<T> c_ele;
  multiply_connectivities(a_cnt, a_con, b_cnt, b_con, c_cnt, c_con, c_ele);
}

/**
* @brief Restrict a connectivity graph to edges with a weight equal or above a
*        certain threshold.
*
* @tparam T   Integer type.
* @param [in] cnt  Row size of graph.
* @param [in] con  Connections of graph edges.
* @param [in] wght Weight associated to the edges.
* @param [in] thr  Threshold for restriction.
*
* @post Graph and weights have been restricted.
*/
template<class T>
void restrict_connectivity(mt_vector<T> & cnt,
                           mt_vector<T> & con,
                           mt_vector<T> & wght,
                           T thr)
{
  size_t ridx=0, widx=0;
  for(size_t i=0; i<cnt.size(); i++)
  {
    T c = 0;
    for(T j=0; j<cnt[i]; j++)
    {
      if(wght[ridx] >= thr )
      {
        c++;
        con [widx] = con [ridx];
        wght[widx] = wght[ridx];
        widx++;
      }
      ridx++;
    }
    cnt[i] = c;
  }
  con .resize(widx);
  wght.resize(widx);
}

/**
* @brief Compute full (n2e, n2n) mesh connectivity information.
*
* @param mesh[in, out] The mesh.
*/
void compute_full_mesh_connectivity(mt_meshdata & mesh, bool verbose = true);

/**
* @brief Try to read (n2e, n2n) from disk, if not possible recompute.
*
* @param mesh [in, out] The mesh.
* @param basename [in]  The basename used for storing (n2e, n2n)  
*/
void compute_full_mesh_connectivity(mt_meshdata & mesh, std::string basename, bool verbose = true);

/**
* @brief Insert elements that hold a certain edge into a set.
*
* @param [in]  mesh    The mesh containing the elements.
* @param [in]  v0      The first edge index.
* @param [in]  v1      The second edge index.
* @param [out] elemset The elemset holding the identified elements.
*/
void elements_with_edge(const mt_meshdata & mesh, mt_int v0, mt_int v1, std::set<mt_int> & elemset);

/**
* @brief Insert elements that hold a certain face into a set.
*
* @param [in]  mesh    The mesh containing the elements.
* @param [in]  v0      The first face index.
* @param [in]  v1      The second face index.
* @param [in]  v3      The third face index.
* @param [out] elemset The elemset holding the identified elements.
*/
void elements_with_face(const mt_meshdata & mesh, mt_int v0, mt_int v1, mt_int v2, std::set<mt_int> & elemset);

/**
* @brief Refine a tet element by placing a new node in its center
*
* @param mesh     The mesh.
* @param ref_elem The indices of the elements to refine
*/
void tet_refine_centerPoint(mt_meshdata & mesh, const std::set<mt_int> & ref_elem);


/**
* @brief Execute uniform refinement on Tri and Tetra mesh elements.
*
* @param mesh      The mesh we refine.
* @param edge_map  The edge definition, maps from an edge (tuple) to an edge index.
*/
void refine_uniform(mt_meshdata & mesh,
                    const MT_MAP<tuple<mt_int>, mt_int> & edge_map);

/**
* @brief Add one edge into the edge datastructs.
*
* @param v1            Start edge index.
* @param v2            End edge index.
* @param tupbuff       Buffer tuple used for accessing the edge_map
* @param edgeidx       The global edge index counter.
* @param ele2edge_con  Pointer into the ele2edge connectivity.
* @param edge_map      The edge definition, maps from an edge (tuple) to an edge index.
*/
void add_edge(const mt_int v1,
              const mt_int v2,
              tuple<mt_int> & tupbuff,
              size_t & edgeidx,
              mt_int* ele2edge_con,
              MT_MAP<tuple<mt_int>, mt_int> & edge_map);

/**
* @brief Add one face into the face datastructs.
*
* @param v1            First face index.
* @param v2            Second face index.
* @param v3            Third face index.
* @param tripbuff      Buffer triple used for accessing the face_map
* @param faceidx       The global face index counter.
* @param ele2face_con  Pointer into the ele2face connectivity.
* @param face_map      The global list of faces.
*/
void add_face(const mt_int v1,
              const mt_int v2,
              const mt_int v3,
              triple<mt_int> & tripbuff,
              size_t & faceidx,
              mt_int* ele2face_con,
              MT_MAP<triple<mt_int>, mt_int> & face_map);

/// get the index of one edge
mt_int get_edge(const mt_int v1,
                const mt_int v2,
                tuple<mt_int> & tupbuff,
                const MT_MAP<tuple<mt_int>, mt_int> & edge_map);

/**
* @brief Add edges of a line into the edge datastructs
*
* @param mesh_con      Pointer to start of current element connectivity.
* @param edgeidx       Edge index counter.
* @param ele2edge_con  Pointer of elements-to-edges map connectivity at current element.
* @param edge_map      Map between edge definitions and their indices.
*/
/**
* @brief Add edges of a line into the edge datastructs
*
* @param mesh_con      Pointer to start of current element connectivity.
* @param edgeidx       Edge index counter.
* @param ele2edge_con  Pointer of elements-to-edges map connectivity at current element.
* @param edge_map      Map between edge definitions and their indices.
*/
void add_edges_line(const mt_int* mesh_con,
                    mt_int* ele2edge_con,
                    size_t & edgeidx,
                    MT_MAP<tuple<mt_int>, mt_int> & edge_map);


/**
* @brief Add edges of a triangle into the edge datastructs
*
* @param mesh_con      Pointer to start of current element connectivity.
* @param edgeidx       Edge index counter.
* @param ele2edge_con  Pointer of elements-to-edges map connectivity at current element.
* @param edge_map      Map between edge definitions and their indices.
*/
void add_edges_tri(const mt_int* mesh_con,
                   mt_int* ele2edge_con,
                   size_t & edgeidx,
                   MT_MAP<tuple<mt_int>, mt_int> & edge_map);


/**
* @brief Add edges of a tetrahedra into the edge datastructs
*
* @param mesh_con      Pointer to start of current element connectivity.
* @param edgeidx       Edge index counter.
* @param ele2edge_con  Pointer of elements-to-edges map connectivity at current element.
* @param edge_map      Map between edge definitions and their indices.
*/
void add_edges_tet(const mt_int* mesh_con,
                   mt_int* ele2edge_con,
                   size_t & edgeidx,
                   MT_MAP<tuple<mt_int>, mt_int> & edge_map);

/**
* @brief Get the edge indices for the edges of a tri.
*
* @param [in]  mesh_con  Pointer to start of tri connectivity
* @param [in]  edge_map  The edge map holding the edge indices
* @param [out] edge_idx  Array with the edge indices.
* @param [out] edge      Array with the edges.
*/
void get_edges_tri(const mt_int* mesh_con,
                   const MT_MAP<tuple<mt_int>, mt_int> & edge_map,
                   mt_int* edge_idx, tuple<mt_int>* edge);

/**
* @brief Get the edge indices for the edges of a quad.
*
* @param [in]  mesh_con  Pointer to start of tri connectivity
* @param [in]  edge_map  The edge map holding the edge indices
* @param [out] edge_idx  Array with the edge indices.
* @param [out] edge      Array with the edges.
*/
void get_edges_quad(const mt_int* mesh_con,
                    const MT_MAP<tuple<mt_int>, mt_int> & edge_map,
                    mt_int* edge_idx, tuple<mt_int>* edge);

/**
* @brief Get the edge indices for the edges of a tetra.
*
* @param [in]  mesh_con  Pointer to start of tri connectivity
* @param [in]  edge_map  The edge map holding the edge indices
* @param [out] edge_idx  Array with the edge indices.
* @param [out] edge      Array with the edges.
*/
void get_edges_tet(const mt_int* mesh_con,
                   const MT_MAP<tuple<mt_int>, mt_int> & edge_map,
                   mt_int* edge_idx, tuple<mt_int>* edge);

int get_edges(const mt_meshdata & mesh, const MT_MAP<tuple<mt_int>, mt_int> & edge_map,
              const mt_int eidx, mt_vector<mt_int> & edge_indices, mt_vector<tuple<mt_int>> & edges);
/**
* @brief Add faces of a tetrahedra into the face datastructs
*
* @param mesh_con      Pointer to start of current element connectivity.
* @param faceidx       Face index counter.
* @param ele2face_con  Pointer of elements-to-faces map connectivity at current element.
* @param face_map      Map between edge definitions and their indices.
*/
void add_faces_tet(const mt_int* mesh_con,
                   mt_int* ele2face_con,
                   size_t & faceidx,
                   MT_MAP<triple<mt_int>, mt_int> & face_map);

/**
* @brief Set up edge datastructs for a mesh.
*
* @param [in]  mesh  Surface mesh.
* @param [out] ele2edge  The element-to-edges graph.
* @param [out] edge_map  A map of all unique edges.
*/
void compute_edges(const mt_meshdata & mesh,
                   mt_mapping<mt_int> & ele2edge,
                   MT_MAP<tuple<mt_int>, mt_int> & edge_map);

/**
* @brief Set up face datastructs for a mesh.
*
* @param [in]  mesh  Surface mesh.
* @param [out] ele2face  The element-to-faces graph.
* @param [out] face_map  A map of all unique faces.
*/
void compute_faces(const mt_meshdata & mesh,
                   mt_mapping<mt_int> & ele2face,
                   MT_MAP<triple<mt_int>, mt_int> & face_map);


/**
* @brief Convert an edge map into a mesh datastruct
*
* @param edgemap
* @param edgemesh
*/
void edgemesh_from_edgemap(const MT_MAP<tuple<mt_int>, mt_int> & edgemap,
                           mt_meshdata & edgemesh);

/**
* @brief Compute a line mesh from the interfaces between surface manifolds,
* and optionally from sharp edges.
*
* Sharp edges are identified on the surface of the whole mesh.
*
* @param [in]  mesh               The mesh.
* @param [in]  surfmesh           A surface on the mesh.
* @param [in]  edge_ang           A threshold angle for sharp edges.
* @param [in]  edges_on_geom_surf Whether the sharp edges should be computed on the
*                                 geometric surface or on surfmesh.
* @param [out] linemesh           The line mesh.
*/
void compute_line_interfaces(const mt_meshdata & mesh,
                             const mt_meshdata & surfmesh,
                             const mt_real edge_ang,
                             const bool edges_on_geom_surf,
                             mt_meshdata & linemesh);

/**
* @brief Find the vertex closest to a reference point via linear search over all vertices
*
* @param [in]  xyz    The vertices.
* @param [in]  ref_pt The reference point.
* @param [out] idx    Index of the closest vertex.
* @param [out] dist   Distance of the closest vertex.
*/
void linear_search_vtx(const mt_vector<mt_real> & xyz, const mt_point<mt_real> & ref_pt,
                       mt_int & idx, mt_real & dist);

/**
* @brief Find the vertex closest to a reference point via linear search over all vertices
*
* @param [in]  xyz     The vertices.
* @param [in]  nod_set The nodes we iterate over.
* @param [in]  ref_pt  The reference point.
* @param [out] idx     Index of the closest vertex.
* @param [out] dist    Distance of the closest vertex.
*/
void linear_search_vtx(const mt_vector<mt_real> & xyz,
                       const MT_USET<mt_int> & nod_set,
                       const mt_point<mt_real> & ref_pt,
                       mt_int & idx, mt_real & dist);

/**
* @brief From all the neighbours of a vertex, find that vertex which
*        is closest to a given coordinate.
*
* @param [in]  mesh       The mesh we search on.
* @param [in]  ref_pt     The given coordinate.
* @param [in]  nidx       The index of the vertex spanning the neighbourhood.
* @param [out] end_idx    The index of the closest vertex.
* @param [out] end_dist   The distance to ref_pt of the closest vertex.
*/
void neighbour_closest_to_vtx(const mt_meshdata & mesh,
                              const mt_point<mt_real> & ref_pt,
                              const mt_int nidx,
                              mt_int & end_idx,
                              mt_real & end_dist);

/**
* @brief Find the vertex closest to a given coordinate.
*
* @param [in]  mesh       The mesh we search on.
* @param [in]  ref_pt     The given coordinate.
* @param [in]  start_idx  The vertex we start searching on.
* @param [out] end_idx    The index of the closest vertex.
* @param [out] end_dist   The distance to ref_pt of the closest vertex.
*/
void find_closest_vtx(const mt_meshdata & mesh,
                      const mt_point<mt_real> & ref_pt,
                      const mt_int start_idx,
                      mt_int & end_idx,
                      mt_real & end_dist);

/**
* @brief  Compute the correspondance between two meshes.
*
* The correspondance corr[i] is a vertex in mesh2 which has the closest coordinate
* to i, a vertex in mesh1.
*
* @param [in]  mesh1     Primary mesh.
* @param [in]  mesh2     Secondary mesh.
* @param [out] corr      Vertex correspondance.
* @param [out] corr_dist Correspondance distance.
*/
void compute_correspondance(const mt_meshdata & mesh1,
                            const mt_meshdata & mesh2,
                            mt_vector<mt_int> & corr,
                            mt_vector<mt_real> & corr_dist);

/**
* @brief Traverse full connectivity.
*
* The starting indices i, have to be set by setting reached[i] = true
* before launching this funciton.
*
* @param [in]  n2n_cnt The connectivity row counts.
* @param [in]  n2n_dsp The displacements.
* @param [in]  n2n_con The column entries.
* @param [out] reached Whether each node has been reached.
*/
void traverse_nodal_connectivity(const mt_vector<mt_int> & n2n_cnt,
                           const mt_vector<mt_int> & n2n_dsp,
                           const mt_vector<mt_int> & n2n_con,
                           mt_vector<bool> & reached);

/**
* @brief Traverse connectivity taking into acount blocking nodes.
*
* The starting indices i, have to be set by setting reached[i] = true
* before launching this funciton.
*
* @param [in]  n2n_cnt The connectivity row counts.
* @param [in]  n2n_dsp The displacements.
* @param [in]  n2n_con The column entries.
* @param [in]  blocked A set specifying nodes that block traversal.
* @param [out] reached Whether each node has been reached.
*/
void traverse_nodal_connectivity(const mt_vector<mt_int> & n2n_cnt,
                           const mt_vector<mt_int> & n2n_dsp,
                           const mt_vector<mt_int> & n2n_con,
                           const MT_USET<mt_int> & blocked,
                           mt_vector<bool> & reached);

/**
* @brief Traverse connectivity in a certain radius.
*
* The starting indices i, have to be set by setting reached[i] = true
* before launching this funciton.
*
* @param [in]  n2n_cnt The connectivity row counts.
* @param [in]  n2n_dsp The displacements.
* @param [in]  n2n_con The column entries.
* @param [in]  xyz     The vertex coordinates.
* @param [in]  ref     The reference vertex.
* @param [in]  rad     The radius allowed to traverse.
* @param [out] reached Whether each node has been reached.
*/
void traverse_nodal_connectivity(const mt_vector<mt_int> & n2n_cnt,
                           const mt_vector<mt_int> & n2n_dsp,
                           const mt_vector<mt_int> & n2n_con,
                           const mt_vector<mt_real> & xyz,
                           mt_point<mt_real> & ref,
                           const mt_real rad,
                           mt_vector<bool> & reached);

/**
* @brief Class encapsulating the is-inside-sphere check.
*
*  Used for generic traversal and selection
*
*/
class sphere_check
{
  public:
  mt_real squared_rad;
  vec3r   ctr;

  sphere_check(const vec3r inp_ctr, mt_real inp_rad) :
                squared_rad(inp_rad * inp_rad),
                ctr(inp_ctr)
  {}

  inline bool operator()(const vec3r pt) const
  {
    vec3r e = pt - ctr;
    return e.length2() < squared_rad;
  }
};


class circle_check
{
  public:
  mt_real squared_rad;
  vec3r   nrml;
  vec3r   ctr;
  bool    check_depth;

  circle_check(const vec3r inp_ctr, const vec3r inp_nrml, const mt_real inp_rad, bool depth) :
                squared_rad(inp_rad * inp_rad),
                nrml(inp_nrml),
                ctr(inp_ctr),
                check_depth(depth)
  {}

  inline bool operator()(const vec3r pt) const
  {
    vec3r e = pt - ctr;

    if(!check_depth)
      e -= nrml * e.scaProd(nrml);

    return e.length2() < squared_rad;
  }
};

class box_check
{
  public:
  vec3r   p0, x, y, z;
  mt_real size;

  box_check(const vec3r inp_ctr, const vec3r inp_x, const vec3r inp_y, const vec3r inp_z,
            const mt_real inp_size) :
            p0(0,0,0),
            x(inp_x), y(inp_y), z(inp_z),
            size(inp_size)
  {
    p0 = inp_ctr - x * (size * 0.5) - y * (size * 0.5) - z * (size * 0.5);
  }

  inline bool operator()(const vec3r pt) const
  {
    vec3r e = pt - p0;
    mt_real sca = e.scaProd(x);
    bool in_x = sca > 0 && sca < size;
    sca = e.scaProd(y);
    bool in_y = sca > 0 && sca < size;
    sca = e.scaProd(z);
    bool in_z = sca > 0 && sca < size;

    return in_x && in_y && in_z;
  }
};


/**
* @brief Nodal traversal using a generic inside_check class. This can be used
*        for selecting nodes in different geometric objects.
*
* @tparam CHECK     The inside_check class.
*/
template<class CHECK>
void traverse_nodal_connectivity(const mt_vector<mt_int> & n2n_cnt,
                           const mt_vector<mt_int> & n2n_dsp,
                           const mt_vector<mt_int> & n2n_con,
                           const mt_vector<mt_real> & xyz,
                           const CHECK & inside_check,
                           const MT_USET<mt_int> & blocked,
                           mt_vector<bool> & reached)
{
  size_t nnodes = n2n_cnt.size();

  MT_USET<mt_int> cur_reached, next_reached;
  for(size_t i=0; i < nnodes; i++)
    if(reached[i]) cur_reached.insert(i);

  size_t nreached = cur_reached.size(), last_nreached = 0;

  while(last_nreached < nreached)
  {
    last_nreached = nreached;

    for(mt_int & nidx : cur_reached)
    {
      mt_int start = n2n_dsp[nidx], stop = start + n2n_cnt[nidx];
      for(mt_int j = start; j<stop; j++)
      {
        mt_int cidx = n2n_con[j];
        bool inside = inside_check(vec3r(xyz.data() + 3*cidx));

        if( (!reached[cidx]) && (blocked.count(cidx) == 0) && inside ) {
          reached[cidx] = true;
          next_reached.insert(cidx);
          nreached++;
        }
      }
    }

    cur_reached.clear();
    cur_reached.insert(next_reached.begin(), next_reached.end());
    next_reached.clear();
  }
}

/**
* @brief Traverse the connectivity between surface elements and edges. Stop at blocking edges.
*
* @param surfmesh     The surface mesh.
* @param ele2edge     The mapping between surface elements and edges.
* @param block_edges  The edges that block traversal.
* @param reached      The elements we have reached.
*/
void traverse_surfelem_connectivity(const mt_meshdata & surfmesh,
                                    const mt_mapping<mt_int> & ele2edge,
                                    const MT_USET<mt_int> & block_edges,
                                    mt_vector<bool> & reached);

/**
* @brief Traverse the connectivity between surface elements and edges. Stop at angle threshold.
*
* @param surfmesh     The surface mesh.
* @param ele2edge     The mapping between surface elements and edges.
* @param elem_normals The element normal vector. Used to check against ref_nrml.
* @param ref_nrml     The the referenc normal we test against.
* @param thr          The angle threshold indegrees where we stop traversing.
* @param reached      The elements we have reached.
*/
void traverse_surfelem_connectivity(const mt_meshdata & surfmesh,
                                    const mt_mapping<mt_int> & ele2edge,
                                    const mt_vector<mt_real> & elem_normals,
                                    const vec3r ref_nrml,
                                    mt_real thr,
                                    mt_vector<bool> & reached);

/**
* @brief Decompose a mesh into several parts based on connectivity traversal.
*
* As such, each final part consists of those elements connected by edge paths.
*
* @param mesh       The mesh to decompose.
* @param blocked    A set with nodes blocking traversal.
* @param parts      The decomposed mesh parts.
* @param parts_eidx The element indices of the decomposed mesh parts.
*/
void nodal_connectivity_decomposition(const mt_meshdata & mesh,
                                      mt_vector<mt_meshdata> & parts,
                                      mt_vector<mt_vector<mt_int> > & parts_eidx);


void surf_connectivity_decomposition(const mt_meshdata & surfmesh,
                                     const mt_mapping<mt_int> & ele2edge,
                                     const MT_USET<mt_int> & block_edges,
                                     mt_vector<mt_meshdata> & parts,
                                     mt_vector<mt_vector<mt_int> > & parts_eidx);

/**
* @brief Class transporing element or node data across a provided connectivity.
*
* @tparam T Integer type.
* @tparam S Floating point type.
*/
template<class T, class S>
class mt_data_transport
{
  private:
  const mt_vector<T>            & _cnt;
  const mt_vector<T>            & _dsp;
  const mt_vector<T>            & _con;
  const mt_vector<mt_point<S> > & _vtx;

  mt_vector<S>     _traversed_dist;
  mt_vector<short> _upd;
  const short _upd_max;


  /**
  * @brief Search the connectivity information for a source.
  *
  * The direction of the source has to align sufficiently with the
  * provided gradient field.
  *
  * @param grad     The gradient field.
  * @param eidx     The element index we are looking for a source.
  * @param grad_thr The correlation threshold.
  *
  * @return The index of best correlating source, -1 if no source was found.
  */
  T source_above_grad_thr(const mt_point<S> & grad,
                          const T eidx,
                          const S grad_thr)
  {
    T max_idx = -1;
    S max_sca = -1.0;
    const mt_point<S> & rp = _vtx[eidx];
    _upd[eidx]++;

    T start = _dsp[eidx], stop = start + _cnt[eidx];
    for(T i=start; i<stop; i++)
    {
      T idx = _con[i];
      if( _traversed_dist[idx] > S(-1) ) {
        mt_point<S> p = rp - _vtx[idx];
        p.normalize();

        S sca = p.scaProd(grad);

        if(sca > grad_thr && max_sca < sca) {
          max_sca = sca;
          max_idx = idx;
        }
      }
    }
    return max_idx;
  }

  /**
  * @brief Add non-traversed neighbours of an element to the active list.
  *        Then remove the element itself from the active list.
  *
  * @param eidx         Element index.
  * @param dist_thr     Threshold for the traveled distance.
  * @param active_nodes Active nodes set.
  */
  void expand_active_nodes(const T eidx,
                           const S dist_thr,
                           std::set<T> & active_nodes)
  {
    T start = _dsp[eidx], stop = start + _cnt[eidx];
    S base_dist = _traversed_dist[eidx];

    for(T i=start; i<stop; i++)
    {
      T idx = _con[i];
      if( _traversed_dist[idx] < S(0) ) {
        mt_point<S> e = _vtx[eidx] - _vtx[idx];
        if( (base_dist + e.length2()) < dist_thr )
          active_nodes.insert(idx);
      }
    }
    active_nodes.erase(eidx);
  }

  public:
  /// constructor
  mt_data_transport(const mt_vector<T> & cnt,
                    const mt_vector<T> & dsp,
                    const mt_vector<T> & con,
                    const mt_vector<mt_point<S> > & vtx) :
                    _cnt(cnt), _dsp(dsp), _con(con), _vtx(vtx),
                    _traversed_dist(cnt.size(), S(-1)),
                    _upd(cnt.size(), 0),
                    _upd_max(10)
  {
    assert(vtx.size() == cnt.size());
  }

  /**
  * @brief Transport some initial data along a given gradient field for a certain distance.
  *
  * @tparam V        Data type.
  * @param grad      The gradient field.
  * @param grad_thr  Accuracy threshold when comparing gradient alignment. 1 = complete alignment, 0 = orthogonality
  * @param dist_thr  Travel distance.
  * @param data      The data.
  * @param init      The initial data to transport.
  */
  template<class V>
  void operator()(const mt_vector<mt_point<S> > & grad,
                  const S grad_thr,
                  S dist_thr,
                  mt_vector<V> & data,
                  std::set<mixed_tuple<T, V> > & init)
  {
    // we always compare squared lengths
    dist_thr *= dist_thr;
    // small sanity test
    size_t nelems = _cnt.size();
    assert(grad.size() == nelems && data.size() == nelems);

    std::set<T> active_list;
    mt_vector<T> active_list_vec(nelems);

    // initialization
    for(auto it=init.begin(); it != init.end(); ++it)
    {
      // test if node is represented in the surface
      data[it->v1] = it->v2;
      _traversed_dist[it->v1] = S(0);
    }
    // add all connected elements into the active list
    for(auto it=init.begin(); it != init.end(); ++it)
    {
      T start = _dsp[it->v1], stop = start + _cnt[it->v1];
      for(T j=start; j<stop; j++)
      {
        T idx = _con[j];
        if( _traversed_dist[idx] < S(0) )
          active_list.insert(idx);
      }
    }

    // computation loop
    while(active_list.size() > 0)
    {
      active_list_vec.assign(active_list.begin(), active_list.end());
      std::cout << "Active list size: " << active_list_vec.size() << std::endl;

      for(size_t i=0; i<active_list_vec.size(); i++)
      {
        T eidx = active_list_vec[i];
        const mt_point<S> & cg = grad[eidx];
        T vidx = -1;

        // check for a source with a direction correlating with the gradient field
        vidx = source_above_grad_thr(cg, eidx, grad_thr);

        if(vidx > -1)
        {
          data[eidx] = data[vidx];
          mt_point<S> e = _vtx[vidx] - _vtx[eidx];
          _traversed_dist[eidx] = _traversed_dist[vidx] + e.length2();
          expand_active_nodes(eidx, dist_thr, active_list);
        }

        if(_upd[eidx] > _upd_max)
          active_list.erase(eidx);
      }
    }
  }
};

/**
* @brief Change node order in tetrahedral element definitions in order to
* have correct normal vector orientation.
*
* @param mesh The mesh.
*/
void correct_insideOut(mt_meshdata & mesh);
void correct_insideOut_tri(mt_meshdata & mesh, const mt_vector<mt_real> & surf_nrml);

/**
* @brief Correct duplicate elements in element definition
*
* @param mesh The mesh we work on.
*/
void correct_duplicate_elements(mt_meshdata & mesh);

/**
* @brief Correct duplicate vertices in a mesh.
*
* @param mesh The mesh we work on.
*/
void correct_duplicate_vertices(mt_meshdata & mesh, bool verbose = false);

/// check if element is bad volume-wise
bool has_bad_vol(const elem_t type,
                 const mt_int* con,
                 const mt_real* xyz,
                 const mt_real bad_thr);
/// get angles of a tet
void tet_get_angles(const mt_int* con,
                    const mt_vector<mt_real> & xyz,
                    std::set<mt_real> & angles);

/// insert edge lengths of a tet into a set
inline void tet_add_edgelengths(const mt_int* con,
                         const mt_vector<mt_real> & xyz,
                         mt_vector<edgeele> & len)
{
  tuple<mt_int> tupbuff;

  const mt_int v0 = con[0];
  const mt_int v1 = con[1];
  const mt_int v2 = con[2];
  const mt_int v3 = con[3];

  vec3r p0(xyz.data()+v0*3);
  vec3r p1(xyz.data()+v1*3);
  vec3r p2(xyz.data()+v2*3);
  vec3r p3(xyz.data()+v3*3);

  vec3r e01 = p1 - p0;
  vec3r e02 = p2 - p0;
  vec3r e12 = p2 - p1;

  vec3r e03 = p3 - p0;
  vec3r e13 = p3 - p1;
  vec3r e23 = p3 - p2;

  len.resize(6);

  sortTuple<mt_int>(v0, v1, tupbuff.v1, tupbuff.v2);
  len[0] = {tupbuff.v1, tupbuff.v2, 0, e01.length2(), 0.0f};

  sortTuple<mt_int>(v0, v2, tupbuff.v1, tupbuff.v2);
  len[1] = {tupbuff.v1, tupbuff.v2, 1, e02.length2(), 0.0f};

  sortTuple<mt_int>(v1, v2, tupbuff.v1, tupbuff.v2);
  len[2] = {tupbuff.v1, tupbuff.v2, 2, e12.length2(), 0.0f};

  sortTuple<mt_int>(v0, v3, tupbuff.v1, tupbuff.v2);
  len[3] = {tupbuff.v1, tupbuff.v2, 3, e03.length2(), 0.0f};

  sortTuple<mt_int>(v1, v3, tupbuff.v1, tupbuff.v2);
  len[4] = {tupbuff.v1, tupbuff.v2, 4, e13.length2(), 0.0f};

  sortTuple<mt_int>(v2, v3, tupbuff.v1, tupbuff.v2);
  len[5] = {tupbuff.v1, tupbuff.v2, 5, e23.length2(), 0.0f};
}

/// insert edge lengths of a triangle into a set
inline void tri_add_edgelengths(const mt_int* con,
                         const mt_vector<mt_real> & xyz,
                         mt_vector<edgeele> & len)
{
  tuple<mt_int> tupbuff;

  const mt_int v0 = con[0];
  const mt_int v1 = con[1];
  const mt_int v2 = con[2];

  vec3r p0(xyz.data()+v0*3);
  vec3r p1(xyz.data()+v1*3);
  vec3r p2(xyz.data()+v2*3);

  vec3r e01 = p1 - p0;
  vec3r e02 = p2 - p0;
  vec3r e12 = p2 - p1;

  len.resize(3);

  sortTuple<mt_int>(v0, v1, tupbuff.v1, tupbuff.v2);
  len[0] = {tupbuff.v1, tupbuff.v2, 0, e01.length2(), 0.0f};

  sortTuple<mt_int>(v0, v2, tupbuff.v1, tupbuff.v2);
  len[1] = {tupbuff.v1, tupbuff.v2, 1, e02.length2(), 0.0f};

  sortTuple<mt_int>(v1, v2, tupbuff.v1, tupbuff.v2);
  len[2] = {tupbuff.v1, tupbuff.v2, 2, e12.length2(), 0.0f};
}

inline void add_edges_lengths(const mt_meshdata & mesh,
                       const mt_int eidx,
                       mt_vector<edgeele> & len)
{
  switch(mesh.etype[eidx]) {
    case Tri:
      tri_add_edgelengths(mesh.e2n_con.data() + mesh.e2n_dsp[eidx], mesh.xyz, len);
      break;

    case Tetra:
      tet_add_edgelengths(mesh.e2n_con.data() + mesh.e2n_dsp[eidx], mesh.xyz, len);
      break;

    default:
      fprintf(stderr, "mesh merging error: Element type not supported yet.\n");
      exit(1);
  }
}


/**
* @brief Write a graph as a PETSc matrix to disc.
*
* @param cnt       The row counts of the graph.
* @param con       The connectivity of the graph.
* @param filename  The output file name.
*/
void write_graph(const mt_vector<mt_int> & cnt,
                    const mt_vector<mt_int> & con,
                    const char* filename);

/**
* @brief Class for setting up a DD partitioning in shared memory
*/
class mt_shmem_partitioner {

  private:
  mt_vector<MT_USET<mt_int> > part;
  mt_vector<bool>             itf;

  void compute_partitions(const mt_meshdata & mesh,
                          const int npart);

  public:
  mt_shmem_partitioner(const mt_meshdata & mesh, const int npart) {
    compute_partitions(mesh, npart);
  }

  const MT_USET<mt_int> & node_set(int pidx) const
  {
    return part[pidx];
  }

  const mt_vector<bool> & get_interface() const
  {
    return itf;
  }
};

/**
* @brief Find the enclosing element of a given coordinate.
*
* The region where element is searched is defined by the
* given set of nodes. All elements connected to the given
* nodes will be checked.
*
* When testing an element, one also gets the resulting interpolation
* weights for free.
*
* @param [in]  mesh   The mesh we search on.
* @param [in]  nod    The elements connected to this nodes will be checked.
* @param [in]  pt     We search for the element enclosing this coordinate.
* @param [out] elem   The enclosing element.
* @param [out] a      First interpolation weight.
* @param [out] b      Second interpolation weight.
* @param [out] c      Third interpolation weight.
*/
void enclosing_element(const mt_meshdata & mesh,
                       const MT_USET<mt_int> & nod,
                       const mt_point<mt_real> & pt,
                       const mt_real edge_len,
                       mt_int & elem,
                       mt_real & a,
                       mt_real & b,
                       mt_real & c);

/**
* @brief Expand the set of nodes by one element layer.
*/
void expand_nodeset(const mt_meshdata & mesh,
                    MT_USET<mt_int> & nod);

/**
* @brief Expand the node set by all reachable vertices that
* are also in a given radius to a reference point.
*
* @param [in]  mesh  The mesh.
* @param [in]  ref   The reference point.
* @param [in]  squared_rad  The squared radius.
* @param [out] nod   The nodes.
*/
void expand_in_radius(const mt_meshdata & mesh,
                      const mt_point<mt_real> ref,
                      const mt_real squared_rad,
                      std::set<mt_int> & nod);

/**
* @brief Estimate edge average length in a mesh
*
* @param mesh The mesh.
* @param full Whether to look at full mesh or just sample a couple elements.
*
* @return The estimated average edge length.
*/
mt_real avrg_edgelength_estimate(const mt_meshdata & mesh, bool full = false);

/**
* @brief Estimate min edge length
*
* @param mesh The mesh.
* @param full Whether to look at full mesh or just sample a couple elements.
* @param nonzeromin Whether to calculate the nonzero estimate or possible zero estimate
*
* @return The estimated min edge length.
*/
mt_real min_edgelength_estimate(const mt_meshdata & mesh, bool full = false, bool nonzeromin = false);

/**
* @brief Interpolate data from one mesh onto another
*
* @tparam V Data type.
* @param imesh       Mesh of the input data.
* @param omesh       Mesh of the output data.
* @param corr        Nodal correspondance between the two meshes.
* @param corr_dist   Distance between 2 corresponding nodes.
* @param idat        Input data.
* @param odat        Output data.
*/
template<class V>
void nodal_interpolation(const mt_meshdata & imesh,
                        const mt_meshdata & omesh,
                        const mt_manifold & imnfld,
                        const mt_manifold & omnfld,
                        const mt_vector<mt_int> & corr,
                        const mt_vector<mt_real> & corr_dist,
                        const mt_vector<V> & idat,
                        mt_vector<V> & odat)
{
  size_t nnodes = corr.size();
  odat.resize(nnodes);

  int direct = 0, vol_int = 0, cloud_int = 0;
  double avg_edge_len = avrg_edgelength_estimate(omesh);

  PROGRESS<size_t> prg(nnodes, "Interpolation progress: ");

  #ifdef OPENMP
  #pragma omp parallel for schedule(guided) reduction(+:direct,vol_int,cloud_int)
  #endif
  for(size_t nidx = 0; nidx < nnodes; nidx++)
  {
    if(sqrt(corr_dist[nidx]) > (avg_edge_len * 1e-3) ) {
      mt_point<mt_real> ref_pt(omesh.xyz.data() + nidx*3);
      mt_int cidx = corr[nidx];
      MT_USET<mt_int> testnod;
      mt_int eidx = -1;
      mt_real a, b, c;

      // if the corrsponding input vertex is on the surface, we use the surface
      // mesh, else the volumetric mesh
      const mt_meshdata* cmesh = omnfld.on_mnfld.count(nidx) && imnfld.on_mnfld.count(cidx) ?
                                 &imnfld.mesh : &imesh;

      // first we test elements surrounding cidx
      testnod.insert(cidx);
      enclosing_element(*cmesh, testnod, ref_pt, avg_edge_len, eidx, a, b, c);
      // if we dont find an enclosing element we extend to the next neighbours
      if(eidx == -1) {
        expand_nodeset(*cmesh, testnod);
        testnod.erase(cidx);
        enclosing_element(*cmesh, testnod, ref_pt, avg_edge_len, eidx, a, b, c);
      }

      if(eidx > -1) {
        // volume or surface interpolation
        vol_int++;

        const mt_int* con = cmesh->e2n_con.data() + cmesh->e2n_dsp[eidx];
        switch(cmesh->etype[eidx]) {
          case Tri:
            {
              // barycentric interpolation
              V fa = idat[con[0]], fb = idat[con[1]], fc = idat[con[2]];
              odat[nidx] = fa*a + fb*b + fc*c;
              break;
            }
          case Tetra:
            {
              // interpolation based on reference element coordinates
              V fa = idat[con[0]], fb = idat[con[1]],
                fc = idat[con[2]], fd = idat[con[3]];
              odat[nidx] = fa + (fb - fa)*a + (fc - fa)*b + (fd - fa)*c;
              break;
            }
          case Prism:
          {
            //Nodal interpolation based on generalized barycentric coordinates (aka wachspass coordinates)
            mt_vector<V> data(6);
            mt_vector<mt_real> weights(6);
            mt_vector<mt_point<mt_real> > prism_points(6);
            for(size_t i=0; i < 6; i++)
            {
              data[i] = idat[con[i]];
              prism_points[i] = mt_point<mt_real>(cmesh->xyz.data() + con[i]*3);
            }
            //build face normals
            //faces of the prisms have to be oriented counterclockwise.
            //faces are given by {0,2,1}, {3,5,4}, {1,2,4,5}, {2,0,3,4}, {0,1,5,4}
            mt_vector<vec3r> normals(5);
            normals[0] = ((prism_points[2]-prism_points[0]).crossProd(prism_points[1]-prism_points[0])); normals[0].normalize();
            normals[1] = ((prism_points[5]-prism_points[3]).crossProd(prism_points[4]-prism_points[3])); normals[1].normalize();
            normals[2] = ((prism_points[2]-prism_points[1]).crossProd(prism_points[4]-prism_points[1])); normals[2].normalize();
            normals[3] = ((prism_points[0]-prism_points[2]).crossProd(prism_points[3]-prism_points[2])); normals[3].normalize();
            normals[4] = ((prism_points[1]-prism_points[0]).crossProd(prism_points[5]-prism_points[0])); normals[4].normalize();
            mt_vector<mt_vector<mt_int> > vertex2normals(6);
            for(size_t i=0; i < 6; i++)
              vertex2normals[i].resize(3);
            vertex2normals[0][0] = 0;
            vertex2normals[0][1] = 3;
            vertex2normals[0][2] = 4;

            vertex2normals[1][0] = 0;
            vertex2normals[1][1] = 2;
            vertex2normals[1][2] = 4;

            vertex2normals[2][0] = 0;
            vertex2normals[2][1] = 2;
            vertex2normals[2][2] = 3;

            vertex2normals[3][0] = 1;
            vertex2normals[3][1] = 3;
            vertex2normals[3][2] = 4;

            vertex2normals[4][0] = 1;
            vertex2normals[4][1] = 2;
            vertex2normals[4][2] = 3;

            vertex2normals[5][0] = 1;
            vertex2normals[5][1] = 2;
            vertex2normals[5][2] = 4;

            //calculate weights
            mt_real weight_sum = 0.;
            for(size_t i=0; i < 6; i++)
            {
              const mt_real nax = normals[vertex2normals[i][0]].x;
              const mt_real nay = normals[vertex2normals[i][0]].y;
              const mt_real naz = normals[vertex2normals[i][0]].z;
              const mt_real nbx = normals[vertex2normals[i][1]].x;
              const mt_real nby = normals[vertex2normals[i][1]].y;
              const mt_real nbz = normals[vertex2normals[i][1]].z;
              const mt_real ncx = normals[vertex2normals[i][2]].x;
              const mt_real ncy = normals[vertex2normals[i][2]].y;
              const mt_real ncz = normals[vertex2normals[i][2]].z;
              const mt_real kappa = std::abs(-(naz*nby*ncx) + nay*nbz*ncx + naz*nbx*ncy - nax*nbz*ncy - nay*nbx*ncz + nax*nby*ncz);
              mt_point<mt_real> test = prism_points[i] - ref_pt;
              mt_real val = 1.0;
              for(size_t j=0; j < 3; j++)
                val *= normals[vertex2normals[i][j]].scaProd(test);
              weights[i] = kappa / val;
              weight_sum += weights[i];
              odat[nidx] += data[i] * weights[i];
            }
            odat[nidx] /= weight_sum;
            break;
          }
          case Hexa:
          {
            //Nodal interpolation based on generalized barycentric coordinates (aka wachspass coordinates)
            mt_vector<V> data(8);
            mt_vector<mt_real> weights(8);
            mt_vector<mt_point<mt_real> > hexa_points(8);
            for(size_t i=0; i < 8; i++)
            {
              data[i] = idat[con[i]];
              hexa_points[i] = mt_point<mt_real>(cmesh->xyz.data() + con[i]*3);
            }
            //build face normals
            //faces of the hex have to be oriented counterclockwise.
            //faces are given by {0,3,2,1}, {1,2,6,7}, {4,7,6,5}, {0,4,5,3}, {2,3,5,6}, {0,1,7,4}
            mt_vector<vec3r> normals(6);
            normals[0] = ((hexa_points[3]-hexa_points[0]).crossProd(hexa_points[2]-hexa_points[0])); normals[0].normalize();
            normals[1] = ((hexa_points[2]-hexa_points[1]).crossProd(hexa_points[6]-hexa_points[1])); normals[1].normalize();
            normals[2] = ((hexa_points[7]-hexa_points[4]).crossProd(hexa_points[6]-hexa_points[4])); normals[2].normalize();
            normals[3] = ((hexa_points[4]-hexa_points[0]).crossProd(hexa_points[5]-hexa_points[0])); normals[3].normalize();
            normals[4] = ((hexa_points[3]-hexa_points[2]).crossProd(hexa_points[5]-hexa_points[2])); normals[4].normalize();
            normals[5] = ((hexa_points[1]-hexa_points[0]).crossProd(hexa_points[7]-hexa_points[0])); normals[5].normalize();

            mt_vector<mt_vector<mt_int> > vertex2normals(8);
            for(size_t i=0; i < 8; i++)
              vertex2normals[i].resize(3);
            vertex2normals[0][0] = 0;
            vertex2normals[0][1] = 3;
            vertex2normals[0][2] = 5;

            vertex2normals[1][0] = 0;
            vertex2normals[1][1] = 1;
            vertex2normals[1][2] = 5;

            vertex2normals[2][0] = 0;
            vertex2normals[2][1] = 1;
            vertex2normals[2][2] = 4;

            vertex2normals[3][0] = 0;
            vertex2normals[3][1] = 3;
            vertex2normals[3][2] = 4;

            vertex2normals[4][0] = 2;
            vertex2normals[4][1] = 3;
            vertex2normals[4][2] = 5;

            vertex2normals[5][0] = 2;
            vertex2normals[5][1] = 3;
            vertex2normals[5][2] = 4;

            vertex2normals[6][0] = 2;
            vertex2normals[6][1] = 4;
            vertex2normals[6][2] = 1;

            vertex2normals[7][0] = 2;
            vertex2normals[7][1] = 1;
            vertex2normals[7][2] = 5;

            //calculate weights
            mt_real weight_sum = 0.;
            for(size_t i=0; i < 8; i++)
            {
              const mt_real nax = normals[vertex2normals[i][0]].x;
              const mt_real nay = normals[vertex2normals[i][0]].y;
              const mt_real naz = normals[vertex2normals[i][0]].z;
              const mt_real nbx = normals[vertex2normals[i][1]].x;
              const mt_real nby = normals[vertex2normals[i][1]].y;
              const mt_real nbz = normals[vertex2normals[i][1]].z;
              const mt_real ncx = normals[vertex2normals[i][2]].x;
              const mt_real ncy = normals[vertex2normals[i][2]].y;
              const mt_real ncz = normals[vertex2normals[i][2]].z;
              const mt_real kappa = std::abs(-(naz*nby*ncx) + nay*nbz*ncx + naz*nbx*ncy - nax*nbz*ncy - nay*nbx*ncz + nax*nby*ncz);
              mt_point<mt_real> test = hexa_points[i] - ref_pt;
              mt_real val = 1.0;
              for(size_t j=0; j < 3; j++)
                val *= normals[vertex2normals[i][j]].scaProd(test);
              weights[i] = kappa / val;
              weight_sum += weights[i];
              odat[nidx] += data[i] * weights[i];
            }
            odat[nidx] /= weight_sum;
            break;
          }
          default:
            fprintf(stderr, "interpolate error: Unsupported elem type.\n");
            exit(EXIT_FAILURE);
            break;
        }
      }
      else {
        // pointcloud interpolation
        cloud_int++;

        std::set<mt_int> nod;
        {
          nod.insert(cidx);
          mt_real rad = avg_edge_len * 3.0;
          expand_in_radius(imesh, ref_pt, rad*rad, nod);
        }

        // in the nodal_interpolation we always should be able to at least map the
        // value from the correspondance vertex. It may not be the thing to do actually
        if(nod.size() == 0) {
          fprintf(stderr, "%s error: Could not cloud-interpolate vertex %ld !\n",
                  __func__, (long int) nidx);
          odat[nidx] = V();
        }
        else {
          mt_vector<mt_real> coeff(nod.size());
          mt_vector<mt_real> con  (nod.size());
          mt_int widx = 0;
          for(auto n : nod) {
            con[widx] = n;
            mt_point<mt_real> e = mt_point<mt_real>(imesh.xyz.data() + n*3) - ref_pt;
            mt_real l = e.length();
            coeff[widx] = 1.0 / l / l; // squared inverse distance
            widx++;
          }
          // normalize coeffs
          mt_real sum = 0;
          for(const mt_real & l : coeff) sum += l;
          for(mt_real & co : coeff) co /= sum;

          // interpolate value
          odat[nidx] = V();
          for(size_t i=0; i<coeff.size(); i++)
            odat[nidx] += idat[con[i]] * coeff[i];
        }
      }
    }
    else {
      // direct mapping
      direct++;

      odat[nidx] = idat[corr[nidx]];
    }

    #ifdef OPENMP
    #pragma omp critical
    #endif
    {
      prg.next();
    }
  }

  prg.finish();

  printf("Interpolation stats:\n");
  printf("Direct mappings: %d, surface/volumetric interp: %d, point-distance interp: %d\n\n",
         direct, vol_int, cloud_int);
}


/**
* @brief Interpolate element data between two meshes.
*
* @tparam V     Data type that implements '=' operator.
* @param imesh               Input mesh.
* @param omesh               Output mesh.
* @param imesh_vert_tree     KDtree with the vertices of input mesh
* @param idat                Input data.
* @param odat                Output data.
*/
template<class V>
void elementwise_interpolation(const mt_meshdata & imesh,
                        const mt_meshdata & omesh,
                        const kdtree & imesh_vert_tree,
                        const mt_vector<V> & idat,
                        mt_vector<V> & odat)
{
  odat.assign(omesh.e2n_cnt.size(), V());

  // int num_layers = (avrg_edgelength_estimate(omesh) / avrg_edgelength_estimate(imesh)) * 1.5;
  // if(num_layers < 1) num_layers = 1;
  int num_layers = 1;
  int vol_int = 0, cloud_int = 0;
  double avg_edge_len = avrg_edgelength_estimate(omesh);

  PROGRESS<size_t> prg(omesh.e2n_cnt.size(), "Interpolation progress: ");

  #ifdef OPENMP
  #pragma omp parallel for schedule(guided) reduction(+:vol_int,cloud_int)
  #endif
  for(size_t eidx = 0; eidx < omesh.e2n_cnt.size(); eidx++)
  {
    int idx;
    vec3r closest;
    mt_real len2;

    vec3f ectr = barycenter(omesh, eidx);
    imesh_vert_tree.closest_vertex(ectr, idx, closest, len2);

    MT_USET<mt_int> testnod;
    mt_int found_eidx = -1;
    mt_real a, b, c;

    // first we test elements surrounding idx
    testnod.insert(idx);
    enclosing_element(imesh, testnod, ectr, avg_edge_len, found_eidx, a, b, c);
    // if we dont find an enclosing element we extend to the next neighbours
    if(found_eidx == -1) {

      for(int i=0; i<num_layers; i++)
        expand_nodeset(imesh, testnod);

      testnod.erase(idx);
      enclosing_element(imesh, testnod, ectr, avg_edge_len, found_eidx, a, b, c);
    }

    if(found_eidx > -1) {
      vol_int++;
      odat[eidx] = idat[found_eidx];
    }
    else {
      // in case we couldnt find an enclosing element, we search for the element with the
      // closest centerpoint
      cloud_int++;

      testnod.clear();
      testnod.insert(idx);

      MT_USET<mt_int> testelem;
      nodeSet_to_elemSet(imesh, testnod, testelem);

      float min_dist = FLT_MAX;
      mt_int min_idx = -1;

      for(const mt_int & tidx : testelem) {
        vec3f tctr = barycenter(imesh, tidx);
        float l = (tctr - ectr).length();
        if(min_dist > l) {
          min_dist = l;
          min_idx  = tidx;
        }
      }

      if(min_idx > -1)
        odat[eidx] = idat[min_idx];
    }

    #ifdef OPENMP
    #pragma omp critical
    #endif
    {
      prg.next();
    }
  }

  prg.finish();

  printf("Interpolation stats:\n");
  printf("Surface/volumetric interp: %d, Centerpoint-distance interp: %d\n\n", vol_int, cloud_int);
}

/**
* @brief Interpolate data from a pointcloud onto a mesh.
*
* @tparam V Data type.
* @param omesh       Mesh of the output data.
* @param ipts        Input data points.
* @param thr         Optional cutoff threshold.
* @param idat        Input data.
* @param odat        Output data.
*/
template<class V>
void pointcloud_interpolation(const mt_meshdata & omesh,
                              const mt_vector<mt_real> & ipts,
                              const mt_real thr,
                              const mt_real edge_est,
                              const mt_vector<V> & idat,
                              mt_vector<V> & odat)
{
  const mt_real eps = 1e-6 * edge_est;
  size_t npts   = ipts.size() / 3;
  size_t nnodes = omesh.xyz.size() / 3;
  odat.resize(nnodes);

  PROGRESS<size_t> prg(nnodes, "Interpolation progress: ");

  #ifdef OPENMP
  #pragma omp parallel
  #endif
  {
    mt_vector<mt_real> coeff(npts);

    #ifdef OPENMP
    #pragma omp for schedule(guided)
    #endif
    for(size_t nidx = 0; nidx < nnodes; nidx++)
    {
      mt_point<mt_real> ref_pt(omesh.xyz.data() + nidx*3);
      int num_used = 0;

      for(size_t widx=0; widx < npts; widx++)
      {
        mt_point<mt_real> e = mt_point<mt_real>(ipts.data() + widx*3) - ref_pt;
        mt_real l = e.length();

        // in case one vertex is really close to reference, we take only that vertex
        if(l < eps) {
          coeff.zero();
          coeff[widx] = 1.0;
          num_used++;
          break;
        }

        if(l < thr) {
          coeff[widx] = 1.0 / (l * l); // squared inverse distance
          num_used++;
        }
        else
          coeff[widx] = 0.0;
      }

      odat[nidx] = V();
      if(num_used) {
        // normalize coeffs
        mt_real sum = 0;
        for(const mt_real & co : coeff) sum += co;
        for(mt_real & co : coeff) co /= sum;

        // interpolate value
        for(size_t i=0; i<coeff.size(); i++)
          odat[nidx] += idat[i] * coeff[i];
      }

      #ifdef OPENMP
      #pragma omp critical
      #endif
      {
        prg.next();
      }
    }
  }

  prg.finish();
}

/**
* @brief Find local minima of a data field, starting at specified seed locations.
*
* @tparam T     Integer datatype.
* @tparam S     Data field datatype.
* @param cnt    Row counts of connectivity graph.
* @param dsp    Row displacements of connectivity graph.
* @param con    Connectivity entries.
* @param data   Data field.
* @param seeds  Seed vertex indices.
*
* @post Every seeds component contains the minimum found from the
* respective seed location.
*/
template<class T, class S>
void find_min(const mt_vector<T> & cnt,
              const mt_vector<T> & dsp,
              const mt_vector<T> & con,
              const mt_vector<S> & data,
              mt_vector<T> & seeds)
{
  #ifdef OPENMP
  #pragma omp parallel for schedule(guided)
  #endif
  for(size_t sidx=0; sidx < seeds.size(); sidx++)
  {
    T s = seeds[sidx];
    S      minval  = data[s],
           lastval = minval;
    T minidx  = s;

    do {
      lastval = minval;

      T start = dsp[minidx], stop = start + cnt[minidx];
      for(T i=start; i<stop; i++)
      {
        T c = con[i];
        S t = data[c];
        if(minval > t) {
          minval = t;
          minidx = c;
        }
      }
    }
    while(minval < lastval);

    seeds[sidx] = minidx;
  }
}

/**
* @brief Find geodesic paths to local minima of a data field,
*        starting at specified seed locations.
*
* @tparam T     Integer datatype.
* @tparam S     Data field datatype.
*
* @param cnt    Row counts of connectivity graph.
* @param dsp    Row displacements of connectivity graph.
* @param con    Connectivity entries.
* @param data   Data field.
* @param geo    The list geodesic starting locations.
*
* @post geo contains the geodesic paths for each seed location.
*/
template<class T, class S>
void find_min(const mt_vector<T> & cnt,
              const mt_vector<T> & dsp,
              const mt_vector<T> & con,
              const mt_vector<S> & data,
              mt_vector<std::list<T> > & geo)
{
  #ifdef OPENMP
  #pragma omp parallel for schedule(guided)
  #endif
  for(size_t sidx=0; sidx < geo.size(); sidx++)
  {
    typename std::list<T>::iterator s = geo[sidx].begin();
    S      minval  = data[*s],
           lastval = minval;
    T minidx  = *s;

    do {
      lastval = minval;

      T start = dsp[minidx], stop = start + cnt[minidx];
      for(T i=start; i<stop; i++)
      {
        T c = con[i];
        S t = data[c];
        if(minval > t) {
          minval = t;
          minidx = c;
          geo[sidx].push_back(minidx);
        }
      }
    }
    while(minval < lastval);
  }
}

/**
* @brief Compute the stats of the edge lengths of an element.
*
* @param type      Element type.
* @param con       Connectivity pointer pointing to start of element.
* @param xyz       Pointer to global vertex array.
* @param min_edge  Minimum edge length.
* @param max_edge  Maximum edge length.
* @param avg_edge  Average edge length.
*/
void element_edges_stats(const elem_t type, const mt_int* con, const mt_real* xyz,
                           mt_real & min_edge, mt_real & max_edge, mt_real & avg_edge);

void generate_surfmesh_sizing_field(const mt_meshdata & surfmesh, mt_vector<mt_real> & sizes);

/**
* @brief Sharp edges on a surface are identified based on the maximum angle difference between the
* surface elements.
*
* Only the vertices in vtxlist are checked if the are located on a sharp edge.
*
* @param [in]  surfmesh   The surface definition.
* @param [in]  xyz        The surface vertex coordinates.
* @param [in]  vtxlist    The vertices to check.
* @param [out] edge_vtx   A set conatining the vertices on sharp edges.
* @param [in]  thr        The threshold angle (in degrees).
*/
void identify_sharp_edge_nodes(const mt_meshdata & surfmesh,
                               const mt_vector<mt_real> & xyz,
                               const mt_vector<mt_int> & vtxlist,
                               MT_USET<mt_int> & edge_vtx,
                               mt_real thr);

void identify_sharp_edges(const mt_meshdata & surfmesh,
                          const mt_vector<mt_real> & xyz,
                          const mt_mapping<mt_int> & ele2edge,
                          const MT_MAP<tuple<mt_int>, mt_int> & edges,
                          mt_real thr,
                          MT_USET<mt_int> & sharp_edges);


/**
* @brief Find nodes at edges that are only connected to one surface element.
*
* @param surfmesh The surface.
* @param border   The identified node indices.
*/
void identify_surface_border_nodes(const mt_meshdata & surfmesh,
                                   MT_USET<mt_int> & border);


/**
* @brief Add the interface between the regions defined by mgA and mgB into a splitlist
*
* @param mesh        The mesh we are operating on.
* @param mgA         The mesh graph A
* @param mgB         The mesh graph B
* @param Nidx        The current largest node index.
* @param splitlist   The splitlist we are populating.
*/
void add_interface_to_splitlist(const mt_meshdata & mesh,
                                const mt_meshgraph & mgA,
                                const mt_meshgraph & mgB,
                                mt_int & Nidx,
                                mt_vector<split_item> & splitlist);

/**
* @brief Add a surface into a splitlist
*
* @param mesh        The mesh we are operating on.
* @param surf        The surface defining the interface we split at
* @param vol_eidx    The volumetric elements we split at
* @param Nidx        The current largest node index.
* @param splitlist   The splitlist we are populating.
*/
void add_surface_to_splitlist(const mt_meshdata & mesh,
                              const mt_meshdata & surf,
                              const MT_USET<mt_int> & vol_eidx,
                              mt_int & Nidx,
                              mt_vector<split_item> & splitlist);

#endif

