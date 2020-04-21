/**
* @file mesh_refining.h
* @brief Mesh refinement classes.
* @author Aurel Neic
* @version
* @date 2017-06-02
*/
#ifndef _MESH_REFINING
#define _MESH_REFINING

#include "mesh_utils.h"
#include "topology_utils.h"

// The mesh quality threshold for resampling operations
#define RESAMPLE_QUAL_THR 0.95

/// inner iteration threshold for recomputing node and edge
/// indices
#define RECOMP_THR 64
// forward declarations for edge manager
class mt_edge_splitter;
class mt_edge_collapser;

/**
* @brief Class managing the edges of a mesh.
*
* It is used to conviniently refresh and access edge data.
*
*/
class mt_edge_manager
{
  friend class mt_edge_splitter;
  friend class mt_edge_collapser;

  protected:
  mt_meshdata           &       _mesh;          ///< the mesh we operate on
  mt_mapping<mt_int>            _ele2edge;      ///< elements-to-edges mapping
  MT_MAP<tuple<mt_int>, mt_int> _edge_map;      ///< edges of the provided tag regions

  public:
  /**
  * @brief Refresh the edge data-structs. Optionally remap mesh node indices.
  *
  * @param remap  Wheter to remap the mesh node indices.
  */
  void refresh_edge_data(bool remap)
  {
    if(remap)
      reindex_nodes(_mesh);

    _edge_map.clear();
    compute_edges(_mesh, _ele2edge, _edge_map);
    _ele2edge.transpose(); _ele2edge.setup_dsp();
  }

  /// getter function for the edges of the initially provided tag regions
  const MT_MAP<tuple<mt_int>, mt_int> & edges() const {
    return _edge_map;
  }

  /// getter function for elements-to-edges mapping
  const mt_mapping<mt_int> & elem2edges() const {
    return _ele2edge;
  }

  /// constructor
  mt_edge_manager(mt_meshdata & inp_mesh) : _mesh(inp_mesh)
  {
    bool remap_nod = false;
    refresh_edge_data(remap_nod);
  }
};

/**
* @brief Class for splitting edges in a mesh.
*/
class mt_edge_splitter
{
  private:
  mt_edge_manager   & _edge_manager;  ///< reference to the edge management class
  PROGRESS<size_t>  * _prg;           ///< a pointer to the current progress
  size_t              _edgeidx;       ///< edge index counter for adding new edges
  size_t              _inner_iter;    ///< sweep index counter

  /**
  * @brief Split a line into 2 lines
  *
  * @param eidx      Index of the element to split.
  * @param ele_start Start index for writing element data.
  * @param con_start Start index for writing connectivity data.
  * @param edge      Edge to split along.
  * @param nv        New vertex index
  */
  void split_line(const mt_int eidx,
                  const mt_int ele_start,
                  const mt_int con_start,
                  const tuple<mt_int> & edge,
                  const mt_int nv);
  /**
  * @brief Split a triangle into 2 triangles along a given edge.
  *
  * @param eidx      Index of the element to split.
  * @param ele_start Start index for writing element data.
  * @param con_start Start index for writing connectivity data.
  * @param edge      Edge to split along.
  * @param nv        New vertex index
  */
  void split_tri(const mt_int eidx,
                 const mt_int ele_start,
                 const mt_int con_start,
                 const tuple<mt_int> & edge,
                 const mt_int nv);

  /**
  * @brief Split a tetra into 2 tetras along a given edge.
  *
  * @param eidx      Index of the element to split.
  * @param ele_start Start index for writing element data.
  * @param con_start Start index for writing connectivity data.
  * @param edge      Edge to split along.
  * @param nv        New vertex index
  */
  void split_tet(const mt_int eidx,
                 const mt_int ele_start,
                 const mt_int con_start,
                 const tuple<mt_int> & edge,
                 const mt_int nv);

  /**
  * @brief Compute one splitting sweep.
  *
  * @param split_edge  The edges we split in the current sweep.
  * @param split_elem  The elems we split in the current sweep.
  */
  void split_iter(MT_MAP<mt_int, mixed_triple<mt_int,mt_int,float> > & split_edge,
                  MT_USET<mt_int> & split_elem);

  public:

  /**
  * @brief Constructor
  *
  * @param [in] inp_mesh The mesh to operate on.
  */
  mt_edge_splitter(mt_edge_manager & inp_manager) :
                   _edge_manager(inp_manager),
                   _inner_iter(0)
  {
    _edgeidx = _edge_manager.edges().size();
  }

  /**
  * @brief Split mesh edges iteratively.
  *
  * @param [in] edge_set  The edges to split.
  * @param [in] inp_prg   Progress output class.
  */
  void operator()(std::set<edgeele, edge_len_dsc> & edge_set, PROGRESS<size_t> & inp_prg);
};

/**
* @brief Class for collapsing edges
*/
class mt_edge_collapser
{
  private:
  mt_edge_manager  & _edge_manager;  ///< reference to the edge management class
  PROGRESS<size_t> * _prg;           ///< a pointer to the current progress

  /**
   * @brief Apply the collapsing on the mesh
   *
   * @param nodmap    A map holding the nodes to remove, the indices that they will be removed with and the new coords.
   * @param rem_elem  The elements to remove because share a collapsed edge.
   */
  void collapse_iter(MT_MAP<mt_int, mixed_tuple<mt_int, mt_point<mt_real> > > & nodmap,
      MT_USET<mt_int> & rem_elem);

  /**
   * @brief Check if tetrahedral elements will deteriorate
   * because of the planned collapsing.
   *
   * @param edgeelems The elements to check.
   * @param nodmask   A nodal mask holding whether a node is in nodmap.
   * @param nodmap    The planned collapsing information.
   *
   * @return Whether the collapsing will yield bad elements.
   */
   bool check_bad_tet(const MT_USET<mt_int> & edgeelems,
          const mt_vector<bool> & nodmask,
          const MT_MAP<mt_int, mixed_tuple<mt_int, mt_point<mt_real> > > & nodmap);

  /**
   * @brief Check if triangle elements will deteriorate
   * because of the planned collapsing.
   *
   * @param edgeelems The elements to check.
   * @param nodmap    The planned collapsing information.
   *
   * @return Whether the collapsing will yield bad elements.
   */
  bool check_bad_tri(const MT_USET<mt_int> & edgeelems,
      MT_MAP<mt_int, mixed_tuple<mt_int, mt_point<mt_real> > > & nodmap);

  public:
  /// constructor
  mt_edge_collapser(mt_edge_manager & inp_manager): _edge_manager(inp_manager)
  {}

  /**
  * @brief Try to collapse a given set of edges on a volumetric (tet) mesh.
  *
  * @param edge_set          The edges to collapse.
  * @param surf_nodes        The surface nodes (internal and external surfaces).
  * @param surf_nrml         The surface node normals.
  * @param fixed_surf_nodes  Surface nodes that cannot be moved.
  * @param surf_corr         The surface normal correlation coefficient. (0 = no correlation, 1 = full correlation).
  * @param qual_thr          Element quality threshold used to preserve mesh quality.
  * @param inp_prg           Progress display class.
  */
  void operator()(std::set<edgeele, edge_len_asc> & edge_set,
                  const mt_vector<bool> & surf_nodes,
                  const mt_vector<mt_real> & surf_nrml,
                  const mt_vector<bool> & fixed_surf_nodes,
                  const mt_real surf_corr,
                  PROGRESS<size_t> & inp_prg);

  /**
  * @brief Try to collapse a given set of edges on a surface (tri) mesh.
  *
  * @param edge_set          The edges to collapse.
  * @param surf_nrml         The surface node normals.
  * @param fixed_surf_nodes  Surface nodes that cannot be moved.
  * @param surf_corr         The surface normal correlation coefficient. (0 = no correlation, 1 = full correlation).
  * @param inp_prg           Progress display class.
  */
  void operator()(std::set<edgeele, edge_len_asc> & edge_set,
                  const mt_vector<mt_real> & surf_nrml,
                  const mt_vector<bool> & fixed_surf_nodes,
                  const mt_real surf_corr,
                  PROGRESS<size_t> & inp_prg);
};

#endif

