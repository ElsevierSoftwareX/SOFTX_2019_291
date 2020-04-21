/**
* @file mesh_smoothing.h
* @brief Mesh smoothing algorithms and classes.
* @author Aurel Neic
* @version
* @date 2016-12-13
*/


#ifndef _MESH_SMOOTHING
#define _MESH_SMOOTHING

#define SMOOTH_ITER_DEFAULT 100
#define SMOOTH_DEFAULT 0.15
#define EDGE_DETECT_DEFAULT 0.0
/// factor between forward and backward smoothing scaling
#define SMOOTH_FREQ_SCA 1.04f
/// default value for the mesh quality threshold
#define SMOOTH_THR_DEFAULT 0.95
#define QUAL_THR_CHK_SCA  0.5
/// every QUAL_THR_UPD, element quality thresholds will be refreshed
#define QUAL_THR_UPD     10
#define VOL_UPD          2


enum smoothing_t {LAPLACE_SMOOTHING = 0, CIRCUMSPHERE_SMOOTHING, CENTROID_SMOOTHING};

inline const char*
get_smoothing_type_string(const smoothing_t type)
{
  const char* ret = NULL;
  switch(type) {
    case LAPLACE_SMOOTHING:      ret = "laplace"; break;
    case CIRCUMSPHERE_SMOOTHING: ret = "circumspere"; break;
    case CENTROID_SMOOTHING:     ret = "centroid"; break;
  }

  return ret;
}

/**
* @brief Class for mesh smoothing while keeping mesh quality in check
*
*/
class quality_aware_smoother
{
  private:
  mt_vector<bool>     _applied;    ///< flag buffer to keep track of applied displacements
  mt_vector<mt_real>  _qual_thr;   ///< per-node quality threshold
  mt_vector<mt_real>  _vols_msh,   ///< main mesh volumes
                      _vols_mnfld; ///< manifold volumes

  /**
  * @brief initialize quality threshold
  *
  * @param checkmesh The mesh we check element quality on.
  * @param nodes     The nodes we are interested in because they will be smoothed.
  * @param inp_thr   The user-specified quality threshold
  */
  void init_qual_thr(const mt_meshdata & checkmesh,
                     const mt_vector<mt_int> & nodes,
                     const mt_real inp_thr);

  /**
  * @brief The forward smoothing iteration
  *
  * @param checkmesh   The mesh we check quality on.
  * @param smoothmesh  The main smoothing mesh.
  * @param manifold    The manifold mesh.
  * @param onManifold  Flag encoding which nodes are on the manifold.
  * @param nodes       Vector holding the nodes to smooth.
  * @param sca         Smoothing scale factor.
  * @param inp_thr     User-spcifed quality threshold.
  */
  void qual_smooth_iter_fwd(mt_meshdata & checkmesh,
                            mt_meshdata & smoothmesh,
                            mt_meshdata & manifold,
                            const mt_mask & onManifold,
                            const mt_vector<mt_int> & nodes,
                            const mt_real sca,
                            const mt_real inp_thr);

  void qual_smooth_iter_fwd(const mt_vector<mt_meshdata*> & manifolds,
                            const mt_vector<int> & mnfld_idx,
                            mt_vector<mt_real> & inp_xyz,
                            const mt_vector<mt_int> & nodes,
                            const mt_real sca,
                            const mt_real inp_thr);

  /**
  * @brief Backward smoothing iteration.
  *
  * @param checkmesh   The mesh we check quality on.
  * @param smoothmesh  The main smoothing mesh.
  * @param manifold    The manifold mesh.
  * @param onManifold  Flag encoding which nodes are on the manifold.
  * @param nodes       Vector holding the nodes to smooth.
  * @param sca         Smoothing scale factor.
  * @param inp_thr     User-spcifed quality threshold.
  */
  void qual_smooth_iter_bwd(mt_meshdata & checkmesh,
                            mt_meshdata & smoothmesh,
                            mt_meshdata & manifold,
                            const mt_mask & onManifold,
                            const mt_vector<mt_int> & nodes,
                            mt_real sca,
                            mt_real inp_thr);

  void qual_smooth_iter_bwd(const mt_vector<mt_meshdata*> & manifolds,
                            const mt_vector<int> & mnfld_idx,
                            mt_vector<mt_real> & inp_xyz,
                            const mt_vector<mt_int> & nodes,
                            const mt_real sca,
                            const mt_real inp_thr);
  public:

  smoothing_t smoothing_type = LAPLACE_SMOOTHING;

  /**
  * @brief Smooth a mesh with quality control.
  *
  * @param checkmesh   The mesh we check quality on.
  * @param smoothmesh  The main smoothing mesh.
  * @param manifold    The manifold mesh.
  * @param onManifold  Flag encoding which nodes are on the manifold.
  * @param nodes       Vector holding the nodes to smooth.
  * @param nsmooth     The number of smoothing operations.
  * @param sca         Smoothing scale factor.
  * @param inp_thr     User-spcifed quality threshold.
  */
  void operator()(mt_meshdata & checkmesh,
                  mt_meshdata & mesh,
                  mt_meshdata & manifold,
                  const mt_mask & onManifold,
                  const mt_vector<mt_int> & nodes,
                  const size_t nsmooth,
                  const mt_real sca,
                  const mt_real inp_thr);

  void operator()(const mt_vector<mt_meshdata*> & manifolds,
                  const mt_vector<int> & mnfld_idx,
                  mt_vector<mt_real> & inp_xyz,
                  const mt_vector<mt_int> & nodes,
                  const size_t  nsmooth,
                  const mt_real sca,
                  const mt_real inp_thr);
};




/**
* @brief Apply a given number of smoothing iterations.
*
* For volume preservation, smoothing is applied in forward and backward directions.
*
* @param [in]  mesh        The mesh.
* @param [in]  manifold    A mesh containing a manifold (either surface or line) of the mesh.
* @param [in]  onManifold  onManifold[i] tells wether node i is on the given manifold.
* @param [in]  nodes       The nodes to smooth.
* @param [in]  nsmooth     The number of smoothing iterations.
* @param [in]  sca         The smoothing coefficient.
*
* @post The vertex coordinates mesh.xyz have been updated.
*/
void smooth_nodes(mt_meshdata & mesh,
                  mt_meshdata & manifold,
                  const mt_mask & onManifold,
                  const mt_vector<mt_int> & nodes,
                  size_t nsmooth,
                  mt_real sca,
                  smoothing_t type = LAPLACE_SMOOTHING);



/**
* @brief Perform one data smoothing operation
*
* @param mesh      The mesh.
* @param nodes     Nodes to smooth data on.
* @param onSubset  Only smooth data using contributions from a nodal subset.
* @param sca       Smoothing scale factor.
* @param data      The data to smooth.
*/
template<class S>
void smooth_data_iter(const mt_meshdata & mesh,
                      const mt_vector<mt_int> & nodes,
                      const mt_vector<bool> & onSubset,
                      const mt_real sca,
                      mt_vector<S> & data)
{
  size_t nnodes = nodes.size();

  #ifdef OPENMP
  #pragma omp parallel for schedule(guided)
  #endif
  for(size_t i=0; i<nnodes; i++)
  {
    S avrg = S();
    mt_int nidx = nodes[i];

    mt_int start = mesh.n2n_dsp[nidx], stop = start + mesh.n2n_cnt[nidx];
    mt_int counter = 0;

    for(mt_int j=start; j<stop; j++)
    {
      mt_int n = mesh.n2n_con[j];
      if(onSubset[n]) {
        avrg += data[n];
        counter++;
      }
    }
    avrg /= counter;

    data[nidx] += (avrg - data[nidx]) * sca;
  }
}



/**
* @brief Smooth an arbitrary data vector.
*
* @tparam S   Data vector data-type.
*
* @param mesh           The mesh.
* @param nodes          Nodes to smooth data on.
* @param onSubset       Only smooth data using contributions from a nodal subset.
* @param data           The data to smooth.
* @param nsmooth        The number of smoothing operations.
* @param sca            Smoothing scale factor.
* @param data_is_nodal  Whether the input data is in nodal representation.
*/
template<class S>
void smooth_data(const mt_meshdata & mesh,
                 const mt_vector<mt_int> & nodes,
                 const mt_vector<bool> & onSubset,
                 const size_t nsmooth,
                 const mt_real sca,
                 mt_vector<S> & data,
                 const bool data_is_nodal)
{
  mt_vector<S> * d = &data;
  mt_vector<S> nodal_data;

  // if the input data is element-wise, we convert it to nodal
  if(!data_is_nodal) {
    elemData_to_nodeData(mesh, data, nodal_data);
    d = &nodal_data;
  }

  for(size_t sidx = 0; sidx < nsmooth; sidx++)
  {
    smooth_data_iter(mesh, nodes, onSubset, sca, *d);
    smooth_data_iter(mesh, nodes, onSubset, -(sca*SMOOTH_FREQ_SCA), *d);
  }

  // if the input data was element-wise, we convert the ouput nodal data back to
  // element data
  if(!data_is_nodal) {
    nodeData_to_elemData(mesh, nodal_data, data);
  }
}

/**
* @brief Function that bundles the usual volumetric smoothing steps
*
* @param mesh         The mesh we smooth
* @param stags        Groups of tags. The surface of each group will be smoothed.
* @param iter         Number of smoothing iterations.
* @param smth         Smoothing coefficient.
* @param edge_ang     Angle threshold defining sharp edges.
* @param max_qual     Maximum allowed elem qual.
* @param skip_lines   Whether to skip lines when smoothing.
* @param verbose      Whether to be verbose.
*/
void volumetric_smooth_from_tags(mt_meshdata & mesh,
                                 const mt_vector<MT_USET<mt_int> > & stags,
                                 const int iter,
                                 const mt_real smth,
                                 const mt_real edge_ang,
                                 const mt_real max_qual,
                                 const bool skip_lines,
                                 const bool verbose = false);

void surface_smooth(mt_meshdata & surfmesh,
                    const int iter,
                    const mt_real smth,
                    const mt_real edge_ang,
                    const mt_real max_qual,
                    const bool skip_lines,
                    const bool verbose);

/**
* @brief Smooth mesh vertices. Each vertex has also a direction that can be
*        removed from the smoothing update motion.
*
* @param mesh            The mesh to smooth.
* @param nod             The nodes to smooth in the mesh.
* @param ortho_dir       The nodal direction vectors used to orthogonalize the smoothing.
* @param ortho_scale     Scaling factor for the orthogonalization. (0 = no ortho, 1 = full ortho)
* @param only_positive   Only motion in direction of the provided vector field.
* @param smth            Smoothing coefficient.
* @param iter            Number of iterations.
*/
void directional_smoothing(mt_meshdata & mesh,
                           const mt_vector<mt_int> & nod,
                           const mt_vector<mt_real> & scale_dir,
                           mt_real scale,
                           mt_real ortho_scale,
                           const bool only_positive,
                           const mt_real smth,
                           const int iter);

void directional_smoothing_elem(mt_meshdata & mesh,
                           const mt_vector<mt_int> & nod,
                           const mt_vector<mt_real> & ortho_dir,
                           mt_real ortho_scale,
                           const bool only_positive,
                           const mt_real smth,
                           const int iter);
#endif

