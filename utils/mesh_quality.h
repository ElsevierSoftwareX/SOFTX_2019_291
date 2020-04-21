/**
* @file mesh_quality.h
* @brief Functions for computing mesh quality.
* @author Aurel Neic
* @version
* @date 2017-07-13
*/

#ifndef _MESH_QUALITY
#define _MESH_QUALITY

#include "data_structs.h"

#if defined QMETRIC_STD
#define QMETRIC tet_qmetric_volume
#define QMETRIC_NAME "tet_qmetric_volume"
#elif defined QMETRIC_COND
#define QMETRIC tet_qmetric_weightedShapeDist
#define QMETRIC_NAME "tet_qmetric_weightedShapeDist"
#elif defined QMETRIC_SINUS
#define QMETRIC tet_qmetric_minimalSinus
#define QMETRIC_NAME "tet_qmetric_minimalSinus"
#else
#define QMETRIC tet_qmetric_volume
#define QMETRIC_NAME "tet_qmetric_volume"
#endif

// forward declaration of the lookup table.
template<class V> class lookup_table;

/**
* @brief Compute tet element quality based on its volume.
*
* @param tet Tet connectivity.
* @param xyz Coorinates
*
* @return Element quality in (0,1) with 0 = worst, 1 = best
*/
mt_real tet_qmetric_volume(const mt_int* tet, const mt_real* xyz);

mt_real tri_qmetric_surface(const mt_int* tet, const mt_real* xyz);

/**
* @brief Compute tet element quality based on its shape distortion.
*
* @param tet Tet connectivity.
* @param xyz Coorinates
*
* @return Element quality in (0,1) with 0 = worst, 1 = best
*/
mt_real tet_qmetric_shapeDist(const mt_int* tet, const mt_real* xyz);

/**
* @brief Compute tet element quality based on its weighted shape distortion.
*
* @param tet Tet connectivity.
* @param xyz Coorinates
*
* @return Element quality in (0,1) with 0 = worst, 1 = best
*/
mt_real tet_qmetric_weightedShapeDist(const mt_int* tet, const mt_real* xyz);


/**
* @brief Compute tet element quality based on the biased sinus of the minimal diheadral angle.
         Biased means that quality is weighted towards small angles and less towards large angles 
*
* @param tet Tet connectivity.
* @param xyz Coorinates
*
* @return biased relative distance to optimal angle of 70.53°
*/

mt_real tet_qmetric_minimalSinus(const mt_int* tet, const mt_real* xyz);

mt_real distance_to_centroid(const mt_meshdata & mesh,
                            const mt_meshdata & manifold,
                            const mt_vector<bool> & onManifold,
                            const mt_int vtx);

mt_real minimal_surf_normal_correlation(const mt_meshdata & surfmesh,
                                       const mt_vector<mt_real> & xyz,
                                       const mt_int vtx);

/**
* @brief Get minimum angle sinus of all elements in a mesh
*
* @param mesh      The mesh.
* @param minangles Vector of minimum sinus angle per element.
*/
void mesh_min_angles(const mt_meshdata & mesh,
                       mt_vector<mt_real> & minangles);

/**
* @brief Compute volume for every mesh element.
*
* @param [in]  mesh The mesh.
* @param [in]  xyz  Vertex coordinates.
* @param [out] volumes Vector holding the quality of each element.
*/
void mesh_volumes(mt_meshdata & mesh,
                  mt_vector<mt_real> & xyz,
                  mt_vector<mt_real> & volumes);

inline mt_real element_quality(const elem_t & type, const mt_int* con, const mt_real* xyz)
{
  switch(type) {
    case Tetra:
      return QMETRIC(con, xyz);
    case Tri:
      return tri_qmetric_surface(con, xyz);
    default:
      fprintf(stderr, "%s error: unsupported element type! Aborting!\n", __func__);
      exit(1);
  }
}

template<class T, class S> inline
void get_thres_subset(const mt_vector<S> & vec,
                      const S threshold,
                      MT_USET<T> & set)
{
  for(size_t i=0; i < vec.size(); i++)
    if(vec[i] > threshold) set.insert(i);
}

/**
* @brief Compute mesh quality for every mesh element.
*
* @param [in]  mesh The mesh.
* @param [out] qual Vector holding the quality of each element.
*/
void mesh_quality(const mt_meshdata & mesh, mt_vector<mt_real> & qual);

void mesh_quality(const mt_meshdata & mesh, mt_real & min, mt_real & max, mt_real & avrg);

void mesh_quality(const mt_meshdata & mesh, mt_real & min, mt_real & max, mt_real & avrg, const mt_real & thr, mt_real & percentage);

vec3r quality_gradient(mt_meshdata & mesh, const mt_int vidx, const mt_real thr, const mt_real delta);

void improve_nodeset_gradient_method(mt_meshdata & mesh, const mt_meshdata & surf,
                                   const mt_mask & is_surf,
                                   const MT_USET<mt_int> & nodes, const mt_real thr,
                                   const mt_real relative_step);

/**
* @brief Compute the quality average of the elements connected to a node.
*
* @param mesh The mesh.
* @param nidx The node of interest.
*
* @return The average element quality.
*/
mt_real nbhd_quality_average(const mt_meshdata & mesh, const mt_int nidx);
/**
* @brief Compute the maximum quality of the elements connected to a node.
*
* @param mesh The mesh.
* @param nidx The node of interest.
*
* @return The average element quality.
*/
mt_real nbhd_quality_max(const mt_meshdata & mesh, const mt_int nidx);

/**
* @brief Compute the quality metric and element index of the worst quality element.
*
* @param [in]  mesh     The mesh.
* @param [in]  nidx     The node of interest.
* @param [out] max_qual The element quality.
* @param [out] max_eidx The element index.
*
* @return The average element quality.
*/
void nbhd_quality_max(const mt_meshdata & mesh, const mt_int nidx,
                      mt_real & max_qual, mt_int & max_eidx);

mt_real nbhd_quality_functional(const mt_meshdata & mesh,
                                const mt_int nidx,
                                const mt_real thr);

mt_real nbhd_smoothness_functional(const mt_meshdata & surfmesh,
                                   const mt_vector<mt_real> & xyz,
                                   const mt_int nidx,
                                   const mt_real thr);

void check_self_intersection(const mt_meshdata & mesh,
                             mt_vector<tuple<mt_int>> & intersec_edges,
                             mt_vector<mt_int> & intersec_eidx);

#endif
