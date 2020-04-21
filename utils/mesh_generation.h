#ifndef _MESH_GENERATION_H
#define _MESH_GENERATION_H

#include "mt_utils_base.h"
#include "mesh_utils.h"
#include "topology_utils.h"
#include "mesh_smoothing.h"
#include "tetgen/tetgen.h"

/**
* @brief Convert a meshtool mesh struct into a tetgenio datastruct
*/
void mesh_to_tetgenio(const mt_meshdata & mesh, tetgenio & io);

/**
* @brief Convert a tetgenio datastruct into a meshtool mesh struct
*/
void tetgenio_to_mesh(const tetgenio & io, mt_meshdata & mesh);

/**
* @brief Wrapper function that call tetgen to mesh a surface into a volume.
*
* @param surfmesh           The input surface mesh.
* @param volmesh            The output volumetric mesh.
* @param holes_xyz          Hole seeding points that get passed to tetgen.
* @param edge_length        Average input mesh edge length. Used for computing a constant elemet size.
* @param tetgen_qual        Desired mesh quality.
* @param size_field         Optional space dependent sizing function (e.g. average edge length in each node).
* @param preserve_boundary  Whether to preserve the boundary discretization.
*/
void mesh_with_tetgen(mt_meshdata & surfmesh,
                      mt_meshdata & volmesh,
                      mt_vector<float> & holes_xyz,
                      float edge_length,
                      float tetgen_qual,
                      mt_vector<mt_real> * size_field,
                      bool preserve_boundary);

/**
* @brief Generate Prism boundary layers from a triangle input surface.
*
* @param outer_surf         The input surface.
* @param num_layers         The number of boundary layers to generate.
* @param transition_scale   The scale to the normal direction tranisition distance.
* @param transition_inc     The (linear) increment in normal transition.
* @param inner_surf         The triangular final inner surface.
* @param output_mesh        The volumetric output mesh holding the boundary layers.
*/
void generate_boundary_layers(const mt_meshdata & outer_surf,
                              const int num_layers,
                              const mt_real transition_scale,
                              const mt_real transition_inc,
                              mt_meshdata & inner_surf,
                              mt_meshdata & output_mesh);

void generate_bbox_mesh(const mt_meshdata & inp_mesh, const mt_real scale, mt_meshdata & out_mesh);

#endif

