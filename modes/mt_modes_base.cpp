/**
* @file mt_modes_base.cpp
* @brief Basic meshtool modes functions
* @author Aurel Neic
* @version 
* @date 2017-08-16
*/

#include "mt_modes_base.h"

void print_usage(const char* exe)
{
  fprintf(stderr, "Error: Wrong usage.\n");
  fprintf(stderr, "Use: %s [help|convert|collect|clean|extract|generate|"
                  "itk|insert|interpolate|split|smooth|map|query|resample|restore|reindex] \n", exe);
  fprintf(stderr, "For help use: %s help\n", exe);
  fprintf(stderr, "\n");
}

void print_general_help(const char* exe)
{
  fprintf(stderr, "\nHELP:\n\n");
  fprintf(stderr, "With %s, the user uses \"modes\" to specify the operation he wants to perform on a mesh.\n", exe);
  fprintf(stderr, "Currently the supported modes are:\n");
  fprintf(stderr,
                  "clean quality\t\t\tdeform mesh elements to reach a certain quality threshold value\n"
                  "clean topology\t\t\tclean the mesh from bad topology definitions.\n"
                  "convert\t\t\t\tconvert between different mesh formats\n"
                  "collect\t\t\t\tmerge a mesh with datasets\n"
                  "extract data\t\t\tdata defined on a mesh is extracted for a given submesh\n"
                  "extract gradient\t\tcompute gradient and gradient magnitude of a scalar function on a mesh\n"
                  "extract mesh\t\t\ta submesh is extracted from a given mesh based on given element tags\n"
                  "extract myocard\t\t\tthe myocardium is extracted from a given mesh\n"
                  "extract surface\t\t\textract a sequence of surfaces defined by set operations on element tags\n"
                  "extract unreachable\t\textract elements that are unreachable through edge-traversal from a given start vertex\n"
                  "extract volume\t\t\textract elements inside a given box volume\n"
                  "generate fibres\t\t\tgenerate default fibers for a given mesh file\n"
                  "generate distancefield\t\tgenerate a distancefield between two surfaces\n"
                  "generate mesh\t\t\tgenerate a tetrahedral mesh from a list of nested triangle surfaces\n"
                  "insert data\t\t\tdata defined on a submesh is inserted back into a mesh\n"
                  "insert meshdata\t\t\tthe fiber and tag data of a mesh is inserted into another mesh\n"
                  "insert submesh\t\t\ta submesh is inserted back into a mesh and written to an output mesh\n"
                  "interpolate clouddata\t\tinterpolate data from a pointcloud onto a mesh\n"
                  "interpolate elemdata\t\tinterpolate element data from one mesh onto another\n"
                  "interpolate elem2node\t\tinterpolate data from elements onto nodes\n"
                  "interpolate node2elem\t\tinterpolate data from nodes onto elements\n"
                  "interpolate nodedata\t\tinterpolate nodal data from one mesh onto another\n"
                  "itk close\t\t\tApply closing (i.e. dilate-erode) algorithm to itk data\n"
                  "itk crop\t\t\tremove surrounding whitespace\n"
                  "itk dtype\t\t\tconvert datatype\n"
                  "itk extract\t\t\textract slices of an itk image stack\n"
                  "itk flip\t\t\tflip the voxel data along given axes\n"
                  "itk normalize\t\t\tNormalize voxel spacing\n"
                  "itk padding\t\t\tadd padding to voxel data\n"
                  "itk refine\t\t\trefine voxel data\n"
                  "itk sample\t\t\tcreate an itk image stack from sampeling surfaces\n"
                  "itk smooth\t\t\tSmooth the voxel data\n"
                  "map\t\t\t\tmap .vtx, .surf and .neubc files to the indexing of a submesh\n"
                  "merge surface\t\t\tmerge the geometry given by a closed surface mesh into a different mesh\n"
                  "merge meshes\t\t\tmerge two meshes, unifying co-located vertices\n"
                  "query bbox\t\t\tprint the bounding box of a given mesh\n"
                  "query curvature\t\t\tcompute the curvature of a surface\n"
                  "query edges\t\t\tprint several statistics related to the mesh edges\n"
                  "query graph\t\t\tprint the nodal connectivity graph\n"
                  "query idx\t\t\tprint indices in a proximity to a given coordinate\n"
                  "query insidepoint\t\tget a point inside a given closed surface\n"
                  "query quality\t\t\tprint mesh quality statistics\n"
                  "query tags\t\t\tprint the tags present in a given mesh\n"
                  "reindex\t\t\t\treindex a mesh to improve matrix bandwidth and cache efficiency\n"
                  "resample mesh\t\t\tresample a tetrahedral mesh to fit a given edge size range\n"
                  "resample purkinje\t\tresample purkinje cables as close as possible to a segment size\n"
                  "resample surfmesh\t\tresample a triangle mesh to fit a given edge size range\n"
                  "restore mapping\t\t\trestore nodal and element mapping for a submesh w.r.t. a reference mesh\n"
                  "smooth data\t\t\tsmooth data defined on a mesh\n"
                  "smooth mesh\t\t\tsmooth surfaces and volume of a mesh\n"
                  "smooth surface\t\t\tsmooth one or multiple surfaces of a mesh\n"
                  "split\t\t\t\tgenerate the split list for given split operations\n\n");
  fprintf(stderr, "For example, to extract a submesh from a given mesh, type\n\n");
  fprintf(stderr, "%s extract mesh <modeopts>\n\n", exe);
  fprintf(stderr, "Note, that <modeopts> denote the options of the \"extract mesh\" mode.\n");
  fprintf(stderr, "To view the options of any mode, just call that mode without any options.\n");
  fprintf(stderr, "\n");
}
