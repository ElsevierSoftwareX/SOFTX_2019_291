/**
* @file mesh_mode_funcs.h
* @brief Main meshtool include file.
* @author Aurel Neic
* @version
* @date 2017-02-13
*/

#include "mt_utils.h"

#ifdef OPENMP
#include <omp.h>
#endif

// global command line option variables
static const std::string mesh_par         = "-msh=";
static const std::string submesh_par      = "-submsh=";
static const std::string outmesh_par      = "-outmsh=";
static const std::string mesh_data_par    = "-msh_data=";
static const std::string submesh_data_par = "-submsh_data=";
static const std::string tags_par         = "-tags=";
static const std::string coord_par        = "-coord=";
static const std::string oper_par         = "-op=";
static const std::string thr_par          = "-thr=";
static const std::string surf_par         = "-surf=";
static const std::string file_par         = "-files=";
static const std::string outdir_par       = "-outdir=";
static const std::string inp_mesh_par     = "-imsh=";
static const std::string out_mesh_par     = "-omsh=";
static const std::string inp_format_par   = "-ifmt=";
static const std::string out_format_par   = "-ofmt=";
static const std::string smooth_par       = "-smth=";
static const std::string iter_par         = "-iter=";
static const std::string split_par        = "-split=";
static const std::string idx_par          = "-idx=";
static const std::string edge_par         = "-edge=";
static const std::string size_par         = "-size=";
static const std::string lower_size_par   = "-lower_size=";
static const std::string corr_par         = "-corr_thr=";
static const std::string idat_par         = "-idat=";
static const std::string odat_par         = "-odat=";
static const std::string min_par          = "-min=";
static const std::string max_par          = "-max=";
static const std::string avrg_par         = "-avrg=";
static const std::string lvl_par          = "-lvl=";
static const std::string fix_par          = "-fix_bnd=";
static const std::string ins_tag_par      = "-ins_tag=";
static const std::string scale_par        = "-scale=";
static const std::string mode_par         = "-mode=";
static const std::string mesh1_par        = "-msh1=";
static const std::string mesh2_par        = "-msh2=";
static const std::string vtx_par          = "-vtx=";
static const std::string info_par         = "-info=";
static const std::string uvc_par          = "-uvc=";
static const std::string out_par          = "-out=";

static const char tag_separator           = ',';



/**
* @brief Enumeration encoding the mode type of extraction and insertion modes.
*/
enum op_type{OP_MESH, OP_DATA, OP_SURF, OP_MYO, OP_VOL, OP_UR, OP_GRAD, OP_OVERLP, OP_TAGS, OP_ISOSURF};

/**
* @brief Basic usage message
*
* @param [in]  argc Arguments count.
* @param [in]  argv Arguments string-array.
*/
void print_usage(const char* exe);

/**
* @brief Help message code
*
* @param argc
* @param argv
*/
void print_general_help(const char* exe);

