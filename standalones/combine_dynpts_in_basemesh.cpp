/**
* @file combine_dynpts_in_basemesh.cpp
* @brief Standalone util to combine dynpts into a mesh.
* @author Elias Karabelas
* @version
* @date 2017-08-25
*/

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>

#include "mt_modes_base.h"
#include "igb_utils.hpp"

static const std::string submesh_input_first_par      = "-submsh_first=";
static const std::string submesh_input_first_dynpts   = "-submsh_first_dynpts=";
static const std::string submesh_input_second_par     = "-submsh_second=";
static const std::string submesh_input_second_dynpts  = "-submsh_second_dynpts=";
static const std::string remove_intersection	      = "-remove_second_submsh_intersection=";

static const std::string basemesh_dynpts     = "-msh_dynpts=";

struct combine_dynpts_options{
  std::string msh_base;
  std::string msh_dynpts_file;
  std::string submsh_input_first_base;
  std::string submsh_input_first_dynpts_file;
  std::string submsh_input_second_base;
  std::string submsh_input_second_dynpts_file;
  std::string remove_intersect;
};

void print_combine_dynpts_help()
{
  fprintf(stderr, "combine_dynpts: combine two dynpts file on submeshes in one dynpts file for base mesh\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the reference mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to basename of first input submesh\n", submesh_input_first_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to first input submesh dynpts file\n", submesh_input_first_dynpts.c_str());
  fprintf(stderr, "%s<path>\t (input) path to basename of second input submesh\n", submesh_input_second_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to second input submesh dynpts file\n", submesh_input_second_dynpts.c_str());
  fprintf(stderr, "%s<bool>\t (input) remove the intersection of the nod file from first mesh and second mesh from second one\n", remove_intersection.c_str());
  fprintf(stderr, "%s<path>\t (output) path to target basename dynpts file\n", basemesh_dynpts.c_str());
  fprintf(stderr, "\n");
  fprintf(stderr, "Note that the files defining the submesh must include a *.nod file.\n"
                  "This file defines how the nodes of the submesh map back into the original mesh.\n"
                  "The *.nod file is generated when using the \"extract mesh\" mode.\n");
  fprintf(stderr, "Note that all meshes need to have consistent units!\n"
                  "Note that the union of the two submeshes must yield the total mesh!\n");
	fprintf(stderr, "\n");
}

int combine_dynpts_parse_options(int argc, char** argv, struct combine_dynpts_options & opts)
{
  if(argc < 2) {
    print_combine_dynpts_help();
    return 1;
  }

  // parse all enclose parameters -----------------------------------------------------------------
  for(int i=1; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh_par, opts.msh_base);
    if(!match) match = parse_param(param, submesh_input_first_par, opts.submsh_input_first_base);
    if(!match) match = parse_param(param, submesh_input_first_dynpts, opts.submsh_input_first_dynpts_file);
    if(!match) match = parse_param(param, submesh_input_second_par, opts.submsh_input_second_base);
    if(!match) match = parse_param(param, submesh_input_second_dynpts, opts.submsh_input_second_dynpts_file);
    if(!match) match = parse_param(param, basemesh_dynpts, opts.msh_dynpts_file);
    if(!match) match = parse_param(param, remove_intersection, opts.remove_intersect);


    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 3;
    }
  }
  fixBasename(opts.msh_base);
  fixBasename(opts.submsh_input_first_base);
  fixBasename(opts.submsh_input_second_base);


  // check if all relevant parameters have been set ---------------------------------------------------
  bool mshok = opts.msh_base.size() > 0,
       inputmsh1ok = opts.submsh_input_first_base.size() > 0,
       inputmsh2ok = opts.submsh_input_second_base.size() > 0,
       dynptsinput1ok = opts.submsh_input_first_dynpts_file.size() > 0,
       dynptsinput2ok = opts.submsh_input_second_dynpts_file.size() > 0,
       dynptstargetok = opts.msh_dynpts_file.size() > 0,
       boolflag = opts.remove_intersect.size() > 0;
  if( !(mshok && inputmsh1ok && inputmsh2ok && dynptsinput1ok && dynptsinput2ok && dynptstargetok && boolflag) )
  {
    std::cerr << "Error: Insufficient parameters provided." << std::endl;
    print_combine_dynpts_help();
    return 2;
  }

  return 0;
}

int main(int argc, char** argv)
{
  struct timeval t1, t2;
  struct combine_dynpts_options opts;
  igb_header igb_head_submsh_one;
  igb_header igb_head_submsh_two;
  igb_header igb_head_basemsh;

  int ret = combine_dynpts_parse_options(argc, argv, opts);
  if(ret != 0) return 1;

  std::vector<std::vector<float> > dynpts_data_submsh_one;
  std::vector<std::vector<float> > dynpts_data_submsh_two;
  std::vector<std::vector<float> > dynpts_data_base;

  mt_vector<mt_real> base_mesh_xyz;
  mt_vector<mt_real> submesh_one_xyz;
  mt_vector<mt_real> submesh_two_xyz;
  mt_vector<float> base_mesh_xyz_float;

  //first insert the data into the reference mesh
  mt_int dpn = 3;

  //Read all the point data------------------------------------------------------
  std::cout << "Reading " << opts.msh_base + ".pts / .bpts" << std::endl;
  gettimeofday(&t1, NULL);
  readPoints_general(base_mesh_xyz, opts.msh_base);
  mt_int base_msh_numpts = base_mesh_xyz.size() / dpn;
  base_mesh_xyz_float.assign(base_mesh_xyz.begin(), base_mesh_xyz.end());
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  std::cout << "Reading " << opts.submsh_input_first_base + ".pts / .bpts" << std::endl;
  gettimeofday(&t1, NULL);
  readPoints_general(submesh_one_xyz, opts.submsh_input_first_base);
  mt_int smsh_one_input_numpts = submesh_one_xyz.size() / dpn;
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  std::cout << "Reading " << opts.submsh_input_second_base + ".pts / .bpts" << std::endl;
  gettimeofday(&t1, NULL);
  readPoints_general(submesh_two_xyz, opts.submsh_input_second_base);
  mt_int smsh_two_input_numpts = submesh_two_xyz.size() / dpn;
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  //------------------------------------------------------------------------------

  //Read the input dynpts file first submesh -------------------------------------
  std::cout << "Reading " << opts.submsh_input_first_dynpts_file << std::endl;
  gettimeofday(&t1, NULL);
  read_igb_data(dynpts_data_submsh_one, igb_head_submsh_one, opts.submsh_input_first_dynpts_file);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
  //-----------------------------------------------------------------------------

  //Read the input dynpts file second submesh -------------------------------------
  std::cout << "Reading " << opts.submsh_input_second_dynpts_file << std::endl;
  gettimeofday(&t1, NULL);
  read_igb_data(dynpts_data_submsh_two, igb_head_submsh_two, opts.submsh_input_second_dynpts_file);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
  //-----------------------------------------------------------------------------
  mt_int tsteps     = dynpts_data_submsh_one.size();
  mt_int tsteps_chk = dynpts_data_submsh_two.size();

  if( tsteps != tsteps_chk ) {
    std::cerr << "Error: Number of tsteps for the two dynpts files is not equal! ("
              << tsteps << " != " << tsteps_chk << ")" << std::endl;
    return -1;
  }

  //Read the nod files---------------------------------------------------------------
  std::cout << "Reading " << opts.submsh_input_first_base + NOD_EXT << std::endl;
  gettimeofday(&t1, NULL);
  mt_vector<mt_int> nod_input_one(smsh_one_input_numpts);
  binary_read(nod_input_one, opts.submsh_input_first_base + NOD_EXT);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  std::cout << "Reading " << opts.submsh_input_second_base + NOD_EXT << std::endl;
  gettimeofday(&t1, NULL);
  mt_vector<mt_int> nod_input_two(smsh_two_input_numpts);
  binary_read(nod_input_two, opts.submsh_input_second_base + NOD_EXT);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
  //--------------------------------------------------------------------------------
  int exclude_flag = std::atoi(opts.remove_intersect.c_str());
  std::vector<mt_int> nod_vec_one(smsh_one_input_numpts);
  for(mt_int i=0; i < smsh_one_input_numpts; i++)
    nod_vec_one[i] = nod_input_one[i];
  std::vector<mt_int> nod_vec_two(smsh_two_input_numpts);
  for(mt_int i=0; i < smsh_two_input_numpts; i++)
    nod_vec_two[i] = nod_input_two[i];
  std::map<mt_int,mt_int> backwards_map_nod_vec_two;
  for(mt_int i=0; i < smsh_two_input_numpts; i++)
    backwards_map_nod_vec_two.emplace(std::make_pair(nod_vec_two[i],i));

  std::vector<mt_int> intersection;
  std::vector<mt_int> nod_vec_two_minus_intersection;
  std::set_intersection(nod_vec_one.begin(), nod_vec_one.end(), nod_vec_two.begin(), nod_vec_two.end(), std::back_inserter(intersection));
  std::set_difference(nod_vec_two.begin(), nod_vec_two.end(), intersection.begin(), intersection.end(), std::inserter(nod_vec_two_minus_intersection, nod_vec_two_minus_intersection.begin()));

  //initialize the dynpts vector for the base mesh
  dynpts_data_base.resize(tsteps);

  //std::cout << "Combine dynpts " << opts.submsh_input_first_dynpts_file << " with " << opts.submsh_input_second_dynpts_file  << " into " << opts.msh_dynpts_file << std::endl;
  gettimeofday(&t1, NULL);
  std::stringstream out;
  out << "Combine dynpts " << opts.submsh_input_first_dynpts_file << " with "
      << opts.submsh_input_second_dynpts_file  << " into " << opts.msh_dynpts_file << " :";
  PROGRESS<mt_int> progress(tsteps, out.str().c_str());

  for(mt_int i=0; i < tsteps; i++) {
    progress.next();
    //Always use the coordinates as initial value

    dynpts_data_base[i].assign(base_mesh_xyz_float.begin(), base_mesh_xyz_float.end());
    //first submesh
    for(mt_int j=0; j < smsh_one_input_numpts; j++)
    {
      const float* read  = dynpts_data_submsh_one[i].data() + j*dpn;
      float *      write = dynpts_data_base[i].data() + nod_input_one[j]*dpn;
      write[0] = read[0];
      write[1] = read[1];
      write[2] = read[2];
    }
    //second submesh
    if(exclude_flag) {
      for(size_t j=0; j < nod_vec_two_minus_intersection.size(); j++)
      {
        auto search = backwards_map_nod_vec_two.find(nod_vec_two_minus_intersection[j]);
        if(search != backwards_map_nod_vec_two.end())
        {
          const float* read  = dynpts_data_submsh_two[i].data() + (search->second)*dpn;
          float *      write = dynpts_data_base[i].data() + nod_vec_two_minus_intersection[j]*dpn;
          write[0] = read[0];
          write[1] = read[1];
          write[2] = read[2];
        }
      }
    }
    else {
      for(mt_int j=0; j < smsh_two_input_numpts; j++)
      {
        const float* read  = dynpts_data_submsh_two[i].data() + j*dpn;
        float *      write = dynpts_data_base[i].data() + nod_input_two[j]*dpn;
        write[0] = read[0];
        write[1] = read[1];
        write[2] = read[2];
      }
    }
  }
  progress.finish();
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  std::cout << "Write " << opts.msh_dynpts_file << std::endl;
  //adapt igb header file
  gettimeofday(&t1, NULL);
  igb_head_basemsh          = igb_head_submsh_one; // does not copy the fileptr
  igb_head_basemsh.filename = opts.msh_dynpts_file;
  igb_head_basemsh.fileptr  = NULL;
  igb_head_basemsh.v_x      = base_msh_numpts;
  write_igb_header(igb_head_basemsh);
  write_igb_data(dynpts_data_base, igb_head_basemsh);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
  return 0;
}
