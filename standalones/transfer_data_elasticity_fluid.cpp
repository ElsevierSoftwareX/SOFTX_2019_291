/**
* @file transfer_data_elasticity_fluid.cpp
* @brief Transfer dynpts file from elasticity submesh to reference mesh.
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

static const std::string submesh_input_par      = "-submsh_input=";
static const std::string submesh_input_dynpts   = "-submsh_input_dynpts=";
static const std::string submesh_target_par     = "-submsh_target=";
static const std::string submesh_target_dynpts  = "-submsh_target_dynpts=";

struct transfer_options{
  std::string msh_base;
  std::string submsh_target_base;
  std::string submsh_target_dynpts_file;
  std::string submsh_input_base;
  std::string submsh_input_dynpts_file;
};


void print_transfer_help()
{
  fprintf(stderr, "transfer_data: transfer dynpts file from elasticity submesh to reference mesh\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the reference mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to basename of input submesh\n", submesh_input_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to input submesh dynpts file\n", submesh_input_dynpts.c_str());
  fprintf(stderr, "%s<path>\t (input) path to basename of target submesh\n", submesh_target_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to target submesh dynpts file\n", submesh_target_dynpts.c_str());
  fprintf(stderr, "\n");
  fprintf(stderr, "Note that the files defining the submesh must include a *.nod file.\n"
                  "This file defines how the nodes of the submesh map back into the original mesh.\n"
                  "The *.nod file is generated when using the \"extract mesh\" mode.\n");
  fprintf(stderr, "Note that all meshes need to have consistent units!");
	fprintf(stderr, "\n");
}

int transfer_parse_options(int argc, char** argv, struct transfer_options & opts)
{
  if(argc < 2) {
    print_transfer_help();
    return 1;
  }

  // parse all enclose parameters -----------------------------------------------------------------
  for(int i=1; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh_par, opts.msh_base);
    if(!match) match = parse_param(param, submesh_input_par, opts.submsh_input_base);
    if(!match) match = parse_param(param, submesh_input_dynpts, opts.submsh_input_dynpts_file);
    if(!match) match = parse_param(param, submesh_target_par, opts.submsh_target_base);
    if(!match) match = parse_param(param, submesh_target_dynpts, opts.submsh_target_dynpts_file);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 3;
    }
  }
  fixBasename(opts.msh_base);
  fixBasename(opts.submsh_input_base);
  fixBasename(opts.submsh_target_base);


  // check if all relevant parameters have been set ---------------------------------------------------
  bool mshok = opts.msh_base.size() > 0,
       inputmshok = opts.submsh_input_base.size() > 0,
       targetmshok = opts.submsh_target_base.size() > 0,
       dynptsinputok = opts.submsh_input_dynpts_file.size() > 0,
       dynptstargetok = opts.submsh_target_dynpts_file.size() > 0;
  if( !(mshok && inputmshok && targetmshok && dynptsinputok && dynptstargetok) )
  {
    std::cerr << "Error: Insufficient parameters provided." << std::endl;
    print_transfer_help();
    return 2;
  }

  return 0;
}

int main(int argc, char** argv)
{
  struct timeval t1, t2;
  struct transfer_options opts;
  igb_header igb_head_from;
  igb_header igb_head_to;
  int ret = transfer_parse_options(argc, argv, opts);
  if(ret != 0) return 1;

  std::vector<std::vector<float> > dynpts_data_input;
  std::vector<std::vector<float> > dynpts_data_target;
  std::vector<std::vector<float> > dynpts_data_base;

  mt_vector<mt_real> target_mesh_xyz;
  mt_vector<mt_real> base_mesh_xyz;
  mt_vector<mt_real> input_mesh_xyz;

  mt_vector<float> target_mesh_xyz_float;
  mt_vector<float> base_mesh_xyz_float;
  mt_vector<float> input_mesh_xyz_float;

  //parse the igb header
  //initialize the header
  init_igb_header(opts.submsh_input_dynpts_file, igb_head_from);
  read_igb_header(igb_head_from);


  //first insert the data into the reference mesh
  size_t dpn = 3;

  //Read all the point data------------------------------------------------------
  std::cout << "Reading " << opts.msh_base + ".pts / .bpts" << std::endl;
  gettimeofday(&t1, NULL);
  readPoints_general(base_mesh_xyz, opts.msh_base);
  base_mesh_xyz_float.resize(base_mesh_xyz.size());
  for(size_t i=0; i < base_mesh_xyz.size(); i++)
    base_mesh_xyz_float[i] = (float)(base_mesh_xyz[i]);
  //size_t base_msh_numpts = base_mesh_xyz.size() / dpn;
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  std::cout << "Reading " << opts.submsh_input_base + ".pts / .bpts" << std::endl;
  gettimeofday(&t1, NULL);
  readPoints_general(input_mesh_xyz, opts.submsh_input_base);
  size_t smsh_input_numpts = input_mesh_xyz.size() / dpn;
  input_mesh_xyz_float.resize(input_mesh_xyz.size());
  for(size_t i=0; i < input_mesh_xyz.size(); i++)
    input_mesh_xyz_float[i] = (float)(input_mesh_xyz[i]);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  std::cout << "Reading " << opts.submsh_target_base + ".pts / .bpts" << std::endl;
  gettimeofday(&t1, NULL);
  readPoints_general(target_mesh_xyz, opts.submsh_target_base);
  size_t smsh_target_numpts = target_mesh_xyz.size() / dpn;
  target_mesh_xyz_float.resize(target_mesh_xyz.size());
  for(size_t i=0; i < target_mesh_xyz.size(); i++)
    target_mesh_xyz_float[i] = (float)(target_mesh_xyz[i]);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
  //------------------------------------------------------------------------------

  //Read the input dynpts file----------------------------------------------------
  std::cout << "Reading " << opts.submsh_input_dynpts_file << std::endl;
  gettimeofday(&t1, NULL);
  read_igb_data(dynpts_data_input, igb_head_from, opts.submsh_input_dynpts_file);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
  //-----------------------------------------------------------------------------

  //Read the nod files---------------------------------------------------------------
  std::cout << "Reading " << opts.submsh_input_base + NOD_EXT << std::endl;
  gettimeofday(&t1, NULL);
  mt_vector<mt_int> nod_input(smsh_input_numpts);
  binary_read(nod_input, opts.submsh_input_base + NOD_EXT);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  std::cout << "Reading " << opts.submsh_target_base + NOD_EXT << std::endl;
  gettimeofday(&t1, NULL);
  mt_vector<mt_int> nod_target(smsh_target_numpts);
  binary_read(nod_target, opts.submsh_target_base + NOD_EXT);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
  //--------------------------------------------------------------------------------

  //initialize the vectors

  size_t tsteps = dynpts_data_input.size();
  dynpts_data_target.resize(tsteps);
  dynpts_data_base.resize(tsteps);
  for(size_t i=0; i < tsteps; i++)
  {
    dynpts_data_target[i].assign(target_mesh_xyz_float.begin(), target_mesh_xyz_float.end());
    dynpts_data_base[i].assign(base_mesh_xyz_float.begin(), base_mesh_xyz_float.end());
  }

  std::stringstream out1;
  out1 << "Insert dynpts data from " << opts.submsh_input_dynpts_file << " to base mesh " << opts.msh_base << ": ";
  gettimeofday(&t1, NULL);
  PROGRESS<size_t> progress1(tsteps, out1.str().c_str());
  for(size_t i=0; i < tsteps; i++)
  {
    progress1.next();
    for(size_t j=0; j < smsh_input_numpts; j++)
    {
      const float* read = dynpts_data_input[i].data() + j*3;
      float * write	= dynpts_data_base[i].data() + nod_input[j]*3;
      write[0] = read[0];
      write[1] = read[1];
      write[2] = read[2];
    }
  }
  progress1.finish();
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  std::stringstream out2;
  out2 << "Extract dynpts data from " << opts.msh_base << " to target mesh " << opts.submsh_target_base << ": ";
  gettimeofday(&t1, NULL);
  PROGRESS<size_t> progress2(tsteps, out2.str().c_str());
  for(size_t i=0; i < tsteps; i++)
  {
    progress2.next();
    for(size_t j=0; j < smsh_target_numpts; j++)
    {
      const float* read = dynpts_data_base[i].data() + nod_target[j]*3;
      float * write	= dynpts_data_target[i].data() + j*3;
      write[0] = read[0];
      write[1] = read[1];
      write[2] = read[2];
    }
  }
  progress2.finish();
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  std::cout << "Write " << opts.submsh_target_dynpts_file << std::endl;
  //adapt igb header file
  gettimeofday(&t1, NULL);
  igb_head_to          = igb_head_from; // does not copy the fileptr
  igb_head_to.filename = opts.submsh_target_dynpts_file;
  igb_head_to.fileptr  = NULL;
  igb_head_to.v_x      = smsh_target_numpts;
  write_igb_header(igb_head_to);
  write_igb_data(dynpts_data_target, igb_head_to);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
  return 0;
}
