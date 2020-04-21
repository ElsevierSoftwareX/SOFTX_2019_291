/**
* @file convert_xdynpts_to_displdynpts.cpp
* @brief Convert vertex positions into a displacement for dynpts files.
* @author Elias Karabelas
* @version
* @date 2017-08-25
*/


#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>

#include "mt_modes_base.h"
#include "igb_utils.hpp"

static const std::string mesh_input_dynpts   = "-msh_input_dynpts=";
static const std::string mesh_output_dynpts  = "-msh_output_dynpts=";
static const std::string convert_x_to_u      = "-convert_to_displacement=";

struct convert_dynpts_options
{
  std::string msh_base;
  std::string msh_input_dynpts_file;
  std::string msh_output_dynpts_file;
  std::string convert_to_displacement;
};

void print_convert_dynpts_help()
 {
  fprintf(stderr, "convert_dynpts_data: convert dynpts files give on a mesh to a displacement dynpts file and vice versa\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to input dynpts file\n", mesh_input_dynpts.c_str());
  fprintf(stderr, "%s<bool>\t (input) convert x.dynpt to displ.dynpt\n", convert_x_to_u.c_str());
  fprintf(stderr, "%s<path>\t (output) path to output dynpts file\n", mesh_output_dynpts.c_str());
  fprintf(stderr, "\n");
  fprintf(stderr, "Note that all files need to be specified in the same physical unit.\n"
                  "If convert x to u is set it is assumed that the input file contains the new points in every time steps.\n"
                  "if not it is assumed we transfer a displacement file to original dynpts file.\n");
  fprintf(stderr, "\n");
}

int convert_dynpts_parse_options(int argc, char** argv, struct convert_dynpts_options & opts)
{
  if(argc < 2) {
    print_convert_dynpts_help();
    return 1;
  }

  // parse all enclose parameters -----------------------------------------------------------------
  for(int i=1; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh_par, opts.msh_base);
    if(!match) match = parse_param(param, mesh_input_dynpts, opts.msh_input_dynpts_file);
    if(!match) match = parse_param(param, mesh_output_dynpts, opts.msh_output_dynpts_file);
    if(!match) match = parse_param(param, convert_x_to_u, opts.convert_to_displacement);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 3;
    }
  }
  fixBasename(opts.msh_base);


  // check if all relevant parameters have been set ---------------------------------------------------
  bool mshok = opts.msh_base.size() > 0,
       dynptsinputok = opts.msh_input_dynpts_file.size() > 0,
       dynptstargetok = opts.msh_output_dynpts_file.size() > 0,
       boolflag = opts.convert_to_displacement.size() > 0;
  if( !(mshok && dynptsinputok && dynptstargetok && boolflag) )
  {
    std::cerr << "Error: Insufficient parameters provided." << std::endl;
    print_convert_dynpts_help();
    return 2;
  }

  return 0;
}

int main(int argc, char** argv)
{
  struct timeval t1, t2;
  struct convert_dynpts_options opts;
  igb_header igb_head_from;
  igb_header igb_head_to;
  int ret = convert_dynpts_parse_options(argc, argv, opts);
  if(ret != 0) return 1;

  //parse the igb header
  //initialize the header
  init_igb_header(opts.msh_input_dynpts_file, igb_head_from);
  read_igb_header(igb_head_from);

  //set the new header
  igb_head_to          = igb_head_from; // does not copy the fileptr
  igb_head_to.filename = opts.msh_output_dynpts_file;
  igb_head_to.fileptr  = NULL;
  write_igb_header(igb_head_to);

  mt_vector<mt_real> base_mesh_xyz;

  //first insert the data into the reference mesh
  mt_int dpn = 3;

  //Read the point data of the mesh ------------------------------------------------------
  std::cout << "Reading " << opts.msh_base + ".pts / .bpts" << std::endl;
  gettimeofday(&t1, NULL);
  readPoints_general(base_mesh_xyz, opts.msh_base);
  //size_t base_msh_numpts = base_mesh_xyz.size() / dpn;
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  //Read the input dynpts file----------------------------------------------------
  //std::cout << "Reading " << opts.msh_input_dynpts_file << std::endl;
  //gettimeofday(&t1, NULL);
  //read_dynpts_(dynpts_data_input, igb_head_from, opts.msh_input_dynpts_file);
  //gettimeofday(&t2, NULL);
  //std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
  //-----------------------------------------------------------------------------

  int convert_x_to_u_flag = atoi(opts.convert_to_displacement.c_str());
  std::stringstream out;
  //build the vectors
  if(convert_x_to_u_flag != 0)
    out << "Convert position dynpts file " << opts.msh_input_dynpts_file << " to displacement dynpts file " << opts.msh_output_dynpts_file << " :";
  else
    out << "Convert displacement dynpts file " << opts.msh_input_dynpts_file << " to position dynpts file " << opts.msh_output_dynpts_file << " :";
  gettimeofday(&t1, NULL);
  mt_int tsteps = igb_head_from.v_t; //dynpts_data_input.size();
  mt_int bsize  = 1;
  std::vector<std::vector<float> > dynpts_data_input(bsize);
  std::vector<std::vector<float> > dynpts_data_output(bsize);
  dynpts_data_output[0].resize(base_mesh_xyz.size());

  assert(int(dynpts_data_output[0].size()) == igb_head_to.v_x * dpn);
  const mt_int nvtx = base_mesh_xyz.size() / dpn;

  PROGRESS<mt_int> progress(tsteps, out.str().c_str());
  for(mt_int i=0; i < tsteps; i++)
  {
    progress.next();
    read_igb_block(dynpts_data_input, bsize, igb_head_from);

    for(mt_int j=0; j < nvtx; j++)
    {
        const float *   read_dynpts  = dynpts_data_input[0].data() + j*dpn;
        const mt_real * read_xyz     = base_mesh_xyz.data() + j*dpn;
        float *         write_dynpts = dynpts_data_output[0].data() + j*dpn;
        if(convert_x_to_u_flag !=0)
        {
          //calc displacement as x^new - x^base
          write_dynpts[0] = read_dynpts[0] - (float)(read_xyz[0]);
          write_dynpts[1] = read_dynpts[1] - (float)(read_xyz[1]);
          write_dynpts[2] = read_dynpts[2] - (float)(read_xyz[2]);
        }
        else
        {
          //calc new position x^new  = x^base + displacement
          write_dynpts[0] = read_dynpts[0] + (float)(read_xyz[0]);
          write_dynpts[1] = read_dynpts[1] + (float)(read_xyz[1]);
          write_dynpts[2] = read_dynpts[2] + (float)(read_xyz[2]);
        }
    }
    write_igb_block(dynpts_data_output, igb_head_to);
  }
  progress.finish();

  fclose(igb_head_from.fileptr);
  fclose(igb_head_to.fileptr);

  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
  return 0;
}
