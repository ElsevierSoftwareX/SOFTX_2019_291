/**
* @file dynpts_fillup.cpp
* @brief Extend dynpts constantly from a given time instance on until the end.
* @author Elias Karabelas
* @version
* @date 2017-08-25
*/


#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>

#include "mt_modes_base.h"
#include "mt_utils.h"

static const std::string input_dynpts   = "-input_dynpts=";
static const std::string output_dynpts  = "-output_dynpts=";
static const std::string time_step           = "-time_step=";

struct fillup_dynpts_options
{
  std::string input_dynpts_file;
  std::string output_dynpts_file;
  std::string tstep;
};

void print_fillup_dynpts_help()
 {
  fprintf(stderr, "fillup_dynpts: Take a dynpts file and continue it from time step given constantly\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to input dynpts file\n", input_dynpts.c_str());
  fprintf(stderr, "%s<int>\t  (input) time step to fill up from\n", time_step.c_str());
  fprintf(stderr, "%s<path>\t (output) path to output dynpts file\n", output_dynpts.c_str());
  fprintf(stderr, "\n");
}

int fillup_dynpts_parse_options(int argc, char** argv, struct fillup_dynpts_options & opts)
{
  if(argc < 2) {
    print_fillup_dynpts_help();
    return 1;
  }

  // parse all enclose parameters -----------------------------------------------------------------
  for(int i=1; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, input_dynpts, opts.input_dynpts_file);
    if(!match) match = parse_param(param, output_dynpts, opts.output_dynpts_file);
    if(!match) match = parse_param(param, time_step, opts.tstep);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 3;
    }
  }

  // check if all relevant parameters have been set ---------------------------------------------------
  bool dynptsinputok = opts.input_dynpts_file.size() > 0,
       dynptstargetok = opts.output_dynpts_file.size() > 0,
       tflag = opts.tstep.size() > 0;
  if( !(dynptsinputok && dynptstargetok && tflag) )
  {
    std::cerr << "Error: Insufficient parameters provided." << std::endl;
    print_fillup_dynpts_help();
    return 2;
  }

  return 0;
}

int main(int argc, char** argv)
{
  struct timeval t1, t2;
  struct fillup_dynpts_options opts;
  igb_header igb_head_from; memset(&igb_head_from,0,sizeof(igb_header));
  igb_header igb_head_to;   memset(&igb_head_to,0,sizeof(igb_header));
  int ret = fillup_dynpts_parse_options(argc, argv, opts);
  if(ret != 0) return 1;

  //parse the igb header
  if(igb_head_from.igb_header_initialized == false)
  {
    //initialize the header
    init_igb_header(opts.input_dynpts_file, igb_head_from);
    //read in stuff from the igb header
    read_igb_header(igb_head_from);
  }
  //set the new header
  igb_head_to = igb_head_from;

  std::vector<std::vector<float> > dynpts_data_input;
  std::vector<std::vector<float> > dynpts_data_output;
  std::vector<float> dynpts_to_insert;
  //first insert the data into the reference mesh

  //Read the input dynpts file----------------------------------------------------
  //std::cout << "Reading " << opts.msh_input_dynpts_file << std::endl;
  //gettimeofday(&t1, NULL);
  //read_dynpts_(dynpts_data_input, igb_head_from, opts.msh_input_dynpts_file);
  //gettimeofday(&t2, NULL);
  //std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
  //-----------------------------------------------------------------------------

  mt_int tstep2start = (mt_int)(atoi(opts.tstep.c_str()));
  if(tstep2start > (igb_head_from.v_t -1) )
  {
    std::cerr << "Error: time_step " << tstep2start << " can not be bigger than "<< igb_head_from.v_t - 1 << std::endl;
    return -1;
  }
  std::stringstream out;
  //build the vectors
  out << "Continuing dynpts file " << opts.input_dynpts_file << " from time step " << tstep2start << "  and out to dynpts file " << opts.output_dynpts_file << " :";
  gettimeofday(&t1, NULL);
  mt_int tsteps = igb_head_from.v_t - 1; //dynpts_data_input.size();
  mt_int bsize  = 1;
  PROGRESS<mt_int> progress(tsteps, out.str().c_str());
  for(mt_int i=0; i < tsteps; i++)
  {
    progress.next();
    dynpts_data_output.resize(bsize);
    if(i < tstep2start)
    {
      read_igb_block(dynpts_data_input, bsize, igb_head_from);
      dynpts_data_output[0].assign(dynpts_data_input[0].begin(), dynpts_data_input[0].end());
      dynpts_to_insert.assign(dynpts_data_input[0].begin(), dynpts_data_input[0].end());
    }
    else
      dynpts_data_output[0].assign(dynpts_to_insert.begin(), dynpts_to_insert.end());
    write_igb_block(dynpts_data_output, igb_head_to);
  }
  progress.finish();
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
  return 0;
}
