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
#include "igb_utils.hpp"

static const std::string input_dynpts   = "-input_dynpts=";
static const std::string output_dynpts  = "-output_dynpts=";
static const std::string factor         = "-factor=";

struct change_dynpts_header_options
{
  std::string input_dynpts_file;
  std::string output_dynpts_file;
  std::string factor;
};

void print_change_dynpts_header_help()
 {
  fprintf(stderr, "change_dynpts_header: Dummy TEST thing\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to input dynpts file\n", input_dynpts.c_str());
  fprintf(stderr, "%s<float>\t(input) factor\n", factor.c_str());
  fprintf(stderr, "%s<path>\t (output) path to output dynpts file\n", output_dynpts.c_str());
  fprintf(stderr, "\n");
}

int change_dynpts_header_parse_options(int argc, char** argv, struct change_dynpts_header_options & opts)
{
  if(argc < 2) {
    print_change_dynpts_header_help();
    return 1;
  }

  // parse all enclose parameters -----------------------------------------------------------------
  for(int i=1; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, input_dynpts, opts.input_dynpts_file);
    if(!match) match = parse_param(param, output_dynpts, opts.output_dynpts_file);
    if(!match) match = parse_param(param, factor, opts.factor);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 3;
    }
  }

  // check if all relevant parameters have been set ---------------------------------------------------
  bool dynptsinputok = opts.input_dynpts_file.size() > 0,
       dynptstargetok = opts.output_dynpts_file.size() > 0,
       tflag = opts.factor.size() > 0;
  if( !(dynptsinputok && dynptstargetok && tflag) )
  {
    std::cerr << "Error: Insufficient parameters provided." << std::endl;
    print_change_dynpts_header_help();
    return 2;
  }

  return 0;
}

int main(int argc, char** argv)
{
  struct timeval t1, t2;
  struct change_dynpts_header_options opts;
  igb_header igb_head_from;
  igb_header igb_head_to;
  int ret = change_dynpts_header_parse_options(argc, argv, opts);
  if(ret != 0) return 1;

  //parse the igb header
  init_igb_header(opts.input_dynpts_file, igb_head_from);
  read_igb_header(igb_head_from);

  mt_real myfac = (mt_real)(atof(opts.factor.c_str()));
  //set the new header
  igb_head_to = igb_head_from;
  igb_head_to.v_dim_t = (float)(igb_head_from.v_dim_t * myfac);
  igb_head_to.v_inc_t = (float)(igb_head_from.v_inc_t * myfac);

  std::stringstream out;
  out << "Changing dynpts file header.v_dim_t" << opts.input_dynpts_file
      << " from " << igb_head_from.v_dim_t  << " to " << igb_head_to.v_dim_t
      << " and out to dynpts file " << opts.output_dynpts_file << " :";
  gettimeofday(&t1, NULL);
  mt_int tsteps = igb_head_from.v_t;
  mt_int bsize  = 1;

  std::vector<std::vector<float> > dynpts_data(bsize);
  PROGRESS<mt_int> progress(tsteps, out.str().c_str());

  for(mt_int i=0; i < tsteps; i++)
  {
    read_igb_block(dynpts_data, bsize, igb_head_from);
    write_igb_block(dynpts_data, igb_head_to);
    progress.next();
  }
  progress.finish();

  fclose(igb_head_from.fileptr);
  fclose(igb_head_to.fileptr);

  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
  return 0;
}
