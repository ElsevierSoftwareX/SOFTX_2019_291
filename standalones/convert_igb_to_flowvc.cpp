/**
* @file convert_igb_to_flowvc.cpp
* @brief convert an igb file to a series of binary flowVC files
* @author Elias Karabelas
* @version
* @date 2018-11-30
*/


#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>

#include "mt_modes_base.h"
#include "mt_utils.h"

static const std::string input_igb   = "-input_igb=";
static const std::string output_flowvc_base  = "-output_basename=";
static const std::string begin_t  = "-begin_time=";
static const std::string end_t  = "-end_time=";


struct convert_flowvc_options
{
  std::string input_igb_file;
  std::string output_basename;
  std::string tbegin;
  std::string tend;
};

void print_convert_flowvc_help()
 {
  fprintf(stderr, "convert_igb_to_flowvc: Convert an igb file to a series of flowVC-compatible binary files\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to input igb file\n", input_igb.c_str());
  fprintf(stderr, "%s<int>\t (optional) time slice to start from. Default 0\n", begin_t.c_str());
  fprintf(stderr, "%s<int>\t (optional) time slice to end. Default -1 meaning all\n", end_t.c_str());
  fprintf(stderr, "%s<path>\t (output) output file basename\n", output_flowvc_base.c_str());
  fprintf(stderr, "\n");
}

int convert_flowvc_parse_options(int argc, char** argv, struct convert_flowvc_options & opts)
{
  if(argc < 2) {
    print_convert_flowvc_help();
    return 1;
  }

  // parse all enclose parameters -----------------------------------------------------------------
  for(int i=1; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, input_igb, opts.input_igb_file);
    if(!match) match = parse_param(param, output_flowvc_base, opts.output_basename);
    if(!match) match = parse_param(param, begin_t, opts.tbegin);
    if(!match) match = parse_param(param, end_t, opts.tend);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 3;
    }
  }


  // check if all relevant parameters have been set ---------------------------------------------------
  bool inputok = opts.input_igb_file.size() > 0,
       outputok = opts.output_basename.size() > 0;
  if( !(inputok && outputok) )
  {
    std::cerr << "Error: Insufficient parameters provided." << std::endl;
    print_convert_flowvc_help();
    return 2;
  }

  return 0;
}

int main(int argc, char** argv)
{
  struct timeval t1, t2;
  struct convert_flowvc_options opts;
  igb_header igb_head;
  int ret = convert_flowvc_parse_options(argc, argv, opts);
  if(ret != 0) return 1;

  int tbegin = opts.tbegin.size() > 0 ? atoi(opts.tbegin.c_str()) : 0;
  int tend = opts.tend.size() > 0 ? atoi(opts.tend.c_str()) : -1;



  std::vector<std::vector<float> > igb_data_input;
  //Read the input igb file----------------------------------------------------
  std::cout << "Reading " << opts.input_igb_file << std::endl;
  gettimeofday(&t1, NULL);
  read_igb_data(igb_data_input, igb_head, opts.input_igb_file);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  //-----------------------------------------------------------------------------
  const double tau = static_cast<double>(igb_head.v_inc_t);
  const int tsteps = igb_data_input.size();

  int real_tend = (tend == -1) ? tsteps : std::max(0,std::min(tsteps,tend));
  int real_tbegin = std::max(0,std::min(tsteps, tbegin));
  int tsteps_to_write = tsteps - tbegin;
  if(tend > 0)
   tsteps_to_write -= std::max(0,std::min(tsteps,tend));

  std::stringstream out;
  FILE * fileptr = NULL;
  size_t num_comp = (size_t)(Num_Components[igb_head.v_type]);
  size_t nvals = (size_t) (igb_head.v_x * igb_head.v_y * igb_head.v_z * num_comp);
  size_t byte_per_entry = sizeof(double);
  size_t checksum = (size_t)(tsteps_to_write * byte_per_entry) + (size_t) ( igb_head.v_x * igb_data_input.size() * num_comp);
  size_t count = 0;
  size_t bytes_written = 0;

  //build the vectors
  out << "Convert igb file " << opts.input_igb_file << " to flowVC series:";
  gettimeofday(&t1, NULL);
  PROGRESS<mt_int> progress(tsteps, out.str().c_str());
  for(mt_int i=real_tbegin; i < real_tend; i++)
  {
    progress.next();
    assert(igb_data_input[i].size() == nvals);
    std::string filename = opts.output_basename + "_vel." + std::to_string(i) + ".bin";
    fileptr = fopen(filename.c_str(), "wb");
    if(fileptr == NULL) {
     fprintf(stderr, "An IO error occured when opening file %s: %s\n\n", filename.c_str(), strerror(errno));
     exit(-1);
    }
    const double ts = (double)(igb_head.v_org_t) + tau * (double)(i);
    std::vector<double> wbuff;
    wbuff.assign(igb_data_input[i].begin(), igb_data_input[i].end());
    wbuff.insert(wbuff.begin(), ts);
    const char* wp = (const char*) wbuff.data();
    bytes_written = fwrite(wp, byte_per_entry, wbuff.size(), fileptr);
    count += bytes_written;
    fclose(fileptr);
  }
  progress.finish();
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  if(count != checksum) {
   fprintf(stderr, "Warning: Number of values written does not match blocksize\n");
   fprintf(stderr, "         written:  %zd\n", count);
   fprintf(stderr, "         expected: %zd\n", checksum);
   return -1;
  } else {
   return 0;
  }
}
