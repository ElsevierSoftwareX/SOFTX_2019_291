#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <string>

#include "mt_modes_base.h"
#include "lookup_table.hpp"


struct spctm_options {

  std::string temp_input;
  std::string spacial_input;
  std::string oigb;
  std::string dt;
  std::string amp;
};

static const std::string spc_par  = "-spc=";
static const std::string tm_par   = "-tm=";
static const std::string oigb_par = "-oigb=";
static const std::string dt_par   = "-dt=";
static const std::string amp_par  = "-amp=";


void print_spctm_help()
{
  fprintf(stderr, "spctm_data: Compute a space-time function by combining spactial data with a scaling trace.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to the spacial data file.\n", spc_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to the temporal scaling trace file.\n", tm_par.c_str());
  fprintf(stderr, "%s<float>\t (input) time step.\n", dt_par.c_str());
  fprintf(stderr, "%s<float>\t (input) amplitude of scaling.\n", amp_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to the output igb file.\n", oigb_par.c_str());
  fprintf(stderr, "\n");
}

int spctm_parse_options(int argc, char** argv, spctm_options & opts)
{
  if(argc < 2) {
    print_spctm_help();
    return 1;
  }

  // parse all parameters -----------------------------------------------------------------
  for(int i=1; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, spc_par, opts.spacial_input);
    if(!match) match = parse_param(param, tm_par, opts.temp_input);
    if(!match) match = parse_param(param, oigb_par, opts.oigb);
    if(!match) match = parse_param(param, dt_par, opts.dt);
    if(!match) match = parse_param(param, amp_par, opts.amp);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 3;
    }
  }

  // check if all relevant parameters have been set ---------------------------------------------------

  if( !(opts.spacial_input.size() && opts.temp_input.size() && opts.oigb.size() &&
        opts.dt.size() && opts.amp.size()) )
  {
    std::cerr << "Error: Insufficient parameters provided." << std::endl;
    print_spctm_help();
    return 2;
  }

  return 0;
}

int main(int argc, char** argv)
{
  struct timeval t1, t2;
  struct spctm_options opts;
  // struct mt_meshdata mesh;

  int ret = spctm_parse_options(argc, argv, opts);
  if (ret != 0) return 1;

  // data_idx may be:
  // 0 : scalar dat file
  // 1 : scalar igb file
  // 2 : vec file
  // 3 : vec3 igb file
  short data_idx = -1;
  mt_vector<mt_real>            idat,     odat;
  mt_vector<mt_point<mt_real> > idat_vec, odat_vec;
  igb_header igb, igb_out;
  setup_data_format(opts.spacial_input, opts.oigb, data_idx, igb, igb_out);

  mt_vector<mt_real> spacial_data;
  std::cout << "Reading spatial data: " << opts.spacial_input << std::endl;
  gettimeofday(&t1, NULL);
  switch(data_idx) {
    case 0:
      read_vector_ascii(spacial_data, opts.spacial_input);
      break;

    case 1: {
      std::vector<std::vector<mt_real> > rbuff;
      read_igb_block(rbuff, 1, igb);
      spacial_data.assign(rbuff[0].begin(), rbuff[0].end());
      break;
    }

    default:
      fprintf(stderr, "Unsupported spatial input data format. Aborting!\n");
      exit(1);
  }

  mt_real maxval = fabs(spacial_data[0]);
  for(size_t i=0; i<spacial_data.size(); i++) {
    mt_real val = fabs(spacial_data[i]);
    if(maxval < val) maxval = val;
  }
  for(size_t i=0; i<spacial_data.size(); i++) spacial_data[i] /= maxval;

  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  std::cout << "Reading trace file: " << opts.temp_input << std::endl;
  gettimeofday(&t1, NULL);
  dyn_lookup_table<mt_real> lookup(opts.temp_input);
  lookup.normalize();

  mt_real start_time = lookup.start_x(), end_time = lookup.end_x();
  mt_real dt  = atof(opts.dt.c_str());
  mt_real amp = atof(opts.amp.c_str());
  printf("Input trace: start = %.2f, end = %.2f, dt = %.2f, amplitude = %.2f\n",
         start_time, end_time, dt, amp);

  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  int num_int = (end_time - start_time) / dt;
  int num_slice = num_int + 1;

  char outbuff[2048];
  sprintf(outbuff, "Writing igb file %s: ", opts.oigb.c_str());
  PROGRESS<size_t> prg(num_slice, outbuff);

  igb_header igbhead;
  init_igb_header(opts.oigb, igbhead);
  set_igb_header_dim(spacial_data.size(), 1, 1, num_slice, igbhead);
  igbhead.v_inc_t = dt, igbhead.bool_inc_t = true;
  set_igb_header_units("um", "ms", "", igbhead);
  set_igb_header_datatype("float", igbhead);
  write_igb_header(igbhead);

  gettimeofday(&t1, NULL);
  std::vector<std::vector<float> > slice(1);
  slice[0].resize(spacial_data.size());
  mt_real tm = start_time;

  for(int k=0; k < num_slice; k++)
  {
    for(size_t i=0; i<spacial_data.size(); i++)
      slice[0][i] = spacial_data[i] * lookup(tm) * amp;

    write_igb_block(slice, igbhead);
    tm += dt;
    prg.next();
  }
  prg.finish();

  fclose(igbhead.fileptr);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
}





