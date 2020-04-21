/**
* @file refine_igb_time.cpp
* @brief Refine an igb file in time by upsampling it with a spline interpolation
* @author Elias Karabelas
* @version
* @date 2017-10-28
*/


#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>

#include "spline.h"
#include "mt_utils.h"

static const std::string mesh_input_dynpts   = "-input_igb=";
static const std::string desired_delta_t     = "-dt=";
static const std::string mesh_output_dynpts  = "-output_igb=";
static const std::string eval_deriv          = "-calc_derivative=";
static const std::string clamp_values        = "-clamp_values=";

struct refine_dynpts_options
{
  std::string msh_input_dynpts_file;
  std::string msh_output_dynpts_file;
  std::string des_dt;
  std::string calc_deriv;
  std::string clamp;
};

void print_refine_dynpts_help()
 {
  fprintf(stderr, "refine_igb_time: refine time sampling of an IGB file. This is achieved by using a spline interpolation in time for each point.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to input dynpts file\n", mesh_input_dynpts.c_str());
  fprintf(stderr, "%s<float>\t (input) new delta_t (needs to be in the same unit as input)\n", desired_delta_t.c_str());
  fprintf(stderr, "%s<path>\t (output) path to output dynpts file\n", mesh_output_dynpts.c_str());
  fprintf(stderr, "%s<bool>\t (optional) evaluate derivative of upsampled file\n", eval_deriv.c_str());
  fprintf(stderr, "%s<list>\t (optional) clamp values to defined range. The list is dpn1_low,dpn1_high:dpn2_low,dpn2_high:...\n\t\t\tClamp values have to be defined for each dpn, no clamping is indicated by inf", clamp_values.c_str());
  fprintf(stderr, "\n");
}

int refine_dynpts_parse_options(int argc, char** argv, struct refine_dynpts_options & opts)
{
  if(argc < 3) {
    print_refine_dynpts_help();
    return 1;
  }

  // parse all enclose parameters -----------------------------------------------------------------
  for(int i=1; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh_input_dynpts, opts.msh_input_dynpts_file);
    if(!match) match = parse_param(param, desired_delta_t, opts.des_dt);
    if(!match) match = parse_param(param, mesh_output_dynpts, opts.msh_output_dynpts_file);
    if(!match) match = parse_param(param, eval_deriv, opts.calc_deriv);
    if(!match) match = parse_param(param, clamp_values, opts.clamp);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 3;
    }
  }


  // check if all relevant parameters have been set ---------------------------------------------------
  bool dynptsinputok = opts.msh_input_dynpts_file.size() > 0,
       dynptstargetok = opts.msh_output_dynpts_file.size() > 0,
       dtok = opts.des_dt.size() > 0;
  if( !(dynptsinputok && dynptstargetok && dtok) )
  {
    std::cerr << "Error: Insufficient parameters provided." << std::endl;
    print_refine_dynpts_help();
    return 3;
  }

  return 0;
}

int main(int argc, char** argv)
{
  struct timeval t1, t2;
  struct refine_dynpts_options opts;
  igb_header igb_head_from;
  igb_header igb_head_to;
  igb_header igb_head_to_deriv;
  int ret = refine_dynpts_parse_options(argc, argv, opts);
  if(ret != 0) return 1;

  mt_real newdt = static_cast<mt_real>(atof(opts.des_dt.c_str()));
  bool deriv = opts.calc_deriv.size() > 0;

  short data_idx = -1;

  setup_data_format(opts.msh_input_dynpts_file, opts.msh_output_dynpts_file, data_idx,
                    igb_head_from, igb_head_to);

  //init_igb_header(opts.msh_input_dynpts_file, igb_head_from);
  //read_igb_header(igb_head_from);


  //adjust dt to really reach endpoint
  const mt_real endpoint     = static_cast<mt_real>(igb_head_from.v_dim_t);
  const mt_real dt_orig      = static_cast<mt_real>(igb_head_from.v_dim_t / (igb_head_from.v_t -1));
  const mt_real newdt_scaled = endpoint / std::floor(endpoint / newdt);
  const mt_int tsteps_orig   = static_cast<mt_int>(igb_head_from.v_t) - 1; //v_t counts number of TIME points. tsteps is one less
  const mt_int tsteps_new    = static_cast<mt_int>(std::floor(endpoint / newdt));
  const mt_int num_points    = static_cast<mt_int>(igb_head_from.v_x);
  mt_int dpn = 0;
  switch(data_idx) {
    case 1:
    {
      dpn = 1;
      break;
    }
    case 3:
    {
      dpn = 3;
      break;
    }
    case 4:
    {
      dpn = 9;
      break;
    }
    default:
    {
      fprintf(stderr, "ERROR in %s: Invalid data_type %d in IGB file %s!\n", __func__, data_idx, opts.msh_input_dynpts_file.c_str());
      exit(1);
    }
  }

  // Figure out clamping
  bool do_clamp = opts.clamp.size() > 0;
  mt_vector<float> clamp_vals(2*dpn);
  for(size_t k=0; k < clamp_vals.size(); k++) {
   if(k % 2) // False if even number
    clamp_vals[k] = std::numeric_limits<float>::max();
   else
    clamp_vals[k] = std::numeric_limits<float>::lowest();
  }

  if(do_clamp) {
   mt_vector<std::string> split_str;
   split_string(opts.clamp, ':', split_str);
   assert(split_str.size() == dpn);
   for(size_t k=0; k < split_str.size(); k++) {
    mt_vector<std::string> tmp;
    split_string(split_str[k], ',', tmp);
    assert(tmp.size() == 2);
    for(short t=0; t < 2; t++)
     if(tmp[t] != "inf")
      clamp_vals[2*k + t] = atof(tmp[t].c_str());
   }
  }

  std::vector<std::vector<float> > dynpts_data_input;
  std::vector<std::vector<float> > dynpts_data_output;

  //for creating the splines
  mt_vector<spline*> spline_vec(dpn * num_points);

  //Read the input dynpts file----------------------------------------------------
  std::cout << "Reading " << opts.msh_input_dynpts_file << std::endl;
  gettimeofday(&t1, NULL);
  read_igb_data(dynpts_data_input, igb_head_from, opts.msh_input_dynpts_file);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
  //-----------------------------------------------------------------------------
  std::stringstream out;
  //build the vectors
  std::cout << "Refine igb file " << opts.msh_input_dynpts_file << " and output to igb file " << opts.msh_output_dynpts_file << " :" << std::endl;
  gettimeofday(&t1, NULL);
  dynpts_data_output.resize(1);
  dynpts_data_output[0].resize(dpn*num_points);

  // FILL THE SPLINES one for each point in the mesh
  PROGRESS<mt_int> progress_splines(num_points, "SPLINE GENERATION: ");
  #ifdef OPENMP
  #pragma omp parallel
  #endif
  {
   mt_vector<mt_vector<mt_real> > X(dpn);
   for(mt_int s=0; s < dpn; s++)
    X[s].resize(tsteps_orig + 1);
   #ifdef OPENMP
   #pragma omp for schedule(static)
   #endif
    for(mt_int nidx = 0; nidx < num_points; nidx++)
    {
     for(mt_int s=0; s < dpn; s++) {
      for(mt_int tidx=0; tidx < tsteps_orig + 1; tidx++) {
       X[s][tidx] = static_cast<mt_real>(dynpts_data_input[tidx][nidx * dpn + s]);
      }
      spline_vec[dpn*nidx+s] = new spline(X[s].begin(), X[s].end(), 0.0, dt_orig, 0.0, 0.0);
     }
     #ifdef OPENMP
     #pragma omp critical
     #endif
     progress_splines.next();
    }
  }
  progress_splines.finish();

  //splines are filled now we have to calculate the interpolated values
  PROGRESS<mt_int> progress_upsample(tsteps_new+1, "UPSAMPLING: ");
  //adapt igb header file
  igb_head_to.fileptr  = NULL;
  igb_head_to.v_t      = tsteps_new + 1;
  igb_head_to.v_inc_t  = newdt_scaled;
  write_igb_header(igb_head_to);
  for(mt_int s=0; s < tsteps_new + 1; s++) {
   const mt_real actT = s * newdt_scaled;
   #ifdef OPENMP
   #pragma omp parallel for schedule(static)
   #endif
   for(mt_int nidx = 0; nidx < num_points; nidx++) {
    for(mt_int k=0; k < dpn; k++) {
     dynpts_data_output[0][nidx*dpn+k] = clamp(static_cast<float>((*spline_vec[nidx*dpn+k])(actT)), clamp_vals[2*k+0], clamp_vals[2*k+1], true);
    }
   }
   write_igb_block(dynpts_data_output, igb_head_to);
   progress_upsample.next();
  }
  progress_upsample.finish();

  if(deriv) {
   PROGRESS<mt_int> progress_derivative(tsteps_new+1, "EVAL DERIVATIVE: ");
   //adapt igb header file
   igb_head_to_deriv = igb_head_to;
   igb_head_to_deriv.fileptr  = NULL;
   igb_head_to_deriv.v_t      = tsteps_new + 1;
   igb_head_to_deriv.v_inc_t  = newdt_scaled;
   igb_head_to_deriv.filename = "derivative.igb";
   write_igb_header(igb_head_to_deriv);
   for(mt_int s=0; s < tsteps_new + 1; s++) {
    const mt_real actT = s * newdt_scaled;
    #ifdef OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for(mt_int nidx = 0; nidx < num_points; nidx++) {
     for(mt_int k=0; k < dpn; k++)
       dynpts_data_output[0][nidx*dpn+k] = static_cast<float>(spline_vec[nidx*dpn+k]->prime(actT));
    }
    write_igb_block(dynpts_data_output, igb_head_to_deriv);
    progress_derivative.next();
   }
   progress_derivative.finish();
   fclose(igb_head_to_deriv.fileptr);
  }
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
  fclose(igb_head_to.fileptr);
  for(size_t s=0; s < spline_vec.size(); s++)
   delete spline_vec[s];

  return 0;
}
