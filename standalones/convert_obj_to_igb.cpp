/**
* @file convert_obj_to_igb.cpp
* @brief Standalone for converting a moving mesh stored as a series of .obj files to a .vtk basemesh + an .igb file storing the displacements
* @author Jana Fuchsberger
* @version
* @date 2019-20-05
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <string>

#include "mt_modes_base.h"

static const std::string obj_base_par = "-obj_base=";
static const std::string first_idx_par = "-first_idx=";
static const std::string last_idx_par = "-last_idx=";
static const std::string end_time_par = "-Tend=";

struct obj_to_igb_options{
  std::string obj_base;
  std::string first_idx;
  std::string last_idx;
  std::string tend_str;
};

void print_obj_to_igb_help()
{
  fprintf(stderr, "convert_obj_igb: Convert a moving mesh stored as a series of .obj files to a .vtk basemesh + an .igb file storing the displacements.");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of .obj series.\n", obj_base_par.c_str());
  fprintf(stderr, "%s<int>\t (input) first index of .obj series.\n",first_idx_par.c_str());
  fprintf(stderr, "%s<int>\t (input) last index of .obj series.\n", last_idx_par.c_str());
  fprintf(stderr, "%s<float>\t (optional) End time of movement in seconds.\n", end_time_par.c_str());
  fprintf(stderr, "\n");
}

int obj_to_igb_parse_options(int argc, char** argv, struct obj_to_igb_options & opts)
{
  if(argc < 3) {
    print_obj_to_igb_help();
    return 1;
  }

  // parse all parameters -----------------------------------------------------------------
  for(int i=1; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, obj_base_par, opts.obj_base);
    if(!match) match = parse_param(param, first_idx_par, opts.first_idx);
    if(!match) match = parse_param(param, last_idx_par, opts.last_idx);
    if(!match) match = parse_param(param, end_time_par, opts.tend_str);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 3;
    }
  }

  // check if all relevant parameters have been set ---------------------------------------------------
  bool obj_base_ok = opts.obj_base.size() > 0, first_idx_ok = opts.first_idx.size() > 0, last_idx_ok = opts.last_idx.size() > 0;

  if(!(obj_base_ok&&first_idx_ok&&last_idx_ok))
  {
    std::cerr << "Error: Insufficient parameters provided." << std::endl;
    print_obj_to_igb_help();
    return 2;
  }

  return 0;
}


int main(int argc, char** argv)
{
  struct timeval t1, t2;
  struct obj_to_igb_options opts;

  int ret = obj_to_igb_parse_options(argc, argv, opts);
  if (ret != 0) return 1;

  const int first_idx  = atoi(opts.first_idx.c_str());
  const int last_idx  = atoi(opts.last_idx.c_str());
  const float tend = opts.tend_str.size() > 0 ? atof(opts.tend_str.c_str()) : 1.0;
  const float dt = tend / (float)(last_idx-first_idx);

  gettimeofday(&t1, NULL);

  std::cout << "Creating base mesh from first .obj file." << std::endl;

  char base_filename[300];
  sprintf(base_filename, "%s_%06d%s", opts.obj_base.c_str(), first_idx, OBJ_EXT);
  mt_meshdata base_mesh;
  mt_vector<mt_int> base_nrml_idx;
  mt_vector<mt_real> base_nrml_xyz;

  read_obj_file(base_mesh, base_nrml_idx, base_nrml_xyz, std::string(base_filename));

  writeVTKmesh(base_mesh, opts.obj_base+VTK_EXT);

  PROGRESS<mt_int> progress(last_idx-first_idx, "Converting .obj series to displacement .igb file:");

  igb_header igb;
  init_igb_header(opts.obj_base+IGB_EXT, igb);
  set_igb_header_datatype("vec3f", igb);
  igb.v_x = base_mesh.xyz.size() / 3;
  igb.v_y = 1;
  igb.v_z = 1;
  igb.v_t = last_idx - first_idx + 1;
  igb.v_inc_t = dt;
  igb.bool_x = true;
  igb.bool_y = true;
  igb.bool_z = true;
  igb.bool_t = true;
  igb.bool_type = true;
  std::string unitx = "m";
  std::string unitt = "s";
  strncpy(igb.v_unites_x, unitx.c_str(), 40);
  igb.bool_unites_x = true;
  strncpy(igb.v_unites_y, unitx.c_str(), 40);
  igb.bool_unites_y = true;
  strncpy(igb.v_unites_z, unitx.c_str(), 40);
  igb.bool_unites_z = true;
  strncpy(igb.v_unites_t, unitt.c_str(), 40);
  igb.bool_unites_t = true;
  igb.v_org_t = 0.0f;
  igb.bool_org_t = true;
  igb.v_dim_t = tend;
  igb.bool_dim_t = true;
  igb.v_facteur = 1.0f;
  igb.bool_facteur = true;
  igb.v_systeme = IGB_LITTLE_ENDIAN;
  write_igb_header(igb);

  std::vector<std::vector<float> > disp(1);
  disp[0].assign(base_mesh.xyz.size(),0.0);
  write_igb_block(disp, igb);

  for(int idx=first_idx+1; idx<=last_idx; idx++)
  {
    progress.next();
    char filename[300];
    sprintf(filename, "%s_%06d%s", opts.obj_base.c_str(), idx, OBJ_EXT);
    mt_meshdata mesh;
    mt_vector<mt_int> nrml_idx;
    mt_vector<mt_real> nrml_xyz;

    read_obj_file(mesh, nrml_idx, nrml_xyz, std::string(filename));

    for(unsigned int j=0; j<mesh.xyz.size()/3; j++)
    {
      disp[0][j*3  ]=mesh.xyz[j*3  ]-base_mesh.xyz[j*3  ];
      disp[0][j*3+1]=mesh.xyz[j*3+1]-base_mesh.xyz[j*3+1];
      disp[0][j*3+2]=mesh.xyz[j*3+2]-base_mesh.xyz[j*3+2];
    }
    write_igb_block(disp, igb);
  }
  progress.finish();

  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  return 0;
}
