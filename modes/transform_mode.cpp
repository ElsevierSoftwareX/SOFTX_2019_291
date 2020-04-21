/**
* @file transform_mode.cpp
* @brief Geometric transformation of meshes.
* @author Aurel Neic
* @version 
* @date 2019-03-04
*/


#include "mt_modes_base.h"
#include "dense_mat.hpp"


struct transform_options {
  mt_filename inp_file;
  mt_filename out_file;
  std::string scale;
  std::string translate;
  std::string rotx;
  std::string roty;
  std::string rotz;
};

static const std::string transl_par = "-trans=";
static const std::string rotx_par = "-rotx=";
static const std::string roty_par = "-roty=";
static const std::string rotz_par = "-rotz=";


void print_transform_help()
{
  fprintf(stderr, "transform: Transform the vertex coordinates of a mesh\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the input mesh\n", inp_mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to basename of the output mesh\n", out_mesh_par.c_str());
  fprintf(stderr, "%s<float>\t (optional) Vertex scaling\n", scale_par.c_str());
  fprintf(stderr, "%sx,y,z\t (optional) Vertex translation\n", transl_par.c_str());
  fprintf(stderr, "%s<float>\t (optional) Rotation (in degreees) along the x axis\n", rotx_par.c_str());
  fprintf(stderr, "%s<float>\t (optional) Rotation (in degreees) along the y axis\n", roty_par.c_str());
  fprintf(stderr, "%s<float>\t (optional) Rotation (in degreees) along the z axis\n", rotz_par.c_str());
  fprintf(stderr, "%s<path>\t (optional) format of the input mesh\n", inp_format_par.c_str());
  fprintf(stderr, "%s<path>\t (optional) format of the output mesh\n\n", out_format_par.c_str());
  fprintf(stderr, "\nThe supported input formats are:\n%s\n", input_formats.c_str());
  fprintf(stderr, "\nThe supported output formats are:\n%s\n", output_formats.c_str());
  fprintf(stderr, "\n");
}

int transform_parse_options(int argc, char** argv, struct transform_options & opts)
{
  if(argc < 3) {
    print_transform_help();
    return 1;
  }

  std::string inp_msh_base;
  std::string out_msh_base;
  std::string inp_format;
  std::string out_format;

  for(int i=2; i<argc; i++) {
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, inp_mesh_par, inp_msh_base);
    if(!match) match = parse_param(param, out_mesh_par, out_msh_base);
    if(!match) match = parse_param(param, inp_format_par, inp_format);
    if(!match) match = parse_param(param, out_format_par, out_format);
    if(!match) match = parse_param(param, scale_par, opts.scale);
    if(!match) match = parse_param(param, transl_par, opts.translate);
    if(!match) match = parse_param(param, rotx_par, opts.rotx);
    if(!match) match = parse_param(param, roty_par, opts.roty);
    if(!match) match = parse_param(param, rotz_par, opts.rotz);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 2;
    }
  }

  opts.inp_file.assign(inp_msh_base, inp_format);
  opts.out_file.assign(out_msh_base, out_format);

  if(! (opts.inp_file.isSet() && opts.out_file.isSet()) ) {
    std::cerr << "Mesh transform error: Insufficient parameters provided." << std::endl;
    print_transform_help();
    return 3;
  }

  return 0;
}


void transform_mode(int argc, char** argv)
{
  struct transform_options opts;
  int ret = transform_parse_options(argc, argv, opts);
  if(ret != 0) return;

  struct timeval t1, t2;
  struct mt_meshdata mesh;

  // first read mesh
  std::cout << "Reading mesh: " << opts.inp_file.base << std::endl;
  gettimeofday(&t1, NULL);
  read_mesh_selected(mesh, opts.inp_file.format, opts.inp_file.base);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  // transforming mesh
  bool have_scale     = opts.scale.size();
  bool have_translate = opts.translate.size();
  bool have_rotx      = opts.rotx.size();
  bool have_roty      = opts.roty.size();
  bool have_rotz      = opts.rotz.size();

  if(have_scale || have_translate || have_rotx || have_roty || have_rotz) {
    // first we convert the points we want to work with to vec3r
    mt_vector<vec3r> pts;
    array_to_points(mesh.xyz, pts);

    if(have_rotx || have_roty || have_rotz) {
      mt_real rotx_angle = 0.0, roty_angle = 0.0, rotz_angle = 0.0;

      if(have_rotx) {
        rotx_angle = atof(opts.rotx.c_str());
        rotx_angle *= (MT_PI / 180.0);  // convert degrees angle to radiant angle
      }
      if(have_roty) {
        roty_angle = atof(opts.roty.c_str());
        roty_angle *= (MT_PI / 180.0);  // convert degrees angle to radiant angle
      }
      if(have_rotz) {
        rotz_angle = atof(opts.rotz.c_str());
        rotz_angle *= (MT_PI / 180.0);  // convert degrees angle to radiant angle
      }

      // rotations are applied to a centered object, thus we first compute the centerpoint
      // compute center point
      bbox box; generate_bbox(pts, box);
      vec3r ctr = (box.bounds[0] + box.bounds[1]) * 0.5;
      // now move object into center
      for(vec3r & p : pts) p -= ctr;

      rotate_points(pts, rotx_angle, roty_angle, rotz_angle);

      // now move object back where it was
      for(vec3r & p : pts) p += ctr;
    }

    if(have_translate) {
      mt_vector<std::string> tr_comp; // translation components
      split_string(opts.translate, ',' , tr_comp);
      check_size(tr_comp, 3, "transform_mode");

      mt_real x = atof(tr_comp[0].c_str()), y = atof(tr_comp[1].c_str()),
              z = atof(tr_comp[2].c_str());

      vec3r t(x,y,z);
      for(vec3r & p : pts) p += t;
    }

    if(have_scale) {
      mt_real s = atof(opts.scale.c_str());
      for(vec3r & p : pts) p *= s;
    }

    points_to_array(pts, mesh.xyz);
  }

  // now write mesh in other format
  std::cout << "Writing mesh: " << opts.out_file.base << std::endl;
  gettimeofday(&t1, NULL);
  write_mesh_selected(mesh, opts.out_file.format, opts.out_file.base);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
}
