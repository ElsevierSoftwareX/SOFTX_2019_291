/**
* @file convert_mode.cpp
* @brief Meshtool convert mode.
* @author Aurel Neic
* @version
* @date 2017-02-13
*/

#include "mt_modes_base.h"

struct convert_options {
  mt_filename inp_file;
  mt_filename out_file;
  std::string scale;
};


void print_convert_help()
{
  fprintf(stderr, "convert: convert between different mesh formats\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the input mesh\n", inp_mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) format of the input mesh\n", inp_format_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to basename of the output mesh\n", out_mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) format of the output mesh\n", out_format_par.c_str());
  fprintf(stderr, "%s<float>\t (optional) Vertex scaling applied to output mesh\n", scale_par.c_str());
  fprintf(stderr, "\nThe supported input formats are:\n%s\n", input_formats.c_str());
  fprintf(stderr, "\nThe supported output formats are:\n%s\n", output_formats.c_str());
  fprintf(stderr, "\n");
}


int convert_parse_options(int argc, char** argv, struct convert_options & opts)
{
  if(argc < 3) {
    print_convert_help();
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

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 2;
    }
  }

  opts.inp_file.assign(inp_msh_base, inp_format);
  opts.out_file.assign(out_msh_base, out_format);

  if(! (opts.inp_file.isSet() && opts.out_file.isSet()) ) {
    std::cerr << "Mesh convert error: Insufficient parameters provided." << std::endl;
    print_convert_help();
    return 3;
  }

  return 0;
}


void convert_mode(int argc, char** argv)
{
  struct convert_options opts;
  int ret = convert_parse_options(argc, argv, opts);
  if(ret != 0) return;

  struct timeval t1, t2;
  struct mt_meshdata mesh;

  // first read mesh
  std::cout << "Reading mesh: " << opts.inp_file.base << std::endl;
  gettimeofday(&t1, NULL);
  read_mesh_selected(mesh, opts.inp_file.format, opts.inp_file.base);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  // we correct inside-out tets since some 3rd party tools dont keep the orientation
  gettimeofday(&t1, NULL);
  std::cout << "Processing mesh .." << std::endl;
  correct_insideOut(mesh);

  if(opts.scale.size()) {
    float scale = atof(opts.scale.c_str());
    for(size_t i=0; i<mesh.xyz.size(); i++) mesh.xyz[i] *= scale;
  }

  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  // now write mesh in other format
  std::cout << "Writing mesh: " << opts.out_file.base << std::endl;
  gettimeofday(&t1, NULL);
  write_mesh_selected(mesh, opts.out_file.format, opts.out_file.base);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
}
