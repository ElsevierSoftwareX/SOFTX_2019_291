/**
* @file correspondance.cpp
* @brief Standalone for computing the correspondance between the nodes of two meshes.
* @author Aurel Neic
* @version 
* @date 2017-09-12
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <string>

#include "mt_modes_base.h"


struct corresp_options{
  std::string msh1_base;
  std::string msh2_base;
  std::string ifmt;
  std::string out;
  std::string ofmt;
};

void print_corresp_help()
{
  fprintf(stderr, "correspance: generate the correspondance between the vertices of two meshes.");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of mesh 1.\n", mesh1_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to basename of mesh 2.\n", mesh2_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh input format. may be: %s\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "%s<path>\t (optional) output file name; default: 'corr.`ofmt`\n", out_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) output file format; choices: txt,bin; default: txt\n", out_format_par.c_str());
  fprintf(stderr, "\n");
  fprintf(stderr, "Output formats:\n");
  fprintf(stderr, " -) `txt` format:\n");
  fprintf(stderr, "    n                     \n");
  fprintf(stderr, "    0:idx[0]:dist[0]      \n");
  fprintf(stderr, "    1:idx[1]:dist[1]      \n");
  fprintf(stderr, "    ...                   \n");
  fprintf(stderr, "    n-2:idx[n-2]:dist[n-2]\n");
  fprintf(stderr, "    n-1:idx[n-1]:dist[n-1]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " -) `bin` format:\n");
  fprintf(stderr, "    n      one unsigned integer of size %ldbytes\n", sizeof(size_t));
  fprintf(stderr, "    idx    block of `n` signed integers of size %ldbytes each\n", sizeof(mt_int));
  fprintf(stderr, "    dist   block of `n` reals of size %ldbytes each\n", sizeof(mt_real));
  fprintf(stderr, "\n");
  fprintf(stderr, "Where `n` is the number of nodes in mesh 1, `idx[i]` is\n");
  fprintf(stderr, "the index of the point in mesh 2 which is closest to the\n");
  fprintf(stderr, "point with index `i` in mesh 1. The corresponding distance\n");
  fprintf(stderr, "between the two points is given by `dist[i].`\n");
  fprintf(stderr, "\n");
}

int corresp_parse_options(int argc, char** argv, struct corresp_options & opts)
{
  if(argc < 2)
    return 1;

  // parse all parameters ---------------------------------------------------------------
  for(int i=1; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh1_par, opts.msh1_base);
    if(!match) match = parse_param(param, mesh2_par, opts.msh2_base);
    if(!match) match = parse_param(param, inp_format_par, opts.ifmt);
    if(!match) match = parse_param(param, out_par, opts.out);
    if(!match) match = parse_param(param, out_format_par, opts.ofmt);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 3;
    }
  }
  fixBasename(opts.msh1_base);
  fixBasename(opts.msh2_base);

  // check if all relevant parameters have been set -------------------------------------
  const bool msh1_ok = opts.msh1_base.size() > 0;
  const bool msh2_ok = opts.msh2_base.size() > 0;
  if( !(msh1_ok && msh2_ok) )
  {
    std::cerr << "Error: Insufficient parameters provided." << std::endl;
    return 2;
  }

  // check output file format -----------------------------------------------------------
  if (opts.ofmt.empty())
    opts.ofmt = "txt";
  if ((opts.ofmt != "txt") && (opts.ofmt != "bin"))
    return 4;

  // check output file name -------------------------------------------------------------
  if (opts.out.empty())
    opts.out = "corr." + opts.ofmt;
 
  return 0;
}


int main(int argc, char** argv)
{
  struct timeval t1, t2;
  struct corresp_options opts;
  struct mt_meshdata mesh1, mesh2;

  int ret = corresp_parse_options(argc, argv, opts);
  if (ret != 0) {
    print_corresp_help();
    return 1;
  }

  std::cout << "Reading primary mesh: " << opts.msh1_base << std::endl;
  gettimeofday(&t1, NULL);
  read_mesh_selected(mesh1, opts.ifmt, opts.msh1_base, CRP_READ_ELEM | CRP_READ_PTS);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  std::cout << "Reading secondary mesh: " << opts.msh2_base << std::endl;
  gettimeofday(&t1, NULL);
  read_mesh_selected(mesh2, opts.ifmt, opts.msh2_base, CRP_READ_ELEM | CRP_READ_PTS);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  mt_vector<mt_int>  corr;
  mt_vector<mt_real> corr_dist;
  std::cout << "Computing correspondance .. " << std::endl;
  gettimeofday(&t1, NULL);
  compute_correspondance(mesh1, mesh2, corr, corr_dist);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  const char* corr_file = opts.out.c_str();
  std::cout << "Writing correspondance " << corr_file << ", format: " << opts.ofmt << std::endl;
  if (opts.ofmt == "txt") {
    FILE* fd = fopen(corr_file, MT_FOPEN_WRITE);  
    fprintf(fd, "%lu\n", corr.size());
    for(size_t i=0; i<corr.size(); i++)
      fprintf(fd, "%ld:%ld:%.3f\n", i, corr[i], sqrt(corr_dist[i]));
    fclose(fd);
  }
  else if (opts.ofmt == "bin") {
    FILE* fd = fopen(corr_file, "wb");
    size_t size = corr.size();
    fwrite(&size, 1, sizeof(size_t), fd);
    fwrite(corr.data(), size, sizeof(mt_int), fd);
    fwrite(corr_dist.data(), size, sizeof(mt_real), fd);
    fclose(fd);
  }

  return 0;
}
