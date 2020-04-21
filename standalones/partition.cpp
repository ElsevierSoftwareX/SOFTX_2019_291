#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "mt_modes_base.h"
#include "kdpart.hpp"

struct partition_options{
  std::string msh_base;
  std::string size;
};

void print_partition_help()
{
  fprintf(stderr, "partition: compute kdtree based mesh partitioning\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<int>\t (input) number of desired partitions\n", size_par.c_str());
  fprintf(stderr, "\n");
}

int partition_parse_options(int argc, char** argv, struct partition_options & opts)
{
  if(argc < 2) {
    print_partition_help();
    return 1;
  }

  // parse all enclose parameters -----------------------------------------------------------------
  for(int i=1; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh_par, opts.msh_base);
    if(!match) match = parse_param(param, size_par, opts.size);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 3;
    }
  }

  // check if all relevant parameters have been set ---------------------------------------------------
  bool mshok = opts.msh_base.size() > 0, sizeok = opts.size.size() > 0;

  if( !(mshok && sizeok) )
  {
    std::cerr << "Error: Insufficient parameters provided." << std::endl;
    print_partition_help();
    return 2;
  }

  return 0;
}

int main(int argc, char** argv)
{
  struct timeval t1, t2;
  struct partition_options opts;


  int ret = partition_parse_options(argc, argv, opts);
  if (ret != 0) return 1;

  struct mt_meshdata mesh;
  mt_filename meshfile(opts.msh_base, "");

  std::cout << "Reading mesh: " << meshfile.base << std::endl;
  gettimeofday(&t1, NULL);
  read_mesh_selected(mesh, meshfile.format, meshfile.base);
  bucket_sort_offset(mesh.e2n_cnt, mesh.e2n_dsp);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  gettimeofday(&t1, NULL);
  std::cout << "Partitioning .. " << std::endl;
  // compute centerpoints for partitioner
  std::vector<float> xyz(mesh.e2n_cnt.size() * 3);
  for(size_t eidx = 0; eidx < mesh.e2n_cnt.size(); eidx++) {
    vec3r ctr = element_centerpoint(mesh, eidx);
    xyz[eidx*3+0] = ctr.x;
    xyz[eidx*3+1] = ctr.y;
    xyz[eidx*3+2] = ctr.z;
  }

  std::vector<int>  part;
  kdpart::sequential_partitioner<int,float> partitioner;

  int np = atoi(opts.size.c_str());
  partitioner(xyz, np, part);

  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  mt_vector<mt_int> prt; prt.assign(part.begin(), part.end());
  write_vector_ascii(prt, "partitioning.dat");

  return 0;
}
