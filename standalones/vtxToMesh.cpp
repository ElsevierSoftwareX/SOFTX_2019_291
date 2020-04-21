/**
 * The main of the enclose executable
 */

#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "mt_modes_base.h"

struct vtxToMesh_options{
  std::string msh_base;
  std::string outmsh_base;
  std::string vtx;
};

void print_vtxToMesh_help()
{
  fprintf(stderr, "vtxToMesh: convert a submesh defined by a set of vertices into a mesh\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to the vtx_file file\n", vtx_par.c_str());
  fprintf(stderr, "%s<path>\t (optional) path to basename of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "\n");
}

int vtxToMesh_parse_options(int argc, char** argv, struct vtxToMesh_options & opts)
{
  if(argc < 2) {
    print_vtxToMesh_help();
    return 1;
  }

  // parse all enclose parameters -----------------------------------------------------------------
  for(int i=1; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh_par, opts.msh_base);
    if(!match) match = parse_param(param, vtx_par, opts.vtx);
    if(!match) match = parse_param(param, outmesh_par, opts.outmsh_base);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 3;
    }
  }
  fixBasename(opts.msh_base);
  fixBasename(opts.outmsh_base);

  // check if all relevant parameters have been set ---------------------------------------------------
  bool mshok = opts.msh_base.size() > 0, vtxok = opts.vtx.size() > 0;

  if( !(mshok && vtxok) )
  {
    std::cerr << "Error: Insufficient parameters provided." << std::endl;
    print_vtxToMesh_help();
    return 2;
  }

  return 0;
}

int main(int argc, char** argv)
{
  struct timeval t1, t2;
  struct vtxToMesh_options opts;
  struct mt_meshdata mesh, surf;

  int ret = vtxToMesh_parse_options(argc, argv, opts);
  if (ret != 0) return 1;

  std::cout << "Reading mesh: " << opts.msh_base << std::endl;
  gettimeofday(&t1, NULL);
  readElements_general(mesh, opts.msh_base);
  readPoints_general(mesh.xyz, opts.msh_base);
  readFibers_general(mesh.lon, mesh.e2n_cnt.size(), opts.msh_base);

  transpose_connectivity(mesh.e2n_cnt, mesh.e2n_con, mesh.n2e_cnt, mesh.n2e_con);
  mesh.e2n_dsp.resize(mesh.e2n_cnt.size());
  mesh.n2e_dsp.resize(mesh.n2e_cnt.size());
  bucket_sort_offset(mesh.e2n_cnt, mesh.e2n_dsp);
  bucket_sort_offset(mesh.n2e_cnt, mesh.n2e_dsp);

  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  mt_vector<mt_int> jun_vtx;
  std::cout << "Reading junction vertices: " << opts.vtx << std::endl;
  read_vtx(jun_vtx, opts.vtx);

  MT_USET<mt_int> nodes, elems;

  std::cout << "Computing mesh .." << std::endl;
  gettimeofday(&t1, NULL);

  nodes.insert(jun_vtx.begin(), jun_vtx.end());
  nodeSet_to_elemSet(mesh, nodes, elems);

  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  mt_vector<bool> keep(mesh.e2n_cnt.size(), false);
  for(MT_USET<mt_int>::iterator it = elems.begin(); it != elems.end(); ++it)
  {
    keep[*it] = true;
    printf("EIDX: %ld\n", *it);
  }

  if(opts.outmsh_base.size() > 0)
  {
    mt_vector<mt_int> nod, eidx;
    restrict_meshdata(keep, mesh, nod, eidx);

    std::cout << "Writing mesh: " << opts.outmsh_base << std::endl;
    gettimeofday(&t1, NULL);
    writeElements(mesh, opts.outmsh_base + CARPTXT_ELEM_EXT);
    writePoints(mesh.xyz, opts.outmsh_base + CARPTXT_PTS_EXT);
    writeFibers(mesh.lon, mesh.e2n_cnt.size(), opts.outmsh_base + CARPTXT_LON_EXT);
    binary_write(nod.begin(), nod.end(), opts.outmsh_base + NOD_EXT);
    binary_write(eidx.begin(), eidx.end(), opts.outmsh_base + EIDX_EXT);
    gettimeofday(&t2, NULL);
    std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
  }

  return 0;
}
