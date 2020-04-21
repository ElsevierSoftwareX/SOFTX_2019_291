/**
 * The main of the enclose executable
 */


#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "mt_modes_base.h"

struct enclose_options{
  std::string msh_base;
  std::string submsh_base;
  std::string surf;
  std::string thr;
};

void print_enclose_help()
{
  fprintf(stderr, "enclose: extract mesh and surface of the volume enclosing a given surface\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to the surface file\n", surf_par.c_str());
  fprintf(stderr, "%s<int>\t (input) threshold defining the number of layers of the enclosure\n", thr_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to basename of the mesh to query\n", submesh_par.c_str());
  fprintf(stderr, "\n");
}

int enclose_parse_options(int argc, char** argv, struct enclose_options & opts)
{
  if(argc < 2) {
    print_enclose_help();
    return 1;
  }

  // parse all enclose parameters -----------------------------------------------------------------
  for(int i=1; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh_par, opts.msh_base);
    if(!match) match = parse_param(param, submesh_par, opts.submsh_base);
    if(!match) match = parse_param(param, surf_par, opts.surf);
    if(!match) match = parse_param(param, thr_par, opts.thr);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 3;
    }
  }
  fixBasename(opts.msh_base);
  fixBasename(opts.submsh_base);

  // check if all relevant parameters have been set ---------------------------------------------------
  bool mshok = opts.msh_base.size() > 0, submshok = opts.submsh_base.size() > 0,
       surfok = opts.surf.size() > 0, throk = opts.thr.size() > 0;

  if( !(mshok && submshok && surfok && throk) )
  {
    std::cerr << "Error: Insufficient parameters provided." << std::endl;
    print_enclose_help();
    return 2;
  }

  return 0;
}

int main(int argc, char** argv)
{
  struct timeval t1, t2;
  struct enclose_options opts;
  struct mt_meshdata mesh, innersurf;

  int ret = enclose_parse_options(argc, argv, opts);
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

  std::cout << "Reading surface: " << opts.surf << std::endl;
  gettimeofday(&t1, NULL);
  readElements(innersurf, opts.surf);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  int num_layers = atoi(opts.thr.c_str());
  MT_USET<mt_int> nodes, elems;

  std::cout << "Computing " << num_layers << " layer enclosure .." << std::endl;
  gettimeofday(&t1, NULL);

  mt_vector<mt_int> surf_vtx(innersurf.e2n_con);
  binary_sort(surf_vtx); unique_resize(surf_vtx);

  nodes.insert(surf_vtx.begin(), surf_vtx.end());
  nodeSet_to_elemSet(mesh, nodes, elems);

  for(int i=1; i<num_layers; i++)
  {
    elemSet_to_nodeSet(mesh, elems, nodes);
    nodeSet_to_elemSet(mesh, nodes, elems);
  }

  mt_vector<bool> keep(mesh.e2n_cnt.size(), false);
  for(MT_USET<mt_int>::iterator it = elems.begin(); it != elems.end(); ++it)
  {
    keep[*it] = true;
  }
  mt_vector<mt_int> nod, eidx;
  restrict_meshdata(keep, mesh, nod, eidx);

  // we also map the enclosed surface to the new indexing
  map_glob2loc(nod, innersurf.e2n_con);

  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  std::cout << "Writing mesh: " << opts.submsh_base << std::endl;
  gettimeofday(&t1, NULL);
  writeElements(mesh, opts.submsh_base + CARPTXT_ELEM_EXT);
  writePoints(mesh.xyz, opts.submsh_base + CARPTXT_PTS_EXT);
  writeFibers(mesh.lon, mesh.e2n_cnt.size(), opts.submsh_base + CARPTXT_LON_EXT);
  binary_write(nod.begin(), nod.end(), opts.submsh_base + NOD_EXT);
  binary_write(eidx.begin(), eidx.end(), opts.submsh_base + EIDX_EXT);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  std::cout << "Computing surface .. " << std::endl;
  gettimeofday(&t1, NULL);
  mt_vector<mt_int> surf_nod;
  mt_meshdata meshsurf;
  compute_surface(mesh.etype, mesh.e2n_cnt, mesh.e2n_con, meshsurf);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  // write surface and vtx of new meshsurf
  std::string wname = opts.submsh_base + SURF_EXT;
  std::cout << "Writing " << wname << std::endl;
  write_surf(meshsurf.e2n_cnt, meshsurf.e2n_con, wname);

  surf_nod.assign(meshsurf.e2n_con.begin(), meshsurf.e2n_con.end());
  binary_sort(surf_nod); unique_resize(surf_nod);
  wname = opts.submsh_base + SURF_EXT + VTX_EXT;
  std::cout << "Writing " << wname << std::endl;
  write_vtx(surf_nod, wname);

  // write surface and vtx of enclosed surf mapped to new idx range
  wname = opts.submsh_base + ".inner" + SURF_EXT;
  std::cout << "Writing " << wname << std::endl;
  write_surf(innersurf.e2n_cnt, innersurf.e2n_con, wname);

  surf_nod.assign(innersurf.e2n_con.begin(), innersurf.e2n_con.end());
  binary_sort(surf_nod); unique_resize(surf_nod);

  wname = opts.submsh_base + ".inner" + SURF_EXT + VTX_EXT;
  std::cout << "Writing " << wname << std::endl;
  write_vtx(surf_nod, wname);

  return 0;
}
