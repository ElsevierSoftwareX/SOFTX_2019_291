/**
* @file reindex_mode.cpp
* @brief Reindex the mesh vertices
* @author Aurel Neic
* @version 
* @date 2017-09-27
*/

#include "mt_modes_base.h"


struct reindex_options {
  mt_filename msh;
  mt_filename outmsh;
};


void print_reindex_help()
{
  fprintf(stderr, "reindex: reindex a mesh to improve matrix bandwidth and cache efficiency\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the input mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to basename of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "%s<path>\t (optional) format of the input mesh\n", inp_format_par.c_str());
  fprintf(stderr, "%s<path>\t (optional) format of the output mesh\n", out_format_par.c_str());
  fprintf(stderr, "\nThe supported input formats are:\n%s\n", input_formats.c_str());
  fprintf(stderr, "\nThe supported output formats are:\n%s\n", output_formats.c_str());
  fprintf(stderr, "\n");
}


int reindex_parse_options(int argc, char** argv, struct reindex_options & opts)
{
  if(argc < 3) {
    print_reindex_help();
    return 1;
  }

  std::string msh_base;
  std::string outmsh_base;
  std::string inp_format;
  std::string out_format;

  for(int i=2; i<argc; i++) {
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh_par, msh_base);
    if(!match) match = parse_param(param, outmesh_par, outmsh_base);
    if(!match) match = parse_param(param, inp_format_par, inp_format);
    if(!match) match = parse_param(param, out_format_par, out_format);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 2;
    }
  }

  opts.msh.assign(msh_base, inp_format);
  opts.outmsh.assign(outmsh_base, out_format);

  if(! (opts.msh.isSet() && opts.outmsh.isSet()) ) {
    std::cerr << "Mesh reindex error: Insufficient parameters provided." << std::endl;
    print_reindex_help();
    return 3;
  }

  return 0;
}


void reindex_mode(int argc, char** argv)
{
  struct reindex_options opts;
  int ret = reindex_parse_options(argc, argv, opts);
  if(ret != 0) return;

  struct timeval t1, t2;
  struct mt_meshdata mesh;

  // first read mesh
  std::cout << "Reading mesh: " << opts.msh.base << std::endl;
  gettimeofday(&t1, NULL);
  read_mesh_selected(mesh, opts.msh.format, opts.msh.base);
  compute_full_mesh_connectivity(mesh);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  matrixgraph_plotter<mt_vector<mt_int> > plotter(30, 60);
  std::cout << "\nConnectivity graph before reindexing ..\n" << std::endl;
  plotter.print(mesh.n2n_cnt, mesh.n2n_con, '*');
  std::cout << std::endl;

  std::cout << "Reindexing .." << std::endl;
  gettimeofday(&t1, NULL);

  /*
   * We use the Cuthill-McKee algorithm to generate a permutation
   * of the original indices. This permutation is then renumbered
   * either forward or reverse.
   *
   * The permutation is generated as following:
   * Initially, an arbitrary node is selected. It is inserted in the permutation
   * and in the imaginary set R (represented by a boolean vector).
   * Then we loop over each node in perm, adding all connected nodes, that are not yet
   * in R, in ascending order to perm and R.
   *
   */
  size_t nnod = mesh.n2n_cnt.size();
  mt_vector<bool>   inR(nnod, false);
  mt_vector<mt_int> perm(nnod);
  perm.resize(1); perm[0] = 0; inR[0] = true;
  mt_int pidx=0;
  while(perm.size() < nnod)
  {
    mt_int nidx;
    if (pidx < mt_int(perm.size())) nidx = perm[pidx++];
    else {
      int i=0;
      while(inR[i] == true) i++;
      nidx = i;
      std::cerr << "Warning: node " << nidx << " seems decoupled from mesh as it could not be reached by traversing the mesh edges!" << std::endl;
    }
    mt_int start = mesh.n2n_dsp[nidx], stop = start + mesh.n2n_cnt[nidx];
    mt_vector<mt_int> adj(stop - start);

    mt_int widx=0;
    for(mt_int j = start; j<stop; j++) {
      mt_int cidx = mesh.n2n_con[j];
      if( ! inR[cidx] ) {
        adj[widx++] = cidx;
        inR[cidx] = true;
      }
    }

    adj.resize(widx);
    binary_sort(adj);
    perm.append(adj.begin(), adj.end());
  }

  MT_MAP<mt_int, mt_int> nodmap; // map old_idx to new_idx
  // Cuthill-McKee
  // for(size_t i=0; i<perm.size(); i++) nodmap[perm[i]] = i;
  // Reverse Cuthill-McKee
  for(size_t i=0, j=perm.size()-1; i<perm.size(); i++, j--) nodmap[perm[j]] = i;

  // map connectivity to new indexing
  for(size_t i=0; i<mesh.e2n_con.size(); i++) {
    mesh.e2n_con[i] = nodmap[mesh.e2n_con[i]];
  }

  // reorder vertices
  mt_vector<mt_real> xyz_old(mesh.xyz);
  for(auto it = nodmap.begin(); it != nodmap.end(); ++it)
  {
    mt_int oidx = it->first, nidx = it->second;
    mesh.xyz[nidx*3+0] = xyz_old[oidx*3+0];
    mesh.xyz[nidx*3+1] = xyz_old[oidx*3+1];
    mesh.xyz[nidx*3+2] = xyz_old[oidx*3+2];
  }

  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  compute_full_mesh_connectivity(mesh);
  plotter.reset();
  std::cout << "\nConnectivity graph after reindexing ..\n" << std::endl;
  plotter.print(mesh.n2n_cnt, mesh.n2n_con, '*');
  std::cout << std::endl;

  std::cout << "Writing mesh: " << opts.outmsh.base << std::endl;
  gettimeofday(&t1, NULL);
  write_mesh_selected(mesh, opts.outmsh.format, opts.outmsh.base);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
}

