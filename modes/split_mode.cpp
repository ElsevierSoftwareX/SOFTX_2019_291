/**
* @file split_mode.h
* @brief Meshtool split mode.
* @author Aurel Neic
* @version
* @date 2017-02-13
*/

#include "mt_modes_base.h"

struct split_options {
  mt_filename msh;
  mt_filename outmsh;
  std::string split_name;
};


void print_split_help()
{
  fprintf(stderr, "split: split the mesh connectivity based on a given split file.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the input mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to the split list file\n", split_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to basename of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) input format. (%s)\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "%s<format>\t (optional) output format. (%s)\n", out_format_par.c_str(), output_formats.c_str());
  fprintf(stderr, "\n");
}


int split_parse_options(int argc, char** argv, struct split_options & opts)
{
  if(argc < 3) {
    print_split_help();
    return 1;
  }

  std::string msh_base, outmsh_base;
  std::string ifmt, ofmt;

  for(int i=2; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh_par, msh_base);
    if(!match) match = parse_param(param, outmesh_par, outmsh_base);
    if(!match) match = parse_param(param, inp_format_par, ifmt);
    if(!match) match = parse_param(param, out_format_par, ofmt);
    if(!match) match = parse_param(param, split_par, opts.split_name);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 3;
    }
  }

  opts.msh.assign(msh_base, ifmt);
  opts.outmsh.assign(outmsh_base, ofmt);

  if(opts.msh.isSet() == false || opts.outmsh.isSet() == false || opts.split_name.size() == 0)
  {
    std::cerr << "Split error: Insufficient parameters provided." << std::endl;
    print_split_help();
    return 4;
  }

  return 0;
}


void split_mode(int argc, char** argv)
{
  struct split_options opts;
  int ret = split_parse_options(argc, argv, opts);
  if(ret != 0) return;

  struct timeval t1, t2;

  mt_vector<split_item> splitlist;
  read_split_file(splitlist, opts.split_name);

  // read mesh --------------------------------------------------------------------------------------
  struct mt_meshdata mesh;
  std::cout << "Reading mesh: " << opts.msh.base << std::endl;
  gettimeofday(&t1, NULL);
  read_mesh_selected(mesh, opts.msh.format, opts.msh.base);
  bucket_sort_offset(mesh.e2n_cnt, mesh.e2n_dsp);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  MT_MAP<mt_int,mt_int> new2old;

  for(size_t i=0; i<splitlist.size(); i++) {
    const split_item & se = splitlist[i];

    mt_int eidx = se.elemIdx;

    if(new2old.count(se.newIdx) == 0)
      new2old[se.newIdx] = se.oldIdx;

    const mt_int start = mesh.e2n_dsp[eidx], stop = start + mesh.e2n_cnt[eidx];
    for(mt_int j=start; j<stop; j++) {
      mt_int c = mesh.e2n_con[j];
      if(c == se.oldIdx) {
        mesh.e2n_con[j] = se.newIdx;
        break;
      }
    }
  }

  size_t old_npts = mesh.xyz.size() / 3;
  mesh.xyz.resize( (old_npts + new2old.size()) * 3);

  for(auto it=new2old.begin(); it != new2old.end(); ++it) {
    mt_int newIdx = it->first, oldIdx = it->second;

    mesh.xyz[newIdx*3+0] = mesh.xyz[oldIdx*3+0];
    mesh.xyz[newIdx*3+1] = mesh.xyz[oldIdx*3+1];
    mesh.xyz[newIdx*3+2] = mesh.xyz[oldIdx*3+2];
  }

  std::cout << "Writing mesh: " << opts.outmsh.base << std::endl;
  gettimeofday(&t1, NULL);
  write_mesh_selected(mesh, opts.outmsh.format, opts.outmsh.base);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
}



//      compute_full_mesh_connectivity(surf_start);
//      split_mesh_at_surface(mesh, surf_start);
//      compute_full_mesh_connectivity(mesh);
//      write_mesh_selected(mesh, "carp_txt", "testsplit");
