/**
* @file map_mode.h
* @brief Meshtool map mode.
* @author Andrew Crozier, Aurel Neic.
* @version
* @date 2017-02-13
*/

#include "mt_modes_base.h"

struct map_options {
  std::string submsh_base;
  std::string files;
  std::string outdir;
  std::string mode;
  std::string info;
};

void print_map_help()
{
  fprintf(stderr, "map: map .vtx, .surf and .neubc files to the indexing of a submesh (obtained from 'meshtool extract mesh')\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the submesh\n", submesh_par.c_str());
  fprintf(stderr, "%spatterns\t (input) files to map to submesh\n", file_par.c_str());
  fprintf(stderr, "%s<path>\t (output) directory to place mapped files\n", outdir_par.c_str());
  fprintf(stderr, "%s<m2s|s2m>\t (optional) Map mesh-to-submesh (m2s) or vice-versa (s2m). "
                  "Default is m2s.\n", mode_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to the information file\n", info_par.c_str());
  fprintf(stderr, "\n");
  fprintf(stderr, "'patterns' is a comma separated list of files or glob patterns, e.g.:\n");
  fprintf(stderr, "    mesh/endo.vtx,mesh/endo.surf,mesh/*.neubc\n");
  fprintf(stderr, "Any mapped files that would overwrite their source file are skipped.\n");
  fprintf(stderr, "\n");
}

int map_parse_options(int argc, char** argv, struct map_options & opts)
{
  if (argc < 3) {
    print_map_help();
    return 1;
  }

  for (int i=2; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, submesh_par, opts.submsh_base);
    if(!match) match = parse_param(param, file_par,    opts.files);
    if(!match) match = parse_param(param, outdir_par,  opts.outdir);
    if(!match) match = parse_param(param, mode_par,  opts.mode);
    if(!match) match = parse_param(param, info_par,  opts.info);
    fixBasename(opts.submsh_base);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 2;
    }
  }
  if (! (opts.submsh_base.size() > 0 &&
         opts.files.size()       > 0 &&
         opts.outdir.size()      > 0 ) ) {
    std::cerr << "Mesh map error: Insufficient parameters provided.\n" << std::endl;
    print_map_help();
    return 3;
  }
  if (opts.info.empty())  
    opts.info = opts.outdir + "/map_info.txt";
  return 0;
}


void map_mode(int argc, char** argv)
{
  struct map_options opts;
  if (map_parse_options(argc, argv, opts) != 0)
    return;

  struct timeval t1, t2;

  // start program timer
  gettimeofday(&t1, NULL);

  // Mappings
  mt_vector<mt_int> nod;
  mt_vector<mt_int> eidx;

  // Expand file list
  std::list<std::string> vtx_files;
  std::list<std::string> surf_files;
  std::list<std::string> nbc_files;
  std::list<std::string>::iterator fit;

  mt_vector<std::string> patternlist;
  split_string(opts.files, ',', patternlist);

  for (size_t i=0; i<patternlist.size(); i++) {
    mt_vector<std::string> matches;
    #ifdef WINBUILD
    matches.push_back(patternlist[i]);
    #else
    mt_glob(patternlist[i], matches);
    #endif

    for (size_t j=0; j<matches.size(); j++) {
      std::string outfile = opts.outdir + "/" + mt_basename(matches[j]);
      if (matches[j].compare(outfile) == 0)
        // Output file same as input
        std::cout << "Output file " << outfile << " would overwrite the source file, skipping." << std::endl;
      else if (endswith(matches[j], ".vtx"))
        vtx_files.push_back(matches[j]);
      else if (endswith(matches[j], ".surf"))
        surf_files.push_back(matches[j]);
      else if (endswith(matches[j], ".neubc"))
        nbc_files.push_back(matches[j]);
      else
        std::cout << "Matched file " << matches[j] << " not of supported type, skipping." << std::endl;
    }
  }

  bool map_m2s = true;
  if(opts.mode.size()) {
    if(opts.mode.compare("m2s") == 0) map_m2s = true;
    if(opts.mode.compare("s2m") == 0) map_m2s = false;
  }

  MT_MAP<mt_int, mt_int> g2l, g2l_ele;

  // Read node mapping
  std::cout << "Reading node map " << opts.submsh_base + NOD_EXT << std::endl;
  binary_read(nod, opts.submsh_base + NOD_EXT);

  if(map_m2s)
    for(size_t i=0; i<nod.size(); i++) g2l[nod[i]] = i;

  // Read element mapping if there are neubc files to map
  if (nbc_files.size() > 0) {
    std::cout << "Reading element map " << opts.submsh_base + EIDX_EXT << std::endl;
    binary_read(eidx, opts.submsh_base + EIDX_EXT);
    if(map_m2s)
      for(size_t i=0; i<eidx.size(); i++) g2l_ele[eidx[i]] = i;
  }

  std::ofstream info(opts.info, std::ios_base::out);
  std::string outfile;

  if (vtx_files.size() > 0) {
    std::cout << "Mapping vtx files:" << std::endl;
    for (fit=vtx_files.begin(); fit!=vtx_files.end(); ++fit) {
      // Determine output file
      outfile = opts.outdir + "/" + mt_basename(*fit);
      // Log message
      std::cout << "  " << *fit << " -> " << outfile << std::endl;
      // add file to information
      info << outfile << std::endl;

      // Read vtx
      mt_vector<mt_int> vtx;
      read_vtx(vtx, *fit);
      // Map node indices
      if(map_m2s) {
        map_glob2loc(g2l, vtx);
      }
      else {
        for(size_t i=0; i<vtx.size(); i++) vtx[i] = nod[vtx[i]];
      }
      // Write mapped nodes
      write_vtx(vtx, outfile);
    }
  }

  if (surf_files.size() > 0) {
    std::cout << "Mapping surf files:" << std::endl;
    for (fit=surf_files.begin(); fit!=surf_files.end(); ++fit) {
      // Determine output file
      outfile = opts.outdir + "/" + mt_basename(*fit);
      // Log message
      std::cout << "  " << *fit << " -> " << outfile << std::endl;
      // add file to information
      info << outfile << std::endl;

      // Read surface
      mt_meshdata mesh; // input + output
      readElements(mesh, *fit);
      // Map node indices
      mt_vector<bool> rem;
      if(map_m2s) {
        map_connectivity_glob2loc(g2l, mesh.e2n_cnt, mesh.e2n_con, rem);
      }
      else {
        for(size_t i=0; i<mesh.e2n_con.size(); i++)
          mesh.e2n_con[i] = nod[mesh.e2n_con[i]];
      }
      // Write mapped surface
      write_surf(mesh.e2n_cnt, mesh.e2n_con, outfile);
    }
  }

  if (nbc_files.size() > 0) {
    std::cout << "Mapping neubc files:" << std::endl;
    for (fit=nbc_files.begin(); fit!=nbc_files.end(); ++fit) {
      // Determine output file
      outfile = opts.outdir + "/" + mt_basename(*fit);
      // Log message
      std::cout << "  " << *fit << " -> " << outfile << std::endl;
      // add file to information
      info << outfile << std::endl;

      // Read neubc
      nbc_data nbc;
      mt_vector<mt_int> nbc_cnt, nbc_con;
      read_nbc(nbc, nbc_con, *fit);
      nbc_cnt.assign(nbc_con.size() / 3, 3);

      if(map_m2s) {
        // Map node indices
        mt_vector<bool> rem;
        map_connectivity_glob2loc(g2l, nbc_cnt, nbc_con, rem);
        bool frth_err = false, eidx_err = false;

        size_t widx=0;
        for(size_t n=0; n<rem.size(); n++)
        {
          if(!rem[n]) {
            // copy tag
            nbc.tag[widx] = nbc.tag[n];

            // try to map fourth vertex index
            auto it = g2l.find(nbc.sp_vtx[n]);
            if(it != g2l.end())
              nbc.sp_vtx[widx] = it->second;
            else {
              nbc.sp_vtx[widx] = -1;
              frth_err = true;
            }

            // try to map element index
            it = g2l_ele.find(nbc.eidx[n]);
            if(it != g2l_ele.end())
              nbc.eidx[widx] = it->second;
            else {
              nbc.eidx[widx] = -1;
              eidx_err = true;
            }

            widx++;
          }
        }
        // output possible errors
        if(frth_err)
          fprintf(stderr, "meshtool map (neubc) : 4th vertex mapping error!\n");
        if(eidx_err)
          fprintf(stderr, "meshtool map (neubc) : element index mapping error!\n");

        // resize datastructs
        nbc.sp_vtx.resize(widx);
        nbc.eidx.resize(widx);
        nbc.tag.resize(widx);
      }
      else {
        for(size_t i=0; i<nbc_con.size(); i++)
          nbc_con[i] = nod[nbc_con[i]];
        for(size_t i=0; i<nbc.sp_vtx.size(); i++)
          nbc.sp_vtx[i] = nod[nbc.sp_vtx[i]];
        for(size_t i=0; i<nbc.eidx.size(); i++)
          nbc.eidx[i] = eidx[nbc.eidx[i]];
      }

      // Write mapped surface
      write_nbc(nbc, nbc_con, nod.size(), eidx.size(), outfile);
    }
  }
  info.close();

  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
}
