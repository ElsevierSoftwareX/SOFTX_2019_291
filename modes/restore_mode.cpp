/**
* @file restore_mode.h
* @brief Restore element information.
* @author Anton Prassl
* @version
* @date 2017-02-13
*/

#include "mt_modes_base.h"

enum rst_type {RST_MAPPING};

struct restore_options {
  rst_type type;
  mt_filename msh;
  mt_filename submsh;
  std::string op;
};

void print_restore_mapping_help()
{
  fprintf(stderr, "restore mapping: restore nodal and element mapping for a submesh w.r.t. a reference mesh\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the reference mesh.\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh we need to restore mapping.\n", submesh_par.c_str());
  fprintf(stderr, "%s<int>\t (optional) restore operation type. 0 = only nodes, 1 = nodes and elem indices. Default is 1.\n", oper_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh input format.\n\n", inp_format_par.c_str());
  fprintf(stderr, "The supported input formats are:\n%s\n", input_formats.c_str());
  fprintf(stderr, "\n");
}

int restore_parse_options(int argc, char** argv, struct restore_options & opts)
{
  if (argc < 3) {
    fprintf(stderr, "Please choose on of the following modes:\n");
    print_restore_mapping_help();
    return 1;
  }

  std::string restortype = argv[2];
  std::string msh_base;
  std::string submsh_base;
  std::string ifmt;

  for (int i=3; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh_par,    msh_base);
    if(!match) match = parse_param(param, submesh_par, submsh_base);
    if(!match) match = parse_param(param, inp_format_par, ifmt);
    if(!match) match = parse_param(param, oper_par, opts.op);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 2;
    }
  }
  opts.msh.assign(msh_base, ifmt);
  opts.submsh.assign(submsh_base, ifmt);

  if(restortype.compare("mapping") == 0)
  {
    opts.type = RST_MAPPING;

    if (opts.msh.isSet() == false || opts.submsh.isSet() == false) {
      std::cerr << "restore mapping error: Insufficient parameters provided.\n" << std::endl;
      print_restore_mapping_help();
      return 3;
    }
  }
  else {
    print_usage(argv[0]);
    return 2;
  }
  return 0;
}



void restore_mode(int argc, char** argv)
{
  struct timeval t1, t2;

  struct restore_options opts;
  if (restore_parse_options(argc, argv, opts) != 0)
    return;

  switch(opts.type) {
    case RST_MAPPING:
    {
      mt_meshdata refmesh, mesh;

      std::cout << std::endl;
      std::cout << "Reading reference mesh: " << opts.msh.base << std::endl;
      gettimeofday(&t1, NULL);
      read_mesh_selected(refmesh, opts.msh.format, opts.msh.base);
      compute_full_mesh_connectivity(refmesh);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Reading submesh mesh: " << opts.submsh.base << std::endl;
      gettimeofday(&t1, NULL);
      read_mesh_selected(mesh, opts.submsh.format, opts.submsh.base);
      bucket_sort_offset(mesh.e2n_cnt, mesh.e2n_dsp);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Generating mapping.." << std::endl;
      gettimeofday(&t1, NULL);

      short oper_idx = opts.op.size() > 0 ? atoi(opts.op.c_str()) : 1;
      mt_vector<mt_int> eidx(mesh.e2n_cnt.size(), -1), nod(mesh.xyz.size() / 3, -1);

      {
        // set up vtx_to_idx mapper
        MT_MAP<triple<mt_int>,mt_int> mesh_vtx_to_idx;

        for(size_t nidx=0; nidx<nod.size(); nidx++) {
          triple<mt_int> t = {mt_int(mesh.xyz[nidx*3+0] * 1000),
                              mt_int(mesh.xyz[nidx*3+1] * 1000),
                              mt_int(mesh.xyz[nidx*3+2] * 1000)};
          mesh_vtx_to_idx[t] = nidx;
        }

        if(mesh_vtx_to_idx.size() != nod.size()) {
          fprintf(stderr, "%s error: not all vertices could be hashed uniquely! Aborting!\n",
                  __func__);
          exit(1);
        }

        // now map vertices by colocated location
        for(size_t nidx=0; nidx < refmesh.xyz.size() / 3; nidx++) {
          triple<mt_int> t = {mt_int(refmesh.xyz[nidx*3+0] * 1000),
                              mt_int(refmesh.xyz[nidx*3+1] * 1000),
                              mt_int(refmesh.xyz[nidx*3+2] * 1000)};

          auto it = mesh_vtx_to_idx.find(t);

          if(it != mesh_vtx_to_idx.end()) {
            nod[it->second] = nidx;
          }
        }

        bool all_mapped = true;
        for(auto n : nod)
          if(n == -1) {
            all_mapped = false;
            break;
          }

        if(!all_mapped) {
          fprintf(stderr, "%s error: not all vertices could be mapped! Aborting!\n",
                  __func__);
          exit(1);
        }
      }

      if(oper_idx == 1) {
        #ifdef OPENMP
        #pragma omp parallel
        #endif
        {
          mt_mask selected_nodes(refmesh.xyz.size() / 3);

          #ifdef OPENMP
          #pragma omp for schedule(guided)
          #endif
          for(size_t idx = 0; idx < eidx.size(); idx++)
          {
            // we select the nodes of the current element
            mt_int estart = mesh.e2n_dsp[idx], estop = estart + mesh.e2n_cnt[idx];
            mt_int v0 = nod[mesh.e2n_con[estart]];
            for(mt_int i=estart; i<estop; i++)
              selected_nodes.insert(nod[mesh.e2n_con[i]]);

            // iterate over all reference mesh elements connected to v0
            mt_int cstart = refmesh.n2e_dsp[v0], cstop = cstart + refmesh.n2e_cnt[v0];
            for(mt_int i=cstart; i<cstop; i++) {
              mt_int ceidx = refmesh.n2e_con[i];
              mt_int ccount = 0;
              mt_int nstart = refmesh.e2n_dsp[ceidx], nstop = nstart + refmesh.e2n_cnt[ceidx];
              // iterate over nodes of an element and check if they are selected
              for(mt_int j=nstart; j<nstop; j++)
                if(selected_nodes.count(refmesh.e2n_con[j])) ccount++;

              // if all nodes are selected, we have found our maching element and
              // set up the mapping
              if(ccount == refmesh.e2n_cnt[ceidx]) {
                eidx[idx] = ceidx;
                break;
              }
            }

            // if the element couldnt be mapped, we issue a warning
            if(eidx[idx] == -1) {
              fprintf(stderr, "%s warning: Element index %ld could not be mapped!\n",
                      __func__, (long int)idx);
            }

            // deselect the nodes. this is faster than selected_nodes.clear()
            for(mt_int i=estart; i<estop; i++)
              selected_nodes.erase(nod[mesh.e2n_con[i]]);
          }
        }
      }

      // output results under a new name
      std::string outfile = opts.submsh.base + NOD_EXT;
      std::cout << "Writing " << outfile << std::endl;
      binary_write(nod.begin(), nod.end(), outfile);
      if(oper_idx == 1) {
        outfile = opts.submsh.base + EIDX_EXT;
        std::cout << "Writing " << outfile << std::endl;
        binary_write(eidx.begin(), eidx.end(), outfile);
      }
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      break;
    }
  }
}
