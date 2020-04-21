#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "mt_modes_base.h"


struct retag_options{
  std::string msh_base;
  std::string outmsh_base;
  std::string ifmt;
  std::string ofmt;
  std::string con;
  std::string neigh;
};


#define CON_DFLT 3
#define NEIGH_DFLT 2

static const std::string con_par = "-ncon=";
static const std::string neigh_par = "-neigh=";

void print_retag_help()
{
  fprintf(stderr, "retag: identify and re-tag elements that are not connected to enough elements of the same tag.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh.\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to basename of the ouput mesh.\n", outmesh_par.c_str());
  fprintf(stderr, "%s<int>\t (optional) Number of common nodes that define a element connection. Default is %d.\n",
      con_par.c_str(), CON_DFLT);
  fprintf(stderr, "%s<int>\t (optional) Number of neighbours threshold. Elements with less neighbours are isolated. Default is %d.\n",
      neigh_par.c_str(), NEIGH_DFLT);
  fprintf(stderr, "%s<format>\t (optional) mesh input format. may be: %s\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh output format. may be: %s\n", out_format_par.c_str(), output_formats.c_str());
  fprintf(stderr, "\n");
}

int retag_parse_options(int argc, char** argv, struct retag_options & opts)
{
  if(argc < 2) {
    print_retag_help();
    return 1;
  }

  // parse all retag parameters -----------------------------------------------------------------
  for(int i=1; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh_par, opts.msh_base);
    if(!match) match = parse_param(param, outmesh_par, opts.outmsh_base);
    if(!match) match = parse_param(param, out_format_par, opts.ofmt);
    if(!match) match = parse_param(param, inp_format_par, opts.ifmt);
    if(!match) match = parse_param(param, con_par, opts.con);
    if(!match) match = parse_param(param, neigh_par, opts.neigh);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 3;
    }
  }
  fixBasename(opts.msh_base);
  fixBasename(opts.outmsh_base);

  // check if all relevant parameters have been set ---------------------------------------------------
  bool mshok = opts.msh_base.size() > 0, outmshok = opts.outmsh_base.size() > 0;

  if( !(mshok && outmshok) )
  {
    std::cerr << "Error: Insufficient parameters provided." << std::endl;
    print_retag_help();
    return 2;
  }

  return 0;
}

int main(int argc, char** argv)
{
  struct timeval t1, t2;
  struct retag_options opts;
  struct mt_meshdata mesh;

  int ret = retag_parse_options(argc, argv, opts);
  if (ret != 0) return 1;

  std::cout << "Reading mesh: " << opts.msh_base << std::endl;
  gettimeofday(&t1, NULL);
  read_mesh_selected(mesh, opts.ifmt, opts.msh_base);
  compute_full_mesh_connectivity(mesh);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  mt_int ncon   = opts.con.size() ? atoi(opts.con.c_str()) : CON_DFLT;
  mt_int nneigh = opts.neigh.size() ? atoi(opts.neigh.c_str()) : NEIGH_DFLT;

  // get unique tag list
  mt_vector<mt_int> tags = mesh.etags;
  binary_sort(tags); unique_resize(tags);
  mt_int update = 0;

  do {
    update = 0;
    std::set<mt_int> bad_elems;
    PROGRESS<size_t> prg(tags.size(), "Identifying isolated elements: ");

    gettimeofday(&t1, NULL);
    for(size_t i=0; i<tags.size(); i++)
    {
      mt_int ctag = tags[i];
      MT_USET<mt_int> tagset; tagset.insert(ctag);
      mt_meshgraph mg;
      extract_tagged_meshgraph(tagset, mesh, mg);
      mt_vector<mt_int> & e2n_cnt = mg.e2n_cnt;
      mt_vector<mt_int> & e2n_con = mg.e2n_con;
      mt_vector<mt_int> n2e_cnt, n2e_con, e2e_cnt, e2e_con, e2e_mlt;

      transpose_connectivity(e2n_cnt, e2n_con, n2e_cnt, n2e_con);
      multiply_connectivities(e2n_cnt, e2n_con, n2e_cnt, n2e_con,
                              e2e_cnt, e2e_con, e2e_mlt);
      restrict_connectivity(e2e_cnt, e2e_con, e2e_mlt, ncon);

      for(size_t eidx = 0; eidx < e2e_cnt.size(); eidx++)
      {
        if( (e2e_cnt[eidx] - 1) < nneigh )
          bad_elems.insert(mg.eidx[eidx]);
      }
      prg.next();
    }
    prg.finish();
    gettimeofday(&t2, NULL);
    std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

    std::cout << "Found " << bad_elems.size() << " isolated elements." << std::endl;
    if(bad_elems.size())
    {
      #if 0
      std::cout << "Atempting retagging .." << std::endl;
      gettimeofday(&t1, NULL);

      // iterate over bad elements and change their tags to the face-neighbour tags
      for(auto it = bad_elems.begin(); it != bad_elems.end(); ++it)
      {
        mt_int eidx = *it;
        std::set<mt_int> face_neigh;
        std::set<mt_int> tag_candidates;

        mt_int v0 = mesh.e2n_con[mesh.e2n_dsp[eidx]+0];
        mt_int v1 = mesh.e2n_con[mesh.e2n_dsp[eidx]+1];
        mt_int v2 = mesh.e2n_con[mesh.e2n_dsp[eidx]+2];
        mt_int v3 = mesh.e2n_con[mesh.e2n_dsp[eidx]+3];

        // add elements which share the same faces to the set face_neigh
        elements_with_face(mesh, v0, v1, v3, face_neigh);
        elements_with_face(mesh, v0, v1, v2, face_neigh);
        elements_with_face(mesh, v0, v3, v2, face_neigh);
        elements_with_face(mesh, v1, v2, v3, face_neigh);
        face_neigh.erase(eidx);

        mt_int choice = -1;
        if(face_neigh.size())
        {
          double tag_avrg = 0.0;
          for(auto cit = face_neigh.begin(); cit != face_neigh.end(); ++cit) {
            tag_avrg += mesh.etags[*cit];
            tag_candidates.insert( mesh.etags[*cit]);
          }
          tag_avrg /= face_neigh.size();

          double dist = 1e200;
          for(auto cit = tag_candidates.begin(); cit != tag_candidates.end(); ++cit)
          {
            if( dist > fabs(double(*cit) - tag_avrg) ) {
              dist = fabs(double(*cit) - tag_avrg);
              choice = *cit;
            }
          }
        }
        else {
          std::cerr << "WARNING: Element " << eidx <<
                       " shares no face to any element! Setting its tag to " << choice << "." << std::endl;
        }

        mesh.etags[eidx] = choice;
      }
      #else
      for(auto eidx : bad_elems)
      {
        mt_int* con = mesh.e2n_con.data() + mesh.e2n_dsp[eidx];
        MT_USET<mt_int> nset, eset;
        nset.insert(con, con + mesh.e2n_cnt[eidx]);
        nodeSet_to_elemSet(mesh, nset, eset);
        eset.erase(eidx);

        MT_MAP<mt_int, mt_int> tagmults;
        for(auto n : eset)
        {
          mt_int ctag = mesh.etags[n];
          mt_int mult = tagmults[ctag];
          tagmults[ctag] = ++mult;
        }

        mt_int maxmult = 0, maxtag = -1;
        for(auto mit = tagmults.begin(); mit != tagmults.end(); ++mit) {
          if(maxmult < mit->second) {
            maxmult = mit->second;
            maxtag = mit->first;
          }
        }
        if(maxtag > -1 && maxtag != mesh.etags[eidx]) {
          mesh.etags[eidx] = maxtag;
          update++;
        }
      }
      #endif
      gettimeofday(&t2, NULL);
      std::cout << "Updated " << update << " tags in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
    }
  } while(update);

  std::cout << "Writing mesh: " << opts.outmsh_base << std::endl;
  gettimeofday(&t1, NULL);
  write_mesh_selected(mesh, opts.ofmt, opts.outmsh_base);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;


  return 0;
}
