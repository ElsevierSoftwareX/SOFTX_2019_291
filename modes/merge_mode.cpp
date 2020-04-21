/**
* @file merge_mode.cpp
* @brief Merging geomtries into other geometries.
* @author Aurel Neic
* @version 
* @date 2018-07-30
*/

#include "mt_modes_base.h"
#include "kdtree.h"


static const std::string sample_par  = "-sample=";
static const std::string final_par   = "-final_split=";
static const std::string msh1_par    = "-msh1=";
static const std::string msh2_par    = "-msh2=";

enum merge_type {MRG_SURF, MRG_MSH};

struct merge_options {

  merge_type type;
  mt_filename msh;
  mt_filename msh1;
  mt_filename msh2;
  mt_filename outmsh;
  mt_filename surf;
  std::string ins_tag;
  std::string prec;
  std::string sample;
  std::string final_split;
  std::string fix;
};

void print_merge_surf_help()
{
  fprintf(stderr, "merge surface: merge the geometry given by a closed surface mesh into a different mesh\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t (input) path to basename of the mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t\t (input) surface mesh to merge.\n", surf_par.c_str());
  fprintf(stderr, "%s<tag>\t\t (input) tag index of the inserted mesh elems\n", ins_tag_par.c_str());
  fprintf(stderr, "%s<float>\t\t (input) desired precision threshold length\n", thr_par.c_str());
  fprintf(stderr, "%s<int>\t\t (optional) type of element sampling. 0: whole element inside surf,\n"
                  "\t\t\t    1: element center inside surf. Default is 0.\n", sample_par.c_str());
  fprintf(stderr, "%s<int>\t (optional) Whether to do a final splitting step. 0: Don't, 1: Do it.\n"
                  "\t\t\t    Doing it will yield a smoother surface at the cost of the element quality. Default is 0.\n", final_par.c_str());
  fprintf(stderr, "%s<int>\t\t (optional) Fix boundary of mesh. 0 = no, 1 = yes. Default is 0.\n",
                  fix_par.c_str());
  fprintf(stderr, "%s<format>\t\t (optional) mesh input format. (%s)\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "%s<format>\t\t (optional) mesh output format. (%s)\n", out_format_par.c_str(), output_formats.c_str());
  fprintf(stderr, "%s<path>\t\t (output) path to basename of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "\n");
}
void print_merge_meshes_help()
{
  fprintf(stderr, "merge meshes: merge two meshes, unifying co-located vertices\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t (input) path to basename of the mesh 1\n", msh1_par.c_str());
  fprintf(stderr, "%s<path>\t\t (input) path to basename of the mesh 2\n", msh2_par.c_str());
  fprintf(stderr, "%s<format>\t\t (optional) mesh input format. (%s)\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "%s<format>\t\t (optional) mesh output format. (%s)\n", out_format_par.c_str(), output_formats.c_str());
  fprintf(stderr, "%s<path>\t\t (output) path to basename of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "\n");
}


int merge_parse_options(int argc, char** argv, struct merge_options & opts)
{
  if(argc < 3) {
    fprintf(stderr, "Please choose on of the following modes:\n");
    print_merge_surf_help();
    print_merge_meshes_help();
    return 1;
  }

  std::string mergetype = argv[2];
  std::string msh_base;
  std::string msh1_base;
  std::string msh2_base;
  std::string outmsh_base;
  std::string surf_base;
  std::string ifmt;
  std::string ofmt;

  // parse parameters -----------------------------------------------------------------
  for(int i=3; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh_par, msh_base);
    if(!match) match = parse_param(param, msh1_par, msh1_base);
    if(!match) match = parse_param(param, msh2_par, msh2_base);
    if(!match) match = parse_param(param, outmesh_par, outmsh_base);
    if(!match) match = parse_param(param, surf_par, surf_base);
    if(!match) match = parse_param(param, out_format_par, ofmt);
    if(!match) match = parse_param(param, inp_format_par, ifmt);
    if(!match) match = parse_param(param, ins_tag_par, opts.ins_tag);
    if(!match) match = parse_param(param, thr_par, opts.prec);
    if(!match) match = parse_param(param, fix_par, opts.fix);
    if(!match) match = parse_param(param, sample_par, opts.sample);
    if(!match) match = parse_param(param, final_par, opts.final_split);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 2;
    }
  }

  opts.msh.assign(msh_base, ifmt);
  opts.msh1.assign(msh1_base, ifmt);
  opts.msh2.assign(msh2_base, ifmt);
  opts.surf.assign(surf_base, "");
  opts.outmsh.assign(outmsh_base, ofmt);

  if(mergetype.compare("surface") == 0)
  {
    opts.type = MRG_SURF;

    // check if all relevant parameters have been set ---------------------------------------------------
    if( !(opts.msh.isSet() && opts.outmsh.isSet() && opts.surf.isSet() &&
          opts.ins_tag.size() && opts.prec.size()) )
    {
      std::cerr << "merge surface error: Insufficient parameters provided." << std::endl;
      print_merge_surf_help();
      return 4;
    }
  }
  else if(mergetype.compare("meshes") == 0)
  {
    opts.type = MRG_MSH;

    // check if all relevant parameters have been set ---------------------------------------------------
    if( !(opts.msh1.isSet() && opts.msh2.isSet() && opts.outmsh.isSet()) )
    {
      std::cerr << "merge meshes error: Insufficient parameters provided." << std::endl;
      print_merge_meshes_help();
      return 4;
    }
  }
  else {
    print_usage(argv[0]);
    return 2;
  }
  return 0;
}



void surface_intersect_splitting(mt_meshdata & mesh,
                                 mt_edge_splitter & splitter,
                                 mt_edge_manager & edge_manager,
                                 const kdtree & tree,
                                 const MT_USET<mt_int> & fixed_nod,
                                 const float prec,
                                 const bool final_split)
{
  timeval t1, t2;
  size_t iter = 0, max_iter = 20, nsel_merge = 0;
  char msg[1024];

  float max_len = prec * prec;
  bool having_fixed_nodes = fixed_nod.size();

  std::set<edgeele, edge_len_dsc> sel_edges, acc_split_edges;

  mt_vector<edgeele> elem_edge_lengths;
  edge_len_asc ascending;

  const MT_MAP<tuple<mt_int>,mt_int> & edges = edge_manager.edges();
  const mt_mapping<mt_int> & elem2edges = edge_manager.elem2edges();

  do {
    sel_edges.clear();

    MT_USET<mt_int> sel_elem, sel_nod;
    // mark elements connected to edges that intersect the surface for splitting
    for(auto eit = edges.begin(); eit != edges.end(); ++eit)
    {
      mt_int v1 = eit->first.v1, v2 = eit->first.v2;
      mt_int edge_idx = eit->second;

      vec3f p1(mesh.xyz.data() + v1*3), p2(mesh.xyz.data() + v2*3);
      vec3f e = p2 - p1;

      tri_elem hit_tri;
      vec3r    hit_pos;
      int did_intersec = tree.closest_intersect(ray(p1, e), 0.05, 0.95, hit_tri, hit_pos);
      mt_real edge_len2 = e.length2();

      if(did_intersec && edge_len2 > max_len) {
        mt_int start = elem2edges.bwd_dsp[edge_idx], stop = start + elem2edges.bwd_cnt[edge_idx];
        for(mt_int j=start; j<stop; j++) {
          mt_int eidx = elem2edges.bwd_con[j];
          sel_elem.insert(eidx);
        }
      }
    }

    // extend elements to split by one layer in order to reduce element anisotropy
    elemSet_to_nodeSet(mesh, sel_elem, sel_nod);
    nodeSet_to_elemSet(mesh, sel_nod, sel_elem);

    // split selected elements along their longest edge
    for(const mt_int eidx : sel_elem)
    {
      add_edges_lengths(mesh, eidx, elem_edge_lengths);

      edgeele longest = *std::max_element(elem_edge_lengths.begin(), elem_edge_lengths.end(), ascending);

      bool edge_is_fixed = false;
      if(having_fixed_nodes) {
        bool v1_fixed = fixed_nod.count(longest.v1),
             v2_fixed = fixed_nod.count(longest.v2);
        // we cannot touch edges going to the fixed surface. first, because sometimes vertices
        // switch and we might brake the logic, second for elem quality reasons
        edge_is_fixed = v1_fixed || v2_fixed;
      }

      if(!edge_is_fixed) {
        longest.split_at = 0.5f;

        auto it = edges.find({longest.v1, longest.v2});
        longest.idx = it->second;

        sel_edges.insert(longest);
      }
    }

    // do splitting
    nsel_merge = sel_edges.size();
    iter++;
    if(nsel_merge) {
      sprintf(msg, "Iter %lu: splitting %lu edges, %.4f %% of total %lu ",
              iter, nsel_merge, float(nsel_merge)/edges.size()*100.0f, edges.size());
      PROGRESS<size_t> prog(nsel_merge, msg);
      gettimeofday(&t1, NULL);

      splitter(sel_edges, prog);
      transpose_connectivity(mesh.e2n_cnt, mesh.e2n_con, mesh.n2e_cnt, mesh.n2e_con);
      mesh.n2e_dsp.resize(mesh.n2e_cnt.size());
      bucket_sort_offset(mesh.n2e_cnt, mesh.n2e_dsp);

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
    }
  }
  while(nsel_merge > 0 && iter < max_iter);

  // We maybe want to have a last pass were we do something smart to improve the
  // surface smoothness of the sampled mesh.
  if(final_split) {
    sel_edges.clear();
    for(auto eit = edges.begin(); eit != edges.end(); ++eit)
    {
      mt_int v1 = eit->first.v1, v2 = eit->first.v2;
      mt_int edge_idx = eit->second;

      vec3r p1(mesh.xyz.data() + v1*3), p2(mesh.xyz.data() + v2*3);
      vec3r e = p2 - p1;

      tri_elem hit_tri;
      vec3r    hit_pos;
      int did_intersec = tree.closest_intersect(ray(p1, e), 0.0, 1.0, hit_tri, hit_pos);

      if(did_intersec) {
        mt_real hit_dist = (hit_pos - p1).length();
        mt_real split_at = hit_dist / e.length();

        if(split_at < 0.05) // for intersection close to p1 we move p1
        {
          p1 = hit_pos;
          p1.set(mesh.xyz.data() + v1*3);
        }
        else if(split_at > 0.95)  // for intersection close to p2 we move p2
        {
          p2 = hit_pos;
          p2.set(mesh.xyz.data() + v2*3);
        }
        else if(split_at > 0.2 || split_at < 0.8) // for intersection towards the center we split
          sel_edges.insert({v1, v2, edge_idx, e.length2(), float(split_at)});
      }
    }

    nsel_merge = sel_edges.size();
    if(nsel_merge) {
      sprintf(msg, "Final merge split: splitting %lu edges, %.4f %% of total %lu ",
              nsel_merge, float(nsel_merge) / edges.size(), edges.size());
      PROGRESS<size_t> prog(nsel_merge, msg);
      gettimeofday(&t1, NULL);

      splitter(sel_edges, prog);

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
    }
  }

  // finish splitting
  correct_insideOut(mesh);
  edge_manager.refresh_edge_data(true);
}



void merge_mode(int argc, char** argv)
{
  merge_options opts;
  int ret = merge_parse_options(argc, argv, opts);
  if(ret != 0) return;

  struct timeval t1, t2;

  switch(opts.type) {
    case MRG_SURF:
    {
      mt_meshdata mesh, surface;

      // read mesh and surface
      gettimeofday(&t1, NULL);
      std::cout << "Reading mesh: " << opts.msh.base << std::endl;
      read_mesh_selected(mesh, opts.msh.format, opts.msh.base);
      std::cout << "Mesh consists of " << mesh.e2n_cnt.size() << " elements, "
                << mesh.xyz.size()/3 << " nodes." << std::endl;
      std::cout << "Reading mesh: " << opts.surf.base << std::endl;
      read_mesh_selected(surface, opts.surf.format, opts.surf.base, CRP_READ_PTS | CRP_READ_ELEM);
      std::cout << "Mesh consists of " << surface.e2n_cnt.size() << " elements, "
                << surface.xyz.size()/3 << " nodes." << std::endl;

      std::cout << "Setting up n2e / n2n graphs for mesh and surface .. " << std::endl;
      compute_full_mesh_connectivity(mesh, opts.msh.base);
      compute_full_mesh_connectivity(surface);

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      kdtree tree(10);
      tree.build_tree(surface);

      mt_edge_manager  manager(mesh);
      mt_edge_splitter splitter(manager);

      gettimeofday(&t1, NULL);
      std::cout << "Splitting edges to fit surface geometry .. " << std::endl;
      float prec = atof(opts.prec.c_str());
      bool final_split = opts.final_split.size() ? atoi(opts.final_split.c_str()) : false;
      bool fix_bnd     = opts.fix.size() ? atoi(opts.fix.c_str()) : false;

      MT_USET<mt_int> fixed_nod;
      if(fix_bnd) {
        mt_meshdata surfmesh;
        compute_surface(mesh.etype, mesh.e2n_cnt, mesh.e2n_con, surfmesh, false);

        // if fix_bnd is set, we insert the surface nodes of the geometric surface into
        // the fixed_nod set, which is tested in the splitting phase
        for(const mt_int & nidx : surfmesh.e2n_con) fixed_nod.insert(nidx);
      }

      surface_intersect_splitting(mesh, splitter, manager, tree, fixed_nod, prec, final_split);
      compute_full_mesh_connectivity(mesh);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      gettimeofday(&t1, NULL);
      std::cout << "Sampling elements of merged mesh .. " << std::endl;
      mt_int newtag = atoi(opts.ins_tag.c_str());
      int sampling_type = opts.sample.size() ? atoi(opts.sample.c_str()) : 0;
      if(final_split && sampling_type == 0) {
        fprintf(stderr, "Warning: Final split was turned on, but sampling type is %d!\n"
                        "This yields bad surfaces, setting sampling type to 1.\n\n",
                sampling_type);
        sampling_type = 1;
      }

      mt_mask tag_found(mesh.e2n_cnt.size());
      sample_elem_tags(mesh, tag_found, tree, newtag, sampling_type);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Writing mesh " << opts.outmsh.base << std::endl;
      write_mesh_selected(mesh, opts.outmsh.format, opts.outmsh.base);
      break;
    }

    case MRG_MSH:
    {
      mt_meshdata mesh1, mesh2, outmesh;
      gettimeofday(&t1, NULL);

      std::cout << "Reading mesh: " << opts.msh1.base << std::endl;
      read_mesh_selected(mesh1, opts.msh1.format, opts.msh1.base);
      compute_full_mesh_connectivity(mesh1, opts.msh1.base);
      std::cout << "Mesh consists of " << mesh1.e2n_cnt.size() << " elements, "
                << mesh1.xyz.size()/3 << " nodes." << std::endl;
      std::cout << "Reading mesh: " << opts.msh2.base << std::endl;
      read_mesh_selected(mesh2, opts.msh2.format, opts.msh2.base);
      compute_full_mesh_connectivity(mesh2, opts.msh2.base);
      std::cout << "Mesh consists of " << mesh2.e2n_cnt.size() << " elements, "
                << mesh2.xyz.size()/3 << " nodes." << std::endl;

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      gettimeofday(&t1, NULL);
      std::cout << "Merging meshes .. " << std::endl;
      mesh_union(mesh1, mesh2, false);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Writing mesh " << opts.outmsh.base << std::endl;
      write_mesh_selected(mesh1, opts.outmsh.format, opts.outmsh.base);
      break;
    }
  }
}
