/**
* @file smooth_mode.h
* @brief Smoothing of meshes.
* @author Aurel Neic
* @version
* @date 2017-02-13
*/

#include "mt_modes_base.h"


/**
* @brief The type of smoothing: Either surface or volume.
*/
enum sm_type {SM_SURF, SM_MESH, SM_DATA};

/**
* @brief Options of smoothing mode.
*/
struct smooth_options {
  mt_filename msh;
  mt_filename outmsh;

  std::string surf;
  std::string tags;
  std::string smooth;
  std::string iter;
  std::string edge;
  std::string thr;
  std::string lvl;

  std::string idat;
  std::string odat;
  std::string nodal;

  sm_type     type;
};


#define LVL_DEFAULT 2

static const std::string nodal_par = "-nodal=";


/**
* @brief Surface smoothing help message.
*/
void print_smooth_surf_help()
{
  fprintf(stderr, "smooth surface: smooth one or multiple surfaces of a mesh\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t (input) path to basename of the mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<path1>,<path2>\t (input) list of surfaces to smooth.\n", surf_par.c_str());
  fprintf(stderr, "%s<float>\t\t (optional) Maximum allowed element quality metric. default is %g.\n"
                  "\t\t\t Set to 0 to disable quality checking.\n", thr_par.c_str(), SMOOTH_THR_DEFAULT);
  fprintf(stderr, "%s<int>\t\t (optional) Number of volumetric element layers to to use\n"
                  "\t\t\t when smoothing surfaces. Default is %d.\n", lvl_par.c_str(), LVL_DEFAULT);
  fprintf(stderr, "%s<int>\t\t (optional) Number of smoothing iter (default %d).\n", iter_par.c_str(), SMOOTH_ITER_DEFAULT);
  fprintf(stderr, "%s<float>\t\t (optional) Smoothing coefficient (default %.2f).\n", smooth_par.c_str(), SMOOTH_DEFAULT);
  fprintf(stderr, "%s<float>\t\t (optional) Normal-vector angle difference (in degrees, default %g) defining a sharp edge.\n"
                  "\t\t\t If set to 0, edge detection is turned off. Negative values let meshtool skip edge vertices.\n", edge_par.c_str(), EDGE_DETECT_DEFAULT);
  fprintf(stderr, "%s<path>\t\t (output) path to basename of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "%s<format>\t\t (optional) mesh input format.\n", inp_format_par.c_str());
  fprintf(stderr, "%s<format>\t\t (optional) mesh output format.\n\n", out_format_par.c_str());
  fprintf(stderr, "The supported input formats are:\n%s\n", input_formats.c_str());
  fprintf(stderr, "The supported output formats are:\n%s\n", output_formats.c_str());
  fprintf(stderr, "\n");
}

/**
* @brief Volume smoothing help message.
*/
void print_smooth_mesh_help()
{
  fprintf(stderr, "smooth mesh: smooth surfaces and volume of a mesh\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t (input) path to basename of the mesh\n", mesh_par.c_str());
  fprintf(stderr, "%sta1,ta2/tb1,tb2\t (input) List of tag sets. The tags in one set have a common surface.\n"
                  "\t\t\t Surfaces between different tag sets will be smoothed. Use * to include all tags, one \n"
                  "\t\t\t tag per set. Use + to include all tags into one set.\n",
                  tags_par.c_str());
  fprintf(stderr, "%s<float>\t\t (optional) Maximum allowed element quality metric. default is %g.\n"
                  "\t\t\t Set to 0 to disable quality checking.\n", thr_par.c_str(), SMOOTH_THR_DEFAULT);
  fprintf(stderr, "%s<int>\t\t (optional) Number of smoothing iter (default %d).\n", iter_par.c_str(), SMOOTH_ITER_DEFAULT);
  fprintf(stderr, "%s<float>\t\t (optional) Smoothing coefficient (default %.2f).\n", smooth_par.c_str(), SMOOTH_DEFAULT);
  fprintf(stderr, "%s<float>\t\t (optional) Normal-vector angle difference (in degrees, default %g) defining a sharp edge.\n"
                  "\t\t\t If set to 0, edge detection is turned off. Negative values let meshtool skip edge vertices.\n", edge_par.c_str(), EDGE_DETECT_DEFAULT);
  fprintf(stderr, "%s<path>\t\t (output) path to basename of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "%s<format>\t\t (optional) mesh input format.\n", inp_format_par.c_str());
  fprintf(stderr, "%s<format>\t\t (optional) mesh output format.\n\n", out_format_par.c_str());
  fprintf(stderr, "The supported input formats are:\n%s\n", input_formats.c_str());
  fprintf(stderr, "The supported output formats are:\n%s\n", output_formats.c_str());
  fprintf(stderr, "\n");
}

/**
* @brief Data smoothing help message.
*/
void print_smooth_data_help()
{
  fprintf(stderr, "smooth data: smooth data defined on a mesh\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to the input data file\n", idat_par.c_str());
  fprintf(stderr, "%s<int>\t (optional) Number of smoothing iter (default %d).\n",
                  iter_par.c_str(), SMOOTH_ITER_DEFAULT);
  fprintf(stderr, "%s<float>\t (optional) Smoothing coefficient (default %.2f).\n",
                  smooth_par.c_str(), SMOOTH_DEFAULT);
  fprintf(stderr, "%s<0|1>\t (optional) Set data representation: 0 = element data, 1 = nodal data. (default is 1).\n",
                  nodal_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to the output data file\n",
                  odat_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh input format.\n\n", inp_format_par.c_str());
  fprintf(stderr, "The supported input formats are:\n%s\n", input_formats.c_str());
  fprintf(stderr, "\n");
}


/**
* @brief Smoothing mode options parser.
*
* @param [in]  argc Arguments count.
* @param [in]  argv Arguments string-array.
* @param [out] opts Options structure.
*
* @return
*/
int smooth_parse_options(int argc, char** argv, struct smooth_options & opts)
{
  if(argc < 3) {
    fprintf(stderr, "Please choose on of the following modes:\n");
    print_smooth_surf_help();
    print_smooth_mesh_help();
    print_smooth_data_help();
    return 1;
  }

  std::string smoothtype = argv[2];
  std::string msh_base;
  std::string outmsh_base;
  std::string ifmt;
  std::string ofmt;

  // parse parameters -----------------------------------------------------------------
  for(int i=3; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh_par, msh_base);
    if(!match) match = parse_param(param, outmesh_par, outmsh_base);
    if(!match) match = parse_param(param, surf_par, opts.surf);
    if(!match) match = parse_param(param, tags_par, opts.tags);
    if(!match) match = parse_param(param, smooth_par, opts.smooth);
    if(!match) match = parse_param(param, iter_par, opts.iter);
    if(!match) match = parse_param(param, edge_par, opts.edge);
    if(!match) match = parse_param(param, thr_par, opts.thr);
    if(!match) match = parse_param(param, lvl_par, opts.lvl);
    if(!match) match = parse_param(param, out_format_par, ofmt);
    if(!match) match = parse_param(param, inp_format_par, ifmt);

    if(!match) match = parse_param(param, idat_par, opts.idat);
    if(!match) match = parse_param(param, odat_par, opts.odat);
    if(!match) match = parse_param(param, nodal_par, opts.nodal);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 2;
    }
  }

  opts.msh.assign(msh_base, ifmt);
  opts.outmsh.assign(outmsh_base, ofmt);

  if(smoothtype.compare("surface") == 0)
  {
    opts.type = SM_SURF;

    // check if all relevant parameters have been set ---------------------------------------------------
    if( ! (opts.msh.isSet() && opts.outmsh.isSet() && opts.surf.size() > 0) )
    {
      std::cerr << "smooth error: Insufficient parameters provided." << std::endl;
      print_smooth_surf_help();
      return 4;
    }
  }
  else if(smoothtype.compare("mesh") == 0)
  {
    opts.type = SM_MESH;

    // check if all relevant parameters have been set ---------------------------------------------------
    if( ! (opts.msh.isSet() && opts.outmsh.isSet() && opts.tags.size() > 0) )
    {
      std::cerr << "smooth error: Insufficient parameters provided." << std::endl;
      print_smooth_mesh_help();
      return 4;
    }
  }
  else if(smoothtype.compare("data") == 0)
  {
    opts.type = SM_DATA;

    // check if all relevant parameters have been set ---------------------------------------------------
    if( ! (opts.msh.isSet() && opts.idat.size() && opts.odat.size()) )
    {
      std::cerr << "smooth error: Insufficient parameters provided." << std::endl;
      print_smooth_data_help();
      return 4;
    }
  }
  else {
    print_usage(argv[0]);
    return 2;
  }
  return 0;
}



/**
* @brief Smoothing mode function.
*
* @param [in]  argc Arguments count.
* @param [in]  argv Arguments string-array.
*/
void smooth_mode(int argc, char** argv)
{
  struct smooth_options opts;
  int ret = smooth_parse_options(argc, argv, opts);
  struct timeval t1, t2;

  if(ret != 0) return;

  struct mt_meshdata mesh;

  // parse smoothing parameters
  mt_real smooth  = opts.smooth.size() > 0 ? atof(opts.smooth.c_str()) : SMOOTH_DEFAULT;
  mt_int  iter    = opts.iter.size()   > 0 ? atoi(opts.iter.c_str())   : SMOOTH_ITER_DEFAULT;
  mt_real ang_thr = opts.edge.size()   > 0 ? atof(opts.edge.c_str())   : EDGE_DETECT_DEFAULT;
  bool skip_lines = false;
  if(ang_thr < 0) {
    skip_lines = true;
    ang_thr *= -1.0f;
  }

  // quality preservation
  mt_real thr = opts.thr.size() > 0 ? atof(opts.thr.c_str()) : SMOOTH_THR_DEFAULT;

  // read mesh
  std::cout << "Reading mesh: " << opts.msh.base << std::endl;
  gettimeofday(&t1, NULL);
  read_mesh_selected(mesh, opts.msh.format, opts.msh.base);
  std::cout << "Mesh consists of " << mesh.e2n_cnt.size() << " elements, " << mesh.xyz.size()/3 << " nodes." << std::endl;

  std::cout << "Setting up n2e / n2n graphs for mesh .. " << std::endl;
  compute_full_mesh_connectivity(mesh, opts.msh.base, true);

  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  bool write_mesh = false;

  switch(opts.type) {
    case SM_SURF:
    {
      mt_meshdata surfmesh;
      mt_meshdata linemesh;
      MT_USET<mt_int> smth_vtx_set, bad_surf_vtx, bad_ln_vtx;

      unified_surface_from_list(opts.surf, ',' , surfmesh);

      std::cout << "Setting up mesh surface data.. " << std::endl;
      gettimeofday(&t1, NULL);
      if( ! set_in_interval(surfmesh.e2n_con.begin(), surfmesh.e2n_con.end(),
                            mt_int(0), mt_int(mesh.n2e_cnt.size())) ) {
        fprintf(stderr, "%s error: surface connectivity does not match mesh connectivity! "
                        "Aborting!\n", __func__);
        exit(1);
      }
      compute_full_mesh_connectivity(surfmesh);
      remove_bad_surface_edges(mesh, surfmesh);

      // select volumetric vertices we want to smooth
      int numlvl = opts.lvl.size() > 0 ? atoi(opts.lvl.c_str()) : LVL_DEFAULT;
      smth_vtx_set.insert(surfmesh.e2n_con.begin(), surfmesh.e2n_con.end());

      // for volumetric meshes we expand the smoothing numlvl layers into the volume
      if(mesh.etype[0] != Tri && mesh.etype[0] != Quad) {
        MT_USET<mt_int> eset;
        for(int i=0; i<numlvl; i++) {
          nodeSet_to_elemSet(mesh, smth_vtx_set, eset);
          elemSet_to_nodeSet(mesh, eset, smth_vtx_set);
        }
      }
      smth_vtx_set.sort();

      // check if we are working on an open surface
      identify_surface_border_nodes(surfmesh, bad_surf_vtx);
      if(bad_surf_vtx.size()) {
        printf("%s warning: The derived surface to smooth is open! "
            "Skipping surface border smoothing.\n", __func__);
      }

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Computing line interfaces between surfaces .." << std::endl;
      gettimeofday(&t1, NULL);
      compute_line_interfaces(mesh, surfmesh, ang_thr, false, linemesh);

      std::cout << "Setting up n2e / n2n graphs for lines .. " << std::endl;
      compute_full_mesh_connectivity(linemesh);

      // line vertices connected to more than two neighbours are corners and not proper
      // line manifold vertices -> smoothing does not work for them. So we remove from the
      // smoothing vertex set later on. Note that we check n2n_cnt > 3 since one of the
      // edges in the n2n graph is the node itself.
      for(size_t i=0; i<linemesh.n2n_cnt.size(); i++)
        if(linemesh.n2n_cnt[i] > 3) bad_ln_vtx.insert(i);

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Smoothing .." << std::endl;
      gettimeofday(&t1, NULL);

      mt_vector<mt_int> surf_nodes(surfmesh.e2n_con);
      mt_vector<mt_int> line_nodes(linemesh.e2n_con);
      binary_sort(surf_nodes); unique_resize(surf_nodes);
      binary_sort(line_nodes); unique_resize(line_nodes);

      // we do not smooth over the line nodes in the volumetric pass
      for(const mt_int & n : line_nodes) smth_vtx_set.erase(n);
      // remove open surface border from smoothing
      for(auto n : bad_surf_vtx) smth_vtx_set.erase(n);

      mt_mask mnfld(mesh.xyz.size() / 3);
      mt_vector<mt_int> sm_nodes;
      quality_aware_smoother qsmooth;
      size_t viter = iter*0.66, siter = iter*0.33;

      // smooth the volume with surface as manifold -----------------------------------------
      sm_nodes.assign(smth_vtx_set.begin(), smth_vtx_set.end());
      mnfld.insert(surf_nodes.begin(), surf_nodes.end());

      if(thr) qsmooth(mesh, mesh, surfmesh, mnfld, sm_nodes, viter, smooth, thr);
      else    smooth_nodes(mesh, surfmesh, mnfld, sm_nodes, viter, smooth);

      // smooth surface with lines as manifold -----------------------------------------------
      mnfld.clear(); smth_vtx_set.clear();
      smth_vtx_set.insert(surf_nodes.begin(), surf_nodes.end());

      if(skip_lines) {
        for(const mt_int & n : line_nodes) smth_vtx_set.erase(n);
      }
      for(auto n : bad_surf_vtx) smth_vtx_set.erase(n);

      sm_nodes.assign(smth_vtx_set.begin(), smth_vtx_set.end());
      mnfld.insert(line_nodes.begin(), line_nodes.end());

      // a shallow copy of coords into surfmesh
      surfmesh.xyz.assign(mesh.xyz.size(), mesh.xyz.data(), false);

      if(thr) qsmooth(mesh, surfmesh, linemesh, mnfld, sm_nodes, siter, smooth, thr);
      else    smooth_nodes(surfmesh, linemesh, mnfld, sm_nodes, siter, smooth);

      surfmesh.xyz.assign(0, NULL, false);

      write_mesh = true;
      break;
    }

    case SM_MESH:
    {
      std::cout << "Processing tag surfaces into unified surface .." << std::endl;
      gettimeofday(&t1, NULL);

      mt_vector<MT_USET<mt_int> > selected_tags;
      if(opts.tags[0] == '*')
      {
        // we take all tags, one surface per tag
        std::cout << "All tags selected. Generating one surface per tag." << std::endl;
        MT_USET<mt_int> alltags;
        alltags.insert(mesh.etags.begin(), mesh.etags.end());
        alltags.sort();

        selected_tags.resize(alltags.size());
        auto it = alltags.begin();

        for(size_t i=0; i<selected_tags.size(); i++, ++it) {
          selected_tags[i].insert(*it);
          printf("{ %ld }\n", long(*it));
        }
      }
      else if(opts.tags[0] == '+')
      {
        std::cout << "All tags selected: " << std::endl;
        selected_tags.resize(1);
        selected_tags[0].insert(mesh.etags.begin(), mesh.etags.end());
        selected_tags[0].sort();
        printf("{ ");
        for(auto t : selected_tags[0]) printf("%ld ", long(t));
        printf("}\n");
      }
      else {
        std::cout << "Generating one surface per tag-set. Selected sets:" << std::endl;
        // take specified tags
        mt_vector<std::string> oper;
        split_string(opts.tags, '/', oper);
        selected_tags.resize(oper.size());

        for(size_t i=0; i<oper.size(); i++) {
          printf("{ ");
          mt_vector<std::string> tags_str;
          split_string(oper[i], ',', tags_str);

          for(size_t j=0; j<tags_str.size(); j++) {
            int t = atoi(tags_str[j].c_str());
            selected_tags[i].insert(t);
            printf("%d ", t);
          }
          printf("}\n");
        }
      }

      volumetric_smooth_from_tags(mesh, selected_tags,
                                  iter, smooth, ang_thr, thr, skip_lines, true);
      write_mesh = true;
      break;
    }

    case SM_DATA:
    {
      // find out what data we are dealing with =======================================
      short data_idx = -1;
      mt_vector<mt_real>            data;
      mt_vector<mt_point<mt_real> > data_vec;
      igb_header igb, igb_out;
      bool nodal_data = opts.nodal.size() == 0 || atoi(opts.nodal.c_str());

      setup_data_format(opts.idat, opts.odat, data_idx, igb, igb_out);

      // smooth data ==================================================================
      mt_vector<mt_int> sm_nod(mesh.n2n_cnt.size());
      mt_vector<bool>   do_smooth(mesh.n2n_cnt.size(), true);
      for(size_t i=0; i<sm_nod.size(); i++) sm_nod[i] = i;

      std::cout << "Smoothing data .." << std::endl;
      gettimeofday(&t1, NULL);
      switch(data_idx) {
        case 0:
        case 2:
          read_vector_ascii(data, opts.idat, true);

          if(data_idx == 0) {
            smooth_data(mesh, sm_nod, do_smooth, iter, smooth, data, nodal_data);
            write_vector_ascii(data, opts.odat, 1);
          }
          else {
            array_to_points(data, data_vec);
            smooth_data(mesh, sm_nod, do_smooth, iter, smooth, data_vec, nodal_data);
            write_vector_ascii(data_vec, opts.odat);
          }
          break;

        case 1:
        case 3:
        {
          write_igb_header(igb_out);

          std::vector<std::vector<mt_real> > rbuff;
          printf("Processing igb time-slices: %s to %s: \n", opts.idat.c_str(), opts.odat.c_str());

          for(int t=0; t<igb.v_t; t++) {
            printf("\rcurrent time-slice %d / %d .. ", t+1, int(igb.v_t));
            fflush(stdout);

            read_igb_block(rbuff, 1, igb);
            data.assign(rbuff[0].begin(), rbuff[0].end());

            if(data_idx == 1) {
              smooth_data(mesh, sm_nod, do_smooth, iter, smooth, data, nodal_data);
            }
            else {
              array_to_points(data, data_vec);
              smooth_data(mesh, sm_nod, do_smooth, iter, smooth, data_vec, nodal_data);
              points_to_array(data_vec, data);
            }

            rbuff[0].assign(data.begin(), data.end());
            write_igb_block(rbuff, igb_out);
          }
          printf("\n");

          fclose(igb.fileptr);
          fclose(igb_out.fileptr);
          break;
        }
        default: break;
      }
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      break;
    }
  }

  if(write_mesh) {
    // we correct inside-out tets since some 3rd party tools dont keep the orientation
    std::cout << "Correcting inside-out tets .." << std::endl;
    correct_insideOut(mesh);
    gettimeofday(&t2, NULL);
    std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

    std::cout << "Writing mesh: " << opts.outmsh.base << std::endl;
    gettimeofday(&t1, NULL);
    write_mesh_selected(mesh, opts.outmsh.format, opts.outmsh.base);
    gettimeofday(&t2, NULL);
    std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
  }
}

