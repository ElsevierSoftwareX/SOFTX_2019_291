/**
* @file generate_mode.h
* @brief Generate a fiber file.
* @author Anton Prassl
* @version
* @date 2017-02-13
*/

#include "mt_modes_base.h"
#include "eikonal_solver.hpp"

#define BND_STEP_DFLT 0.1
#define BND_INC_DFLT 1.25

enum gen_type {GEN_FIB, GEN_MESH, GEN_DIST, GEN_BBOX, GEN_SPLIT, GEN_SURFSPLIT};
const std::string holes_par = "-holes=";

const std::string ssurf_par = "-ssurf=";
const std::string esurf_par = "-esurf=";
const std::string bdry_par  = "-bdry_layers=";
const std::string bdry_step_par      = "-bdry_step=";
const std::string bdry_step_inc_par  = "-bdry_inc=";
const std::string bdry_pres_par      = "-prsv_bdry=";

/**
* @brief Options for generate mode
*/
struct generate_options {
  gen_type type;
  std::string msh_base;
  std::string surf;
  std::string ssurf;
  std::string esurf;
  std::string outmsh_base;
  std::string ifmt;
  std::string ofmt;
  std::string holes;
  std::string ins_tag;
  std::string scale;
  std::string odat;
  std::string op;
  std::string boundary;
  std::string boundary_step;
  std::string boundary_inc;
  std::string boundary_pres;
  std::string mode;
  std::set<int> tags;
};

/**
* @brief Help for generate mode.
*/
void print_generate_fibres_help()
{
  fprintf(stderr, "generate fibres: generate default fibers for a given mesh file. the optional element"
                  " tags identify bath regions.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the input mesh\n", mesh_par.c_str());
  fprintf(stderr, "%stag1%ctag2\t (input) sequence of \",\" separated tags.\n", tags_par.c_str(), tag_separator);
  fprintf(stderr, "%s<path>\t (output) path to the basename of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "%s<int>\t (optional) number of fiber directions. Default is 2.\n", oper_par.c_str());
  fprintf(stderr, "\n");
}

/**
* @brief Help for generate mode.
*/
void print_generate_mesh_help()
{
  fprintf(stderr, "generate mesh: generate a tetrahedral mesh from a list of nested triangle surfaces.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<surf1,surf2,..>\t (input) list of paths to triangle meshes of closed surfaces.\n"
                  "\t\t\t The surface ordering must be most outer to most inner.\n", surf_par.c_str());
  fprintf(stderr, "%s<path>\t\t (output) path to the basename of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "%s<tag1,tag2,..>\t (optional) tag indices for meshed regions. Should match order of surfaces.\n", ins_tag_par.c_str());
  fprintf(stderr, "%s<x,y,z>/<x,y,z>.. (optional) Specify coordinates of hole seed points.\n", holes_par.c_str());
  fprintf(stderr, "%s<float>\t\t (optional) Scaling for average requested edge length.\n", scale_par.c_str());
  fprintf(stderr, "%s<int>\t (optional) Number of boundary element layers to include in mesh.\n"
                  "\t\t\t If negative, only the boundary mesh part will be generated.\n", bdry_par.c_str());
  fprintf(stderr, "%s<float>\t (optional) Size of a boundary layer step, relative to the average edge length. Default is %.2f.\n",
          bdry_step_par.c_str(), BND_STEP_DFLT);
  fprintf(stderr, "%s<float>\t (optional) Relative boundary layer step increment. Default is %.2f\n",
          bdry_step_inc_par.c_str(), BND_INC_DFLT);
  fprintf(stderr, "%s<int>\t\t (optional) Volumetric element sizing mode. 0 = uniform,\n"
                  "\t\t\t 1 = based on surface element size. Default is 1.\n", mode_par.c_str());
  fprintf(stderr, "%s<int>\t\t (optional) 0 = remesh surface, 1 = preserve surface. Default is 1.\n", bdry_pres_par.c_str());
  fprintf(stderr, "%s<format>\t\t (optional) mesh input format. (%s)\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "%s<format>\t\t (optional) mesh output format. (%s)\n", out_format_par.c_str(), output_formats.c_str());
  fprintf(stderr, "\n");
}

void print_generate_distancefield_help()
{
  fprintf(stderr, "generate distancefield: generate a distance field from one or between two surfaces.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t (input) path to basename of the input mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<surf1,surf2,..>\t (input) .surf surface files. Their union defines the distance field start.\n", ssurf_par.c_str());
  fprintf(stderr, "%s<surf1,surf2,..>\t (optional) .surf surface files. Their union defines the end of a relative distance field.\n", esurf_par.c_str());
  fprintf(stderr, "%stag1%ctag2\t\t (optional) sequence of tags. The associated regions well not contribute to distances.\n", tags_par.c_str(), tag_separator);
  fprintf(stderr, "%s<path>\t\t (output) path to output .dat file.\n", odat_par.c_str());
  fprintf(stderr, "%s<format>\t\t (optional) mesh input format. (%s)\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "\n");
}

void print_generate_bboxmesh_help()
{
  fprintf(stderr, "generate bboxmesh: generate the bounding box mesh of a given mesh.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the input mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to the basename of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "%s<float>\t (optional) scale applied to bbox\n", scale_par.c_str());
  fprintf(stderr, "\n");
}

void print_generate_split_help()
{
  fprintf(stderr, "generate split: generate a split file for given split operations.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the input mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) sequence of split operations.\n", oper_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) input format. (%s)\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "%s<path>\t (output) path to the split list file\n", split_par.c_str());
  fprintf(stderr, "\n");
  fprintf(stderr, "The format of the split operations is:\n");
  fprintf(stderr, "tagA1,tagA2,..:tagB1,../tagA1,..:tagB1../..\n");
  fprintf(stderr, "\",\" separates tags, \":\" separates tag groups to split, \"/\" separates split operations.\n");
  fprintf(stderr, "The splitting is applied to the elements defined by the \"tagA\" tags.\nThese MUST NOT repeat between several split operations!\n");
  fprintf(stderr, "\n");
}

void print_generate_surfsplit_help()
{
  fprintf(stderr, "generate surfsplit: generate a split file from a given surface.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t (input) path to basename of the input mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<surf1,surf2>\t (input) list of surfaces.\n", surf_par.c_str());
  fprintf(stderr, "%s<format>\t\t (optional) input format. (%s)\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "%s<path>\t\t (output) path to the split list file\n", split_par.c_str());
  fprintf(stderr, "\n");
}


/**
* @brief Generate mode argument parser.
*
* @param [in]  argc Arguments count.
* @param [in]  argv Arguments string-array.
* @param [out] opts Options structure.
*
* @return 0 for success, >0 otherwise.
*/
int generate_parse_options(int argc, char** argv, struct generate_options & opts)
{
  if(argc < 3) {
    print_generate_fibres_help();
    print_generate_mesh_help();
    print_generate_distancefield_help();
    print_generate_bboxmesh_help();
    print_generate_split_help();
    print_generate_surfsplit_help();
    return 1;
  }

  std::string gen_type_str = argv[2];

  for (int i=3; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh_par, opts.msh_base);
    if(!match) match = parse_param(param, surf_par, opts.surf);
    if(!match) match = parse_param(param, ssurf_par, opts.ssurf);
    if(!match) match = parse_param(param, esurf_par, opts.esurf);
    if(!match) match = parse_param(param, outmesh_par, opts.outmsh_base);
    if(!match) match = parse_param(param, holes_par, opts.holes);
    if(!match) match = parse_param(param, ins_tag_par, opts.ins_tag);
    if(!match) match = parse_param(param, scale_par, opts.scale);
    if(!match) match = parse_param(param, odat_par, opts.odat);
    if(!match) match = parse_param(param, oper_par, opts.op);
    if(!match) match = parse_param(param, bdry_par, opts.boundary);
    if(!match) match = parse_param(param, bdry_step_par, opts.boundary_step);
    if(!match) match = parse_param(param, bdry_step_inc_par, opts.boundary_inc);
    if(!match) match = parse_param(param, inp_format_par, opts.ifmt);
    if(!match) match = parse_param(param, out_format_par, opts.ofmt);
    if(!match) match = parse_param(param, mode_par, opts.mode);
    if(!match) match = parse_param(param, bdry_pres_par, opts.boundary_pres);
    if(!match) match = parse_param(param, split_par, opts.odat);
    fixBasename(opts.msh_base);
    fixBasename(opts.outmsh_base);

    if( (!match) && param.compare(0, tags_par.size(), tags_par) == 0) {
      std::string tags_str;
      tags_str.assign(param.begin()+tags_par.size(), param.end());
      mt_vector<std::string> taglist;
      split_string(tags_str, tag_separator, taglist);
      for (size_t t=0; t<taglist.size(); t++)
        opts.tags.insert(atoi(taglist[t].c_str()));
      match = true;
    }

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 3;
    }
  }

  if(gen_type_str.compare("fibres") == 0)
  {
    opts.type = GEN_FIB;
    if (opts.msh_base.size() == 0 || opts.outmsh_base.size() == 0)
    {
      std::cerr << "generate fibers error: Insufficient parameters provided." << std::endl;
      print_generate_fibres_help();
      return 4;
    }
  }
  else if(gen_type_str.compare("mesh") == 0) {
    opts.type = GEN_MESH;
    if (opts.surf.size() == 0 || opts.outmsh_base.size() == 0)
    {
      std::cerr << "generate mesh error: Insufficient parameters provided." << std::endl;
      print_generate_mesh_help();
      return 4;
    }
  }
  else if(gen_type_str.compare("distancefield") == 0) {
    opts.type = GEN_DIST;
    if (opts.ssurf.size() == 0 || opts.odat.size() == 0 || opts.msh_base.size() == 0)
    {
      std::cerr << "generate distancefield error: Insufficient parameters provided." << std::endl;
      print_generate_distancefield_help();
      return 4;
    }
  }
  else if(gen_type_str.compare("bboxmesh") == 0) {
    opts.type = GEN_BBOX;
    if (opts.msh_base.size() == 0 || opts.outmsh_base.size() == 0)
    {
      std::cerr << "generate bboxmesh error: Insufficient parameters provided." << std::endl;
      print_generate_bboxmesh_help();
      return 4;
    }
  }
  else if(gen_type_str.compare("split") == 0) {
    opts.type = GEN_SPLIT;
    if (opts.msh_base.size() == 0 || opts.odat.size() == 0 || opts.op.size() == 0)
    {
      std::cerr << "generate split error: Insufficient parameters provided." << std::endl;
      print_generate_split_help();
      return 4;
    }
  }
  else if(gen_type_str.compare("surfsplit") == 0) {
    opts.type = GEN_SURFSPLIT;
    if (opts.msh_base.size() == 0 || opts.odat.size() == 0 || opts.surf.size() == 0)
    {
      std::cerr << "generate surfsplit error: Insufficient parameters provided." << std::endl;
      print_generate_surfsplit_help();
      return 4;
    }
  }
  else {
    std::cerr << "Unknown mode: generate " << gen_type_str << std::endl;
    std::cerr << "Use one of:" << std::endl << std::endl;
    print_generate_fibres_help();
    print_generate_mesh_help();
    print_generate_distancefield_help();
    return 2;
  }

  return 0;
}

void parse_holes(std::string & holes_str,
                 mt_vector<float> & holes_xyz)
{
  mt_vector<std::string> holes;
  split_string(holes_str, ':', holes);

  holes_xyz.resize(holes.size()*3);

  for(size_t h=0; h<holes.size(); h++)
  {
    mt_vector<std::string> coords;
    split_string(holes[h], ',' , coords);

    if(coords.size() != size_t(3)) {
      fprintf(stderr, "%s error: Hole definition \" %s \" malformed! Aborting!\n",
              __func__, holes[h].c_str());
      exit(1);
    }

    holes_xyz[h*3+0] = atof(coords[0].c_str());
    holes_xyz[h*3+1] = atof(coords[1].c_str());
    holes_xyz[h*3+2] = atof(coords[2].c_str());
  }
}

void generate_field_from_surface(mt_meshdata & mesh, mt_meshdata & surface,
                                 mt_vector<vec3r> & field)
{
  MT_USET<mt_int> surf_nodes, ele;
  surf_nodes.insert(surface.e2n_con.begin(), surface.e2n_con.end());

  size_t num_field_nodes = mesh.xyz.size() / 3;

  // now split elements
  mt_vector<mt_real> nodal_surf_nrmls;
  compute_nodal_surface_normals(surface, mesh.xyz, nodal_surf_nrmls);

  // generate a averaged surface normal
  vec3r avrg_nrml;
  mt_real avrg_cnt = 0.0;

  for(mt_int surf_nidx : surf_nodes) {
    vec3r n = vec3r(nodal_surf_nrmls.data() + surf_nidx*3);
    if(n.length2() > 0.0) {
      avrg_nrml += n;
      avrg_cnt += 1.0;
    }
  }
  avrg_nrml /= avrg_cnt;

  field.assign(num_field_nodes, avrg_nrml);

  for(mt_int surf_nidx : surf_nodes) {
    vec3r n = vec3r(nodal_surf_nrmls.data() + surf_nidx*3);
    field[surf_nidx] = n;
  }

  #if 1
  for(int i=0; i<5; i++) {
    nodeSet_to_elemSet(mesh, surf_nodes, ele);
    elemSet_to_nodeSet(mesh, ele, surf_nodes);
  }

  // the nodes we smooth
  mt_vector<mt_int> sm_nodes; sm_nodes.assign(surf_nodes.begin(), surf_nodes.end());
  // we take contributions from all nodes
  mt_vector<bool> contribute(num_field_nodes, true);

  smooth_data(mesh, sm_nodes, contribute, 10, 0.25, field, true);
  #endif
}

void parse_scaling_factors(std::string & scale_str, size_t num_surfs, mt_vector<mt_real> & factors)
{
  if(scale_str.size() == 0) {
    factors.assign(num_surfs, 1.0);
    return;
  }

  mt_vector<std::string> scale_list;
  split_string(scale_str, ',' , scale_list);

  if(scale_list.size() == 1) {
    mt_real sz = atof(scale_list[0].c_str());
    factors.assign(num_surfs, sz);
  }
  else if(scale_list.size() == num_surfs) {
    factors.resize(num_surfs);

    size_t idx = 0;
    for(std::string & s : scale_list)
      factors[idx++] = atof(s.c_str());
  }
  else {
    fprintf(stderr, "%s error: The number of mesh size scaling factors must be either one,\n"
                    "or equal the number of input surfaces. Aborting!\n", __func__);
    exit(EXIT_FAILURE);
  }
}



/**
* @brief The generate mode function.
*
* @param [in]  argc Arguments count.
* @param [in]  argv Arguments string-array.
*/
void generate_mode(int argc, char** argv)
{
  struct timeval t1, t2;
  struct generate_options opts;

  int ret = generate_parse_options(argc, argv, opts);
  if (ret != 0) return;

  switch(opts.type) {
    case GEN_FIB:
    {
      struct mt_meshdata mesh;
      std::cout << "Reading mesh: " << opts.msh_base << std::endl;
      gettimeofday(&t1, NULL);
      readElements_general(mesh, opts.msh_base);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      int numfibers = 2;
      if(opts.op.size()) {
        int inp_fib = atoi(opts.op.c_str());
        if(inp_fib == 1 || inp_fib == 2) numfibers = inp_fib;
      }

      // creating fibers
      size_t numelem = mesh.e2n_cnt.size();
      mesh.lon.resize(numelem * numfibers * 3);

      mt_real *lon_ptr = mesh.lon.data();
      mt_int  *tag_ptr = mesh.etags.data();

      for (size_t i=0; i<numelem; i++) {
        if ( opts.tags.count(tag_ptr[i]) == 0 ) {
          // tag is not in set of bath tags
          lon_ptr[0] = 1.; lon_ptr[1] = 0.; lon_ptr[2] = 0.;
          if(numfibers == 2) {
            lon_ptr[3] = 0.; lon_ptr[4] = 1.; lon_ptr[5] = 0.;
          }
        } else {
          for(int j=0; j<numfibers*3; j++) lon_ptr[j] = 0.0;
        }

        lon_ptr += numfibers * 3;
      }

      std::cout << "Writing mesh: " << opts.outmsh_base << CARPTXT_LON_EXT << std::endl;
      gettimeofday(&t1, NULL);
      writeFibers(mesh.lon, mesh.e2n_cnt.size(), opts.outmsh_base + CARPTXT_LON_EXT);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      break;
    }

    case GEN_MESH:
    {
      mt_vector<std::string> surf_list, tags_list;
      split_string(opts.surf, ',' , surf_list);
      if(opts.ins_tag.size())
        split_string(opts.ins_tag, ',' , tags_list);

      int num_tags  = tags_list.size();
      int num_surfs = surf_list.size();
      bool tags_provided = true;

      if(num_tags != num_surfs) {
        printf("No tag list provided. Using the tag of the first element in each surface insted.\n");
        tags_provided = false;
      }

      int boundary_layer = 0;
      if(opts.boundary.size()) boundary_layer = atoi(opts.boundary.c_str());

      if(boundary_layer != 0 && num_surfs > 1) {
        fprintf(stderr, "Warning: boundary layer generation supports only one input surface!\n");
      }

      mt_meshdata surfmesh;
      mt_vector<mt_int> tags(num_surfs);
      mt_vector<kdtree*> surface_trees(num_surfs);

      gettimeofday(&t1, NULL);

      for(int sidx=0; sidx < num_surfs; sidx++)
      {
        mt_filename sfile(surf_list[sidx], opts.ifmt);
        surface_trees[sidx] = new kdtree(10);

        if(sidx == 0) {
          std::cout << "Reading surface: " << sfile.base << std::endl;
          read_mesh_selected(surfmesh, sfile.format, sfile.base);

          if(tags_provided)
            tags[sidx] = atoi(tags_list[sidx].c_str());
          else
            tags[sidx] = surfmesh.etags[0];

          // we set the tag of each surface to sidx, in this way we can identify
          // each final surface element with a source surface
          surfmesh.etags.assign(surfmesh.e2n_cnt.size(), sidx);

          correct_duplicate_vertices(surfmesh);
          correct_duplicate_elements(surfmesh);

          std::cout << "Building surface kdtree .." << std::endl;
          surface_trees[sidx]->build_tree(surfmesh);
          compute_full_mesh_connectivity(surfmesh);
        }
        else {
          mt_meshdata actsurf;
          std::cout << "Reading surface: " << sfile.base << std::endl;
          read_mesh_selected(actsurf, sfile.format, sfile.base);
          correct_duplicate_vertices(actsurf);
          correct_duplicate_elements(actsurf);

          if(tags_provided)
            tags[sidx] = atoi(tags_list[sidx].c_str());
          else
            tags[sidx] = actsurf.etags[0];

          // we set the tag of each surface to sidx, in this way we can identify
          // each final surface element with a source surface
          actsurf.etags.assign(actsurf.e2n_cnt.size(), sidx);

          compute_full_mesh_connectivity(actsurf);
          std::cout << "Unifying surfaces .." << std::endl;
          mesh_union(surfmesh, actsurf, false);
          std::cout << "Building surface kdtree .." << std::endl;
          surface_trees[sidx]->build_tree(actsurf);
        }
      }

      compute_full_mesh_connectivity(surfmesh);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      // compute average edge length and apply user-set scale
      mt_vector<mt_real> scale_factors;
      parse_scaling_factors(opts.scale, num_surfs, scale_factors);

      // the first scale factor is also the global one
      mt_real scale = scale_factors[0];
      mt_real avg_edge_length = avrg_edgelength_estimate(surfmesh) * scale;

      mt_vector<mt_real> sizes;
      mt_real boundary_step = opts.boundary_step.size() ? atof(opts.boundary_step.c_str()) :
                              BND_STEP_DFLT;
      mt_real boundary_inc  = opts.boundary_inc.size()  ? atof(opts.boundary_inc.c_str()) :
                              BND_INC_DFLT;
      int sizing_mode = 1;
      if(opts.mode.size()) sizing_mode = atoi(opts.mode.c_str());

      boundary_step *= scale;

      mt_meshdata output_mesh, boundary_mesh;
      bool only_boundary_layer = false;

      if(boundary_layer != 0) {
        if(boundary_layer < 0) {
          only_boundary_layer = true;
          boundary_layer = -boundary_layer;
        }

        mt_meshdata inner_surf;
        generate_boundary_layers(surfmesh, boundary_layer, boundary_step, boundary_inc,
                                 inner_surf, boundary_mesh);
        surfmesh.xyz = inner_surf.xyz;
      }

      if(only_boundary_layer) {
        output_mesh = boundary_mesh;
      }
      else {
        mt_vector<float> holes_xyz;
        if(opts.holes.size())
          parse_holes(opts.holes, holes_xyz);

        std::cout << "\n\nMeshing with TETGEN\n\n" << std::endl;
        gettimeofday(&t1, NULL);

        if(sizing_mode == 1) {
          generate_surfmesh_sizing_field(surfmesh, sizes);
          // for(mt_real & s : sizes) s *= scale;

          mt_vector<mt_real> final_scaling_ele(surfmesh.e2n_cnt.size()), final_scaling_nod;
          for(size_t eidx = 0; eidx < surfmesh.e2n_cnt.size(); eidx++)
            final_scaling_ele[eidx] = scale_factors[surfmesh.etags[eidx]];

          elemData_to_nodeData(surfmesh, final_scaling_ele, final_scaling_nod);

          for(size_t i=0; i<final_scaling_nod.size(); i++)
            sizes[i] *= final_scaling_nod[i];
        }

        bool preserve_boundary = true;
        if(opts.boundary_pres.size() > 0)
          preserve_boundary = atoi(opts.boundary_pres.c_str()) == 1;

        mesh_with_tetgen(surfmesh, output_mesh, holes_xyz, avg_edge_length, 1.8,
                         sizing_mode == 1 ? &sizes : nullptr, preserve_boundary);

        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

        if(boundary_layer > 0)
          mesh_union(output_mesh, boundary_mesh, true);

        // for debug purposes we sometimes dont want to overwrite the elem tags
        bool recompute_tag = true;
        if(recompute_tag) {
          std::cout << "Sampling elements " << std::endl;
          gettimeofday(&t1, NULL);

          if(num_surfs > 1) {
            mt_mask tag_found(output_mesh.e2n_cnt.size());
            // for a nested surface input, we sample the output mesh elements w.r.t the surfaces
            for(int sidx=num_surfs - 1; sidx >= 0; sidx--)
              sample_elem_tags(output_mesh, tag_found, *surface_trees[sidx], tags[sidx], 1);
          }
          else {
            // for one surf we just set the tags
            output_mesh.etags.assign(output_mesh.e2n_cnt.size(), tags[0]);
          }
          gettimeofday(&t2, NULL);
          std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
        }
      }

      std::cout << "Writing mesh " << opts.outmsh_base << " .." << std::endl;
      gettimeofday(&t1, NULL);

      mt_filename outfile(opts.outmsh_base, opts.ofmt);
      write_mesh_selected(output_mesh, outfile.format, outfile.base);

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      break;
    }

    case GEN_DIST:
    {
      mt_meshdata mesh, surf_start, surf_end;

      std::cout << "Reading mesh: " << opts.msh_base << std::endl;
      gettimeofday(&t1, NULL);
      mt_filename file(opts.msh_base, opts.ifmt);
      read_mesh_selected(mesh, file.format, file.base);
      compute_full_mesh_connectivity(mesh, file.base);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Reading start surface(s): " << opts.ssurf << std::endl;
      unified_surface_from_list(opts.ssurf, ',' , surf_start);

      bool have_end_surf = opts.esurf.size();
      if(have_end_surf) {
        std::cout << "Reading end surface(s): " << opts.esurf << std::endl;
        unified_surface_from_list(opts.esurf, ',' , surf_end);
      }

      mt_vector<mt_real> sol_start, sol_end, rel_dist;
      ek_initdata<mt_int,mt_real> ekinit;
      ekinit.vf = 1, ekinit.vs = 1, ekinit.vn = 1;

      ekinit.nod = surf_start.e2n_con;
      binary_sort(ekinit.nod); unique_resize(ekinit.nod);
      ekinit.val.assign(ekinit.nod.size(), mt_real(0.0));

      if(opts.tags.size()) {
        std::cout << "Setting fiber directions to field gradient " << std::endl;
        gettimeofday(&t1, NULL);
        eikonal_solver<mt_int,mt_real,mt_real> ekslv(mesh, ekinit);
        ekslv(sol_start);
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

        ekinit.vf = 1, ekinit.vs = 1, ekinit.vn = 1;

        mt_vector<vec3r> grad;
        mt_vector<mt_real> mag;
        compute_gradient(mesh, sol_start, true, false, grad, mag);

        mesh.lon.resize(mesh.e2n_cnt.size() * 3);
        for(size_t i=0; i<grad.size(); i++) {
          grad[i].normalize();
          mesh.lon[i*3+0] = grad[i].x;
          mesh.lon[i*3+1] = grad[i].y;
          mesh.lon[i*3+2] = grad[i].z;
        }

        for(auto t : opts.tags)
          ekinit.tag_velo[t] = {1000.0, 0.001, 0.001};
      }

      // run eikonal from start
      std::cout << "Solving eikonal from start surface " << std::endl;
      gettimeofday(&t1, NULL);
      {
        eikonal_solver<mt_int,mt_real,mt_real> ekslv(mesh, ekinit);
        ekslv(sol_start);
      }
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      mt_vector<mt_real> * rvec;

      if(have_end_surf) {
        ekinit.nod = surf_end.e2n_con;
        binary_sort(ekinit.nod); unique_resize(ekinit.nod);
        ekinit.val.assign(ekinit.nod.size(), mt_real(0.0));

        std::cout << "Solving eikonal from end surface " << std::endl;
        gettimeofday(&t1, NULL);
        {
          eikonal_solver<mt_int,mt_real,mt_real> ekslv(mesh, ekinit);
          ekslv(sol_end);
        }
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

        rel_dist.resize(sol_start.size());
        for(size_t nidx = 0; nidx < rel_dist.size(); nidx++) {
          mt_real d = sol_start[nidx] + sol_end[nidx];
          if(d > 1e-6)
            rel_dist[nidx] = sol_start[nidx] / d;
          else
            rel_dist[nidx] = 0.0;
        }

        rvec = &rel_dist;
      }
      else {
        rvec = &sol_start;
      }

      write_vector_ascii(*rvec, opts.odat, 1);
      break;
    }

    case GEN_BBOX:
    {
      mt_meshdata mesh;
      std::cout << "Reading mesh: " << opts.msh_base << std::endl;
      gettimeofday(&t1, NULL);
      mt_filename file(opts.msh_base, opts.ifmt);
      read_mesh_selected(mesh, file.format, file.base);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      mt_real scale = 1.0;
      if(opts.scale.size())
        scale = atof(opts.scale.c_str());

      mt_meshdata bbox_surf;
      generate_bbox_mesh(mesh, scale, bbox_surf);

      file.assign(opts.outmsh_base, opts.ofmt);
      std::cout << "Writing mesh: " << file.base << std::endl;
      write_mesh_selected(bbox_surf, file.format, file.base);

      break;
    }

    case GEN_SPLIT:
    {
      mt_vector<std::string> splitlist;
      mt_vector< MT_USET<mt_int> > tagsA, tagsB;

      split_string(opts.op, '/', splitlist);
      size_t numsplits = splitlist.size();

      tagsA.resize(numsplits);
      tagsB.resize(numsplits);

      // first load tag groups for all splits into tagsA and tagsB ----------------
      for(size_t i=0; i<numsplits; i++) {
        // std::cout << splitlist[i] << std::endl;
        mt_vector<std::string> tlist;
        split_string(splitlist[i], ':', tlist);
        if(tlist.size() == 2) {
          mt_vector<std::string> alist, blist;
          split_string(tlist[0], ',', alist);
          for(size_t j=0; j<alist.size(); j++) tagsA[i].insert(atoi(alist[j].c_str()));
          split_string(tlist[1], ',', blist);
          for(size_t j=0; j<blist.size(); j++) tagsB[i].insert(atoi(blist[j].c_str()));
        }
        else {
          std::cerr << "Split error: Operation " << splitlist[i] << " has bad syntax."
                    << std::endl;
          return;
        }
      }

      // we have to make sure that tags in tagA[i] do not repeat for all i
      {
        std::set<mt_int> allTagA;
        size_t numTagsA = 0;
        for(size_t i=0; i<numsplits; i++) {
          allTagA.insert(tagsA[i].begin(), tagsA[i].end());
          numTagsA += tagsA[i].size();
        }
        if(allTagA.size() < numTagsA) {
          std::cerr << "Split error: Some element regions are split several times."
                       "This is not allowed!! Check your input.." << std::endl;
          return;
        }
      }

      // read mesh ------------------------------------------------------------------
      struct mt_meshdata mesh;
      mt_filename file(opts.msh_base, opts.ifmt);

      std::cout << "Reading mesh: " << file.base << std::endl;
      gettimeofday(&t1, NULL);
      read_mesh_selected(mesh, file.format, file.base, CRP_READ_ELEM);
      transpose_connectivity(mesh.e2n_cnt, mesh.e2n_con, mesh.n2e_cnt, mesh.n2e_con);
      bucket_sort_offset(mesh.e2n_cnt, mesh.e2n_dsp);
      bucket_sort_offset(mesh.n2e_cnt, mesh.n2e_dsp);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      // index of the new indices
      mt_int Nidx = mesh.n2e_cnt.size();
      mt_vector<split_item> fin_list;

      // try to open split file
      // Now process splits ---------------------------------------------------------
      std::cout << "Computing splitlist .." << std::endl;
      gettimeofday(&t1, NULL);
      for(size_t i=0; i<numsplits; i++) {
        mt_meshgraph mgA, mgB;
        extract_tagged_meshgraph(tagsA[i], mesh, mgA);
        extract_tagged_meshgraph(tagsB[i], mesh, mgB);

        add_interface_to_splitlist(mesh, mgA, mgB, Nidx, fin_list);
      }
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      // write splitlist -------------------------------------------------------------
      std::cout << "Writing splitlist: " << opts.odat << std::endl;
      gettimeofday(&t1, NULL);
      write_split_file(fin_list, opts.odat);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      break;
    }

    case GEN_SURFSPLIT:
    {
      struct mt_meshdata mesh, surf;
      mt_filename file(opts.msh_base, opts.ifmt);

      std::cout << "Reading mesh: " << file.base << std::endl;
      gettimeofday(&t1, NULL);
      read_mesh_selected(mesh, file.format, file.base);
      compute_full_mesh_connectivity(mesh, file.base, true);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Reading surfaces .." << std::endl;
      unified_surface_from_list(opts.surf, ',', surf);
      compute_full_mesh_connectivity(surf, true);

      MT_USET<mt_int> surf_nodes;
      surf_nodes.insert(surf.e2n_con.begin(), surf.e2n_con.end());

      mt_vector<vec3r> surf_alignment_field;
      generate_field_from_surface(mesh, surf, surf_alignment_field);
      for(vec3r & v : surf_alignment_field) v *= -1.0;
      write_vector_ascii(surf_alignment_field, "surface_alignment.vec");

      mt_vector<mt_real> nod_selection(surf_alignment_field.size(), 0.0);

      std::set<mixed_tuple<mt_int, mt_real>> init;
      for(mt_int v : surf_nodes) init.insert({v, mt_real(1.0)});

      mt_vector<vec3r> vtx_nod;
      array_to_points(mesh.xyz, vtx_nod);

      mt_data_transport<mt_int, mt_real> trsp(mesh.n2n_cnt, mesh.n2n_dsp, mesh.n2n_con, vtx_nod);
      mt_real edge_len = avrg_edgelength_estimate(mesh, false);

      trsp(surf_alignment_field, mt_real(-0.05), edge_len*100.0, nod_selection, init);
      write_vector_ascii(nod_selection, "nod_selection.dat");

      MT_USET<mt_int> sel_nidx, sel_eidx;
      sel_nidx.insert(surf_nodes.begin(), surf_nodes.end());
      for(int i=0; i<3; i++) {
        nodeSet_to_elemSet(mesh, sel_nidx, sel_eidx);
        elemSet_to_nodeSet(mesh, sel_eidx, sel_nidx);
      }

      mt_vector<mt_int> idxA, idxB;

      for(mt_int eidx : sel_eidx) {
        mt_int start = mesh.e2n_dsp[eidx], stop = start + mesh.e2n_cnt[eidx];
        float selcnt = 0.0f;
        for(mt_int j=start; j<stop; j++) {
          mt_int c = mesh.e2n_con[j];
          if(nod_selection[c] > 0.0)
            selcnt += 1.0f;
        }

        if(selcnt / float(mesh.e2n_cnt[eidx]) > 0.9f)
          idxA.push_back(eidx);
        else
          idxB.push_back(eidx);
      }

      mt_meshgraph mgA, mgB;
      extract_meshgraph(idxA, mesh, mgA);
      extract_meshgraph(idxB, mesh, mgB);

      mt_vector<split_item> splitlist;
      mt_int Nidx = mesh.xyz.size() / 3;
      add_interface_to_splitlist(mesh, mgA, mgB, Nidx, splitlist);

      mt_vector<mt_int> oldidx, newidx;
      oldidx.reserve(splitlist.size()), newidx.reserve(splitlist.size());
      for(split_item & sp : splitlist) {
        oldidx.push_back(sp.oldIdx);
        newidx.push_back(sp.newIdx);
      }

      binary_sort(oldidx); unique_resize(oldidx);
      binary_sort(newidx); unique_resize(newidx);

      std::string filename = opts.odat;
      write_split_file(splitlist, filename);

      filename = opts.odat + ".new" + VTX_EXT;
      write_vtx(newidx, filename, true);
      filename = opts.odat + ".old" + VTX_EXT;
      write_vtx(oldidx, filename, true);

      break;
    }
  }
}
