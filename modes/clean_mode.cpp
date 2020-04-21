/**
* @file clean_mode.cpp
* @brief Meshtool clean mode.
* @author Aurel Neic
* @version
* @date 2017-02-13
*/

#include "mt_modes_base.h"

// #define WRITE_BEFORE_AFTER
#define SIGMA_DEFAULT 0.05

enum clean_type {CLN_QUAL, CLN_TOP, CLN_NRML};


struct clean_options {
  mt_filename msh;
  mt_filename outmsh;
  std::string surf;
  std::string thr;
  std::string smooth;
  std::string iter;
  std::string edge;
  std::string sigma;
  clean_type  type;
};

static const std::string sigma_par = "-sigma=";

void print_clean_quality_help()
{
  fprintf(stderr, "clean quality: deform mesh elements to reach a certain quality threshold value. Provided surfaces will be preserved.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t (input) path to basename of the mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<float>\t\t (input) threshold for mesh quality\n", thr_par.c_str());
  fprintf(stderr, "%s<path1>,<path2>\t (optional) list of surfaces which must be preserved.\n", surf_par.c_str());
  fprintf(stderr, "%s<float>\t\t (optional) Smoothing coefficient (default %.2f). If set to 0, mesh smoothing is turned off.\n", smooth_par.c_str(), double(0));
  fprintf(stderr, "%s<int>\t\t (optional) Number of smoothing iterations (default %d).\n", iter_par.c_str(), SMOOTH_ITER_DEFAULT);
  fprintf(stderr, "%s<float>\t\t (optional) Edge detection threshold angle in degrees (default %.2f). If set to 0, edge detection is turned off.\n", edge_par.c_str(), EDGE_DETECT_DEFAULT);
  fprintf(stderr, "%s<float>\t\t (optional) Relative vertex translation step size (default %.2f).\n", sigma_par.c_str(), SIGMA_DEFAULT);
  fprintf(stderr, "%s<format>\t\t (optional) mesh output format. may be: %s\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "%s<format>\t\t (optional) mesh output format. may be: %s\n", out_format_par.c_str(), output_formats.c_str());
  fprintf(stderr, "%s<path>\t\t (output) path to basename of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "\n");
}

void print_clean_topology_help()
{
  fprintf(stderr, "clean topology: clean the mesh from bad topology definitions.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to basename of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh output format. may be: %s\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh output format. may be: %s\n", out_format_par.c_str(), output_formats.c_str());
  fprintf(stderr, "\n");
}

void print_clean_normals_help()
{
  fprintf(stderr, "clean normals: make normal directions of a surface consistent.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to basename of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "%s<surface>\t (optional) surface to compute the normals for\n", surf_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh input format. may be: %s\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "\n");
}

int clean_parse_options(int argc, char** argv, struct clean_options & opts)
{
  if(argc < 3) {
    print_clean_quality_help();
    print_clean_topology_help();
    // print_clean_normals_help();
    return 1;
  }

  std::string cleantype = argv[2];

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
    if(!match) match = parse_param(param, thr_par, opts.thr);
    if(!match) match = parse_param(param, smooth_par, opts.smooth);
    if(!match) match = parse_param(param, iter_par, opts.iter);
    if(!match) match = parse_param(param, edge_par, opts.edge);
    if(!match) match = parse_param(param, sigma_par, opts.sigma);
    if(!match) match = parse_param(param, inp_format_par, ifmt);
    if(!match) match = parse_param(param, out_format_par, ofmt);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 2;
    }
  }

  opts.msh.assign(msh_base, ifmt);
  opts.outmsh.assign(outmsh_base, ofmt);

  if(cleantype.compare("quality") == 0)
  {
    opts.type = CLN_QUAL;

    if( ! (opts.msh.isSet() && opts.outmsh.isSet() && opts.thr.size() > 0) )
    {
      std::cerr << "clean quality error: Insufficient parameters provided." << std::endl;
      print_clean_quality_help();
      return 4;
    }
  }
  else if (cleantype.compare("topology") == 0)
  {
    opts.type = CLN_TOP;

    if( ! (opts.msh.isSet() && opts.outmsh.isSet()) )
    {
      std::cerr << "clean topology error: Insufficient parameters provided." << std::endl;
      print_clean_topology_help();
      return 4;
    }
  }
  else if (cleantype.compare("normals") == 0)
  {
    opts.type = CLN_NRML;

    if( ! (opts.msh.isSet() && opts.outmsh.isSet()) )
    {
      std::cerr << "clean normals error: Insufficient parameters provided." << std::endl;
      print_clean_normals_help();
      return 4;
    }
  }
  else {
    print_usage(argv[0]);
    return 4;
  }

  return 0;
}

template<class T, class S>
void insert_elem_nbh(const mt_meshdata & mesh,
                     const T eidx,
                     std::set<T> & nbh)
{
  T estart = mesh.e2n_dsp[eidx], estop = estart + mesh.e2n_cnt[eidx];
  for(T i = estart; i<estop; i++)
  {
    T nidx = mesh.e2n_con[i];
    T nstart = mesh.n2e_dsp[nidx], nstop = nstart + mesh.n2e_cnt[nidx];
    for(T j=nstart; j<nstop; j++)
      nbh.insert(mesh.n2e_con[j]);
  }
}

void print_quality_stats(mt_vector<mt_real> elemqual)
{
  float min = elemqual[0], max = elemqual[0], sum = 0.;

  for(size_t i=0; i<elemqual.size(); i++) {
    mt_real e = elemqual[i];
    if(min > e) min = e;
    if(max < e) max = e;
    sum += e;
  }

  std::cout << "Mesh quality statistics:" << std::endl;
  std::cout << "min: " << min << ", max: " << max << ", avg: " << sum / elemqual.size() << std::endl;
}

template<class T>
void find_enclosed_elems(const mt_meshdata & mesh,
                         const std::set<T> & surf_nodes,
                         std::set<T> & encl_elems)
{
  size_t nnodes = mesh.xyz.size() / 3;
  mt_vector<bool> isSurf(nnodes, false);

  for(auto it=surf_nodes.begin(); it != surf_nodes.end(); ++it)
    isSurf[*it] = true;

  for(size_t eidx = 0, idx=0; eidx < mesh.e2n_cnt.size(); eidx++)
  {
    short c = 0;
    // increment c for every node of eidx which is on a surface
    for(T i=0; i<mesh.e2n_cnt[eidx]; i++, idx++) {
      if(isSurf[mesh.e2n_con[idx]]) c++;
    }

    if ( c == mesh.e2n_cnt[eidx] )
    {
      // element is full enclosed in a surface -> needs fixing
      encl_elems.insert(eidx);
    }
  }
}

template<class T>
void print_edge_problem(const mt_meshdata & mesh, T v0, T v1)
{
  std::set<T> first_level_elems; // elements sharing problematic edge
  std::set<T> second_level_elems; // elements connected to first-level-elems via a face

  elements_with_edge(mesh, v0, v1, first_level_elems);

  for(typename std::set<T>::iterator it = first_level_elems.begin(); it != first_level_elems.end(); ++it)
  {
    T eidx = *it;
    T ev0 = mesh.e2n_con[eidx*4+0];
    T ev1 = mesh.e2n_con[eidx*4+1];
    T ev2 = mesh.e2n_con[eidx*4+2];
    T ev3 = mesh.e2n_con[eidx*4+3];

    elements_with_face(mesh, ev0, ev1, ev3, second_level_elems);
    elements_with_face(mesh, ev0, ev1, ev2, second_level_elems);
    elements_with_face(mesh, ev0, ev3, ev2, second_level_elems);
    elements_with_face(mesh, ev1, ev2, ev3, second_level_elems);
  }
  for(typename std::set<T>::iterator it = first_level_elems.begin(); it != first_level_elems.end(); ++it)
  {
    typename std::set<T>::iterator f = second_level_elems.find(*it);
    if(f != second_level_elems.end())
      second_level_elems.erase(f);
  }

  for(typename std::set<T>::iterator it = first_level_elems.begin(); it != first_level_elems.end(); ++it)
  {
    T estart = mesh.e2n_dsp[*it];
    const T* con = mesh.e2n_con.data()+estart;
    printf("Tt %d %d %d %d 1\n", con[0], con[1], con[2], con[3]);
  }
  for(typename std::set<T>::iterator it = second_level_elems.begin(); it != second_level_elems.end(); ++it)
  {
    T estart = mesh.e2n_dsp[*it];
    const T* con = mesh.e2n_con.data()+estart;
    printf("Tt %d %d %d %d 2\n", con[0], con[1], con[2], con[3]);
  }

}


template<class T>
void inner_surf_without_outer_surf(const mt_meshdata & mesh,
                                   const MT_MAP<triple<T>, tri_sele> & outer_tri_surf,
                                   const MT_MAP<quadruple<T>, quad_sele> & outer_quad_surf,
                                   mt_meshdata & inner_surf)
{
  mt_vector<T> selected_tags(mesh.etags);
  binary_sort(selected_tags); unique_resize(selected_tags);
  {
    MT_MAP<triple<T>, tri_sele>     full_trimap;
    MT_MAP<quadruple<T>, quad_sele> full_quadmap;

    for(size_t i=0; i<selected_tags.size(); i++)
    {
      MT_USET<T> curtag;
      mt_meshgraph mg;
      MT_MAP<triple<T>, tri_sele>     cur_trimap;
      MT_MAP<quadruple<T>, quad_sele> cur_quadmap;

      curtag.insert(selected_tags[i]);
      extract_tagged_meshgraph(curtag, mesh, mg);
      compute_surface(mg.etype, mg.e2n_cnt, mg.e2n_con, cur_trimap, cur_quadmap);
      surface_union(full_trimap, full_quadmap, cur_trimap, cur_quadmap);
    }
    surface_difference(full_trimap, full_quadmap, outer_tri_surf, outer_quad_surf);
    mt_vector<T> elem_orig;
    surfmap_to_vector(full_trimap, full_quadmap, inner_surf.e2n_con, elem_orig);
    inner_surf.e2n_cnt.assign(inner_surf.e2n_con.size() / 3, 3);
  }
}


void clean_mode(int argc, char** argv)
{
  struct clean_options opts;
  int ret = clean_parse_options(argc, argv, opts);
  struct timeval t1, t2;

  if(ret != 0) return;

  struct mt_meshdata mesh;
  bool write_mesh = false; // we will set this to true if we improved something

  // read mesh
  std::cout << "Reading mesh: " << opts.msh.base << std::endl;
  gettimeofday(&t1, NULL);
  read_mesh_selected(mesh, opts.msh.format, opts.msh.base);
  compute_full_mesh_connectivity(mesh, opts.msh.base);

  std::cout << "Mesh consists of " << mesh.e2n_cnt.size() << " elements, " << mesh.xyz.size()/3 << " nodes." << std::endl;
  // we correct inside-out tets since some 3rd party tools dont keep the orientation
  std::cout << "Correcting inside-out tets .." << std::endl;
  correct_insideOut(mesh);

  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  // TODO: Distinguish bad elements as either collapsed edge or collapsed angles.
  // Write them to hd for inspection. Treat bad elements with collapsing or splitting.

  switch(opts.type)
  {
    default: break;

    case CLN_QUAL:
    {
      // threshold option required
      mt_real thr = atof(opts.thr.c_str());

      // rest is optional
      const float smooth = opts.smooth.size() > 0 ? atof(opts.smooth.c_str()) : 0;
      const int   iter   = opts.iter.size()   > 0 ? atoi(opts.iter.c_str())   : SMOOTH_ITER_DEFAULT;
      const float edge   = opts.edge.size()   > 0 ? atof(opts.edge.c_str())   : EDGE_DETECT_DEFAULT;
      const float sigma  = opts.sigma.size()  > 0 ? atof(opts.sigma.c_str())  : SIGMA_DEFAULT;
      const mt_int num_levels = 5;

      MT_USET<mt_int> surf_nodes, line_nodes;
      mt_meshdata surfmesh, linemesh;
      // retrieve the set of surface nodes from the given list of surface files
      if(opts.surf.size() > 0)
      {
        unified_surface_from_list(opts.surf, ',' , surfmesh);
        surf_nodes.insert(surfmesh.e2n_con.begin(), surfmesh.e2n_con.end());

        std::cout << "Setting up n2e / n2n graphs for mesh surface .. " << std::endl;
        compute_full_mesh_connectivity(surfmesh);

        std::cout << "Setting up line interfaces for mesh surface .. " << std::endl;
        compute_line_interfaces(mesh, surfmesh, mt_real(edge), true, linemesh);
        compute_full_mesh_connectivity(linemesh);
        line_nodes.insert(linemesh.e2n_con.begin(), linemesh.e2n_con.end());

        surfmesh.xyz.assign(mesh.xyz.size(), mesh.xyz.data(), false);
      }
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;


      std::cout << std::endl << "Checking mesh quality .. " << std::endl;
      gettimeofday(&t1, NULL);
      // compute element qualities
      MT_USET<mt_int> badelems;
      mt_vector<mt_real> elemqual;

      mesh_quality(mesh, elemqual);
      get_thres_subset(elemqual, thr, badelems);
      size_t numbadelem = badelems.size();

      print_quality_stats(elemqual);
      std::cout << "Number of elements above threshold: " << numbadelem << std::endl;
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      if(numbadelem) {
        std::cout << std::endl << "Improving mesh .. " << std::endl;
        gettimeofday(&t1, NULL);

        mt_mask isSurf(mesh.n2e_cnt.size()), isLine(mesh.n2e_cnt.size());
        isSurf.insert(surf_nodes.begin(), surf_nodes.end());
        isLine.insert(line_nodes.begin(), line_nodes.end());

        #ifdef WRITE_BEFORE_AFTER
        // write a submesh consising of the bad elements and their neighbours
        std::cout << "Writing " << opts.outmsh.base + "before.vtk" << std::endl;
        mt_meshdata badmesh = mesh;
        mt_vector<mt_int> nod, eidx;
        {
          mt_vector<bool> keep(mesh.e2n_cnt.size(), false);
          // we change the element tags to distinguish bad elements from their neigbours
          badmesh.etags.assign(badmesh.etags.size(), 0);
          std::set<mt_int> writeelem;
          std::set<mt_int>::iterator wit = badelems.begin();
          while(wit != badelems.end()) {
            insert_elem_nbh(mesh, *wit, writeelem);
            badmesh.etags[*wit] = 100;
            ++wit;
          }
          wit = writeelem.begin();
          while(wit != writeelem.end()) {
            keep[*wit] = true;
            ++wit;
          }
          restrict_meshdata(keep, badmesh, nod, eidx);
        }
        writeVTKmesh(badmesh, opts.outmsh.base + "before.vtk");
        #endif

        mt_vector<mt_int> smth_vtx_vol;
        mt_vector<mt_int> smth_vtx_surf;
        if(smooth)
        {
          MT_USET<mt_int> smoothelems(badelems);
          MT_USET<mt_int> smoothnodes;

          // we add a certain number of layers of elements to the area we will smooth
          for(mt_int i=0; i<num_levels; i++) {
            elemSet_to_nodeSet(mesh, smoothelems, smoothnodes);
            nodeSet_to_elemSet(mesh, smoothnodes, smoothelems);
          }

          // define volumetric smoothing nodes
          size_t idx=0;
          smth_vtx_vol.resize(smoothnodes.size());
          for(auto it=smoothnodes.begin(); it != smoothnodes.end(); ++it)
            if(line_nodes.count(*it) == 0)
              smth_vtx_vol[idx++] = *it;
          smth_vtx_vol.resize(idx);

          // define surface smoothing nodes
          idx=0;
          smth_vtx_surf.resize(smoothnodes.size());
          for(auto it=smoothnodes.begin(); it != smoothnodes.end(); ++it)
            if(surf_nodes.count(*it) != 0)
              smth_vtx_surf[idx++] = *it;
          smth_vtx_surf.resize(idx);
        }


        // iterate over bad elements and try to improve them
        size_t curit = 0;
        size_t numbadelem_after = badelems.size();

        mt_vector<mt_real> old_xyz;
        if(smooth) {
          size_t viter = iter*0.66, siter = iter*0.33;
          std::cout << "Smoothing .." << std::endl;
          smooth_nodes(mesh, surfmesh, isSurf, smth_vtx_vol, viter, smooth);
          smooth_nodes(surfmesh, linemesh, isLine, smth_vtx_surf, siter, smooth);
        }

        std::cout << "Shifting vertices .." << std::endl;
        while( numbadelem_after > 0 )
        {
          std::cout << "Iter " << curit++ << ": Number of elements left above threshold: "
                    << numbadelem_after << std::endl;

          numbadelem = badelems.size();
          old_xyz = mesh.xyz;

          // improve mesh
          MT_USET<mt_int> badnodes;
          elemSet_to_nodeSet(mesh, badelems, badnodes);
          for(mt_int i=1; i<num_levels; i++) {
            nodeSet_to_elemSet(mesh, badnodes, badelems);
            elemSet_to_nodeSet(mesh, badelems, badnodes);
          }

          improve_nodeset_gradient_method(mesh, surfmesh, isSurf, badnodes, thr, sigma);

          // get updated mesh quality
          mesh_quality(mesh, elemqual);
          badelems.clear();
          get_thres_subset(elemqual, thr, badelems);
          numbadelem_after = badelems.size();

          if(numbadelem_after >= numbadelem) {
            mesh.xyz = old_xyz;
            break;
          }
        }
        std::cout << "Iter " << curit << ": Number of elements left above threshold: " << numbadelem_after << std::endl;
        print_quality_stats(elemqual);

        #ifdef WRITE_BEFORE_AFTER
        // write the initial bad elements after the mesh improvement
        std::cout << "Writing " << opts.outmsh.base + "after.vtk" << std::endl;
        for(size_t i=0; i<nod.size(); i++) {
          badmesh.xyz[i*3+0] = mesh.xyz[nod[i]*3+0];
          badmesh.xyz[i*3+1] = mesh.xyz[nod[i]*3+1];
          badmesh.xyz[i*3+2] = mesh.xyz[nod[i]*3+2];
        }
        writeVTKmesh(badmesh, opts.outmsh.base + "after.vtk");
        #endif

        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

        write_mesh = true;
      }

      if(surfmesh.xyz.size())
        surfmesh.xyz.assign(0, NULL, false);

      break;
    }

    case CLN_TOP:
    {
      size_t start_nodes = mesh.xyz.size() / 3;
      size_t start_elems = mesh.e2n_cnt.size();

      correct_duplicate_vertices(mesh, true);
      correct_duplicate_elements(mesh);

      // Remove zero volume tets / tris
      mt_vector<bool> keep(mesh.e2n_cnt.size(), false);
      mt_vector<mt_int> elems;
      size_t empty = 0;

      for(size_t eidx = 0; eidx < mesh.e2n_cnt.size(); eidx++)
      {
        mt_int estart = mesh.e2n_dsp[eidx];
        switch(mesh.etype[eidx]) {
          case Tetra:
          case Tri:
          {
            const mt_int* con = mesh.e2n_con.data() + estart;
            mt_real tetvol = volume(mesh.etype[eidx], con, mesh.xyz.data());
            if(tetvol > 0.0) {
              keep[eidx] = true;
            }
            else empty++;
            break;
          }
          default: break;
        }
      }

      if(empty) {
        printf("Removing %lu zero-volume elements..\n", empty);
        restrict_meshdata(keep, mesh, elems);
      }

      if(mesh.xyz.size() / 3 != start_nodes) write_mesh = true;
      if(mesh.e2n_cnt.size() != start_elems) write_mesh = true;
      break;
    }

    case CLN_NRML:
    {
      mt_meshdata surfmesh;

      if(opts.surf.size() == 0) {
        if(mesh.etype[0] == Tri) {
          // we assume we have a surface mesh, thus we copy the connectivity data (without
          // the coords) from mesh into surfmesh
          surfmesh = mesh;
          surfmesh.xyz.resize(0); surfmesh.xyz.reallocate();
        }
        else {
          std::cerr << "Error. Volumetric meshes require a surface file input. Aborting." << std::endl;
          exit(1);
        }
      }
      else {
        std::cout << "Reading surface.. " << std::endl;
        gettimeofday(&t1, NULL);
        readElements(surfmesh, opts.surf + SURF_EXT);
        compute_full_mesh_connectivity(surfmesh);
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      }

      surfmesh.xyz.assign(mesh.xyz.size(), mesh.xyz.data(), false);
      bool have_flipped = apply_normal_orientation(surfmesh, true);
      surfmesh.xyz.assign(0, NULL, false);

/*
      size_t flip_cnt = 0;

      for(size_t eidx = 0; eidx < surfmesh.e2n_cnt.size(); eidx++)
      {
        mt_int* con = surfmesh.e2n_con.data() + eidx*3;
        mt_int v0 = con[0], v1 = con[1], v2 = con[2];

        vec3r p0(mesh.xyz.data() + v0*3);
        vec3r p1(mesh.xyz.data() + v1*3);
        vec3r p2(mesh.xyz.data() + v2*3);
        mt_real surf = signed_tri_surf(p0, p1, p2);

        if(surf < 0) {
          con[0] = v2, con[1] = v1, con[2] = v0;
          flip_cnt++;
        }
      }
*/
      if(have_flipped) {
        if(opts.surf.size() > 0) {
          std::string filename;
          filename = opts.outmsh.base + SURF_EXT;
          std::cout << "Writing surface " << filename << std::endl;
          writeElements(surfmesh, filename);
        }
        else {
          mesh.e2n_con = surfmesh.e2n_con;
          write_mesh = true;
        }
      }

      break;
    }
  }

  if(write_mesh) {
    // we correct inside-out tets since some 3rd party tools dont keep the orientation
    std::cout << "Correcting inside-out tets .." << std::endl;
    correct_insideOut(mesh);

    std::cout << "Writing mesh: " << opts.outmsh.base << std::endl;
    gettimeofday(&t1, NULL);
    write_mesh_selected(mesh, opts.outmsh.format, opts.outmsh.base);
    gettimeofday(&t2, NULL);
    std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
  }
}

