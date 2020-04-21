/**
* @file query_mode.h
* @brief Meshtool query mode.
* @author Aurel Neic
* @version
* @date 2017-02-13
*/

#include "mt_modes_base.h"
#include "kdtree.h"

enum QRY_MODE {QRY_IDX, QRY_IDXLIST, QRY_TAGS, QRY_BBOX, QRY_EDGE, QRY_GRPH, QRY_QUAL,
               QRY_SMTH, QRY_CURV, QRY_INSIDE, QRY_NRMLS, QRY_IDXLIST_UVC};

struct query_options{
  enum QRY_MODE type;
  mt_filename msh;
  std::string coord;
  std::string uvc;
  std::string thr;
  std::string ifmt;
  std::string surf;
  std::string tags;
  std::string rad;
  std::string odat;
  std::string mode;
};


void print_query_idx_help()
{
  fprintf(stderr, "query idx: print indices in a proximity to a given coordinate\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh to query\n", mesh_par.c_str());
  fprintf(stderr, "%sx,y,z\t (input) triple of x,y,z coordinates\n", coord_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh input format. may be: %s\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "%s<surface>\t (optional) index query is restricted to this surface\n", surf_par.c_str());
  fprintf(stderr, "%s<float>\t (optional input) threshold defining the proximity to the coordinate\n", thr_par.c_str());
  fprintf(stderr, "\n");
}
void print_query_idxlist_help()
{
  fprintf(stderr, "query idxlist: apply 'query idx' for all coordinates in a file\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh to query\n", mesh_par.c_str());
  fprintf(stderr, "%s<file>\t (input) the coordinates file\n", coord_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh input format. may be: %s\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "%s<surface>\t (optional) index query is restricted to this surface\n", surf_par.c_str());
  fprintf(stderr, "\n");
}

void print_query_idxlist_uvc_help()
{
  fprintf(stderr, "query idxlist_uvc: apply 'query idx' for all coordinates in a file based on uvc coordinates\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh to perform query on\n", mesh_par.c_str());
  fprintf(stderr, "%s<file>\t (input) the uvc coordinates file holding the query points in uvc\n", coord_par.c_str());
  fprintf(stderr, "%s<file>\t (optional input) path to the ucv points file (default `basename`.uvc_pts)\n", uvc_par.c_str());
  fprintf(stderr, "\n");
}

void print_query_tags_help()
{
  fprintf(stderr, "query tags: print the tags present in a given mesh\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh to query\n", mesh_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh input format. may be: %s\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "\n");
}
void print_query_bbox_help()
{
  fprintf(stderr, "query bbox: print the bounding box of a given mesh\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh to query\n", mesh_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh input format. may be: %s\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "\n");
}
void print_query_edges_help()
{
  fprintf(stderr, "query edges: print several statistics related to the mesh edges\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh to query\n", mesh_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh input format. may be: %s\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "%stag1,tag2  (optional) List of region tags to compute the edge statistics on.\n", tags_par.c_str());
  fprintf(stderr, "%s<path>\t (optional) If set, the query data will be written to disk.\n", odat_par.c_str());
  fprintf(stderr, "\n");
}
void print_query_quality_help()
{
  fprintf(stderr, "query quality: print mesh quality statistics\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh to query\n", mesh_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh input format. may be: %s\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "%s<float>\t (optional) output exact number of elements with quality above this threshold\n", thr_par.c_str());
  fprintf(stderr, "%s<int>\t (optional) 0 = query only quality, 1 = query also self-intersection. Default is 0.\n", mode_par.c_str());
  fprintf(stderr, "%s<path>\t (optional) If set, the query data will be written to disk.\n", odat_par.c_str());
  fprintf(stderr, "\n");
}
void print_query_graph_help()
{
  fprintf(stderr, "query graph: print the nodal connectivity graph\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh to query\n", mesh_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh input format. may be: %s\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "\n");
}
void print_query_smth_help()
{
  fprintf(stderr, "query smoothness: compute the nodal smoothness\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh to query\n", mesh_par.c_str());
  fprintf(stderr, "%s<surface>\t (input) surface to compute the smoothness for\n", surf_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh input format. may be: %s\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "\n");
}
void print_query_curv_help()
{
  fprintf(stderr, "query curvature: compute the curvature of a surface.\n"
                  "In each node, a patch of adjustable size will be used to determine local curvature.\n\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh to query\n", mesh_par.c_str());
  fprintf(stderr, "%s<int>\t (input) patch radius (number of elements) for curvature estimation\n", size_par.c_str());
  fprintf(stderr, "%s<surface>\t (optional) surface to compute the curvature for\n", surf_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh input format. may be: %s\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "\n");
}
void print_query_inside_help()
{
  fprintf(stderr, "query insidepoint: get a point inside a given closed surface.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to the mesh of a closed surface\n", mesh_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh input format. may be: %s\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "\n");
}
void print_query_normals_help()
{
  fprintf(stderr, "query normals: output normal vectors for a given mesh and surface.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh to query\n", mesh_par.c_str());
  fprintf(stderr, "%s<surface>\t (optional) surface to compute the normals for\n", surf_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh input format. may be: %s\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "\n");
}

int query_parse_options(int argc, char** argv, struct query_options & opts)
{
  if(argc < 3) {
    std::cerr << "Please choose one of the following query modes: " << std::endl << std::endl;
    print_query_bbox_help();
    print_query_edges_help();
    print_query_graph_help();
    print_query_idx_help();
    print_query_idxlist_help();
    print_query_idxlist_uvc_help();
    print_query_quality_help();
    print_query_tags_help();
    print_query_smth_help();
    print_query_curv_help();
    print_query_inside_help();
    print_query_normals_help();
    return 1;
  }

  std::string query_type = argv[2];
  std::string msh_base, ifmt;

  // parse all query parameters -----------------------------------------------------------------
  for(int i=3; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh_par, msh_base);
    if(!match) match = parse_param(param, coord_par, opts.coord);
    if(!match) match = parse_param(param, uvc_par, opts.uvc);
    if(!match) match = parse_param(param, thr_par, opts.thr);
    if(!match) match = parse_param(param, surf_par, opts.surf);
    if(!match) match = parse_param(param, inp_format_par, ifmt);
    if(!match) match = parse_param(param, tags_par, opts.tags);
    if(!match) match = parse_param(param, size_par, opts.rad);
    if(!match) match = parse_param(param, odat_par, opts.odat);
    if(!match) match = parse_param(param, mode_par, opts.mode);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 2;
    }
  }
  fixBasename(opts.surf);
  opts.msh.assign(msh_base, ifmt);

  // check if all relevant parameters have been set ---------------------------------------------------
  if(query_type.compare("idx") == 0)
  {
    opts.type = QRY_IDX;

    if( !(opts.msh.isSet() && opts.coord.size() > 0) )
    {
      std::cerr << "Query index error: Insufficient parameters provided." << std::endl;
      print_query_idx_help();
      return 3;
    }
  }
  else if(query_type.compare("idxlist") == 0)
  {
    opts.type = QRY_IDXLIST;

    if( !(opts.msh.isSet() && opts.coord.size() > 0) )
    {
      std::cerr << "Query idxlist error: Insufficient parameters provided." << std::endl;
      print_query_idxlist_help();
      return 3;
    }
  }
  else if(query_type.compare("idxlist_uvc") == 0)
  {
    opts.type = QRY_IDXLIST_UVC;
    if( !(opts.msh.isSet() && opts.coord.size() > 0) )
    {
      std::cerr << "Query idxlist_uvc error: Insufficient parameters provided." << std::endl;
      print_query_idxlist_uvc_help();
      return 3;
    }
  }
  else if(query_type.compare("tags") == 0)
  {
    opts.type = QRY_TAGS;

    if(!opts.msh.isSet())
    {
      std::cerr << "Query tags error: Insufficient parameters provided." << std::endl;
      print_query_tags_help();
      return 3;
    }
  }
  else if(query_type.compare("bbox") == 0)
  {
    opts.type = QRY_BBOX;

    if(!opts.msh.isSet())
    {
      std::cerr << "Query bbox error: Insufficient parameters provided." << std::endl;
      print_query_bbox_help();
      return 3;
    }
  }
  else if(query_type.compare("edges") == 0)
  {
    opts.type = QRY_EDGE;

    if(!opts.msh.isSet())
    {
      std::cerr << "Query edges error: Insufficient parameters provided." << std::endl;
      print_query_edges_help();
      return 3;
    }
  }
  else if(query_type.compare("quality") == 0)
  {
    opts.type = QRY_QUAL;

    if(!opts.msh.isSet())
    {
      std::cerr << "Query quality error: Insufficient parameters provided." << std::endl;
      print_query_quality_help();
      return 3;
    }
  }
  else if(query_type.compare("graph") == 0)
  {
    opts.type = QRY_GRPH;

    if(!opts.msh.isSet())
    {
      std::cerr << "Query graph error: Insufficient parameters provided." << std::endl;
      print_query_graph_help();
      return 3;
    }
  }
  else if(query_type.compare("smoothness") == 0)
  {
    opts.type = QRY_SMTH;

    if( !(opts.msh.isSet() && opts.surf.size()) )
    {
      std::cerr << "Query smoothness error: Insufficient parameters provided." << std::endl;
      print_query_smth_help();
      return 3;
    }
  }
  else if(query_type.compare("curvature") == 0)
  {
    opts.type = QRY_CURV;

    if( !(opts.msh.isSet() && opts.rad.size() ) )
    {
      std::cerr << "Query curvature error: Insufficient parameters provided." << std::endl;
      print_query_curv_help();
      return 3;
    }
  }
  else if(query_type.compare("insidepoint") == 0)
  {
    opts.type = QRY_INSIDE;

    if( !opts.msh.isSet() )
    {
      std::cerr << "Query insidepoint error: Insufficient parameters provided." << std::endl;
      print_query_inside_help();
      return 3;
    }
  }
  else if(query_type.compare("normals") == 0)
  {
    opts.type = QRY_NRMLS;

    if(! opts.msh.isSet())
    {
      std::cerr << "Query normals error: Insufficient parameters provided." << std::endl;
      print_query_normals_help();
      return 3;
    }
  }
  else {
    print_usage(argv[0]);
    return 4;
  }

  return 0;
}

/// read coordinates and thresholds
template<class S>
void read_coord_file(mt_vector<S> & coords, mt_vector<S> & thr, std::string file)
{
  FILE* fd = fopen(file.c_str(), MT_FOPEN_READ);
  if(!fd) {
    treat_file_open_error(file, errno);
    exit(1);
  }

  int numcoord;
  fscanf(fd, "%d", &numcoord);

  coords.resize(numcoord*3);
  thr.resize(numcoord);

  float c[3];
  float t;

  for(int i=0; i<numcoord; i++)
  {
    fscanf(fd, "%f %f %f %f", c, c+1, c+2, &t);
    coords[i*3+0] = c[0];
    coords[i*3+1] = c[1];
    coords[i*3+2] = c[2];
    thr[i]        = t;
  }

  fclose(fd);
}

/// read coordinates and thresholds for UVC with ordering of the points in the file (z, rho, phi) -> (rho, phi, z)
template<class S>
void read_coord_file_uvc(mt_vector<S> & uvc_coords_lv, mt_vector<S> & uvc_coords_rv,
                         mt_vector<mt_int> & pos_lv, mt_vector<mt_int> & pos_rv,
                         mt_vector<S> & thr_lv, mt_vector<S> &thr_rv, mt_vector<S> & eps_lv, mt_vector<S> & eps_rv, std::string file)
{
  FILE* fd = fopen(file.c_str(), MT_FOPEN_READ);
  if(!fd) {
    treat_file_open_error(file, errno);
    exit(1);
  }

  int numcoord;
  fscanf(fd, "%d", &numcoord);

  uvc_coords_lv.resize(0);
  uvc_coords_lv.reserve(3*numcoord);
  uvc_coords_rv.resize(0);
  uvc_coords_rv.reserve(3*numcoord);
  thr_lv.resize(0);
  thr_lv.reserve(numcoord);
  thr_rv.resize(0);
  thr_rv.reserve(numcoord);
  eps_lv.resize(0);
  eps_lv.reserve(numcoord);
  eps_rv.resize(0);
  eps_rv.reserve(numcoord);
  pos_lv.resize(0);
  pos_lv.reserve(numcoord);
  pos_rv.resize(0);
  pos_rv.reserve(numcoord);

  float c[3];
  float v;
  float t;
  float eps;

  for(int i=0; i<numcoord; i++)
  {
    fscanf(fd, "%f %f %f %f %f %f", c, c+1, c+2, &v, &t, &eps);
    if(eps <= 0.0f)
      eps = std::numeric_limits<float>::max();
    if(v > 0)
    {
      uvc_coords_lv.push_back(c[1]);
      uvc_coords_lv.push_back(c[2]);
      uvc_coords_lv.push_back(c[0]);
      thr_lv.push_back(t);
      eps_lv.push_back(eps);
      pos_lv.push_back(i);
    } else {
      uvc_coords_rv.push_back(c[1]);
      uvc_coords_rv.push_back(c[2]);
      uvc_coords_rv.push_back(c[0]);
      thr_rv.push_back(t);
      eps_rv.push_back(eps);
      pos_rv.push_back(i);
    }
  }
  fclose(fd);
}

/// write computed vertex indices
template<class T>
void write_vtxlist_file(mt_vector< std::list<T> > & vtxlist, std::string file)
{
  FILE* fd = fopen(file.c_str(), MT_FOPEN_WRITE);
  if(!fd) {
    treat_file_open_error(file, errno);
    exit(1);
  }

  fprintf(fd, "%d\n", (int)vtxlist.size());

  for(size_t i=0; i<vtxlist.size(); i++)
  {
    typename std::list<T>::iterator it = vtxlist[i].begin();
    while(it != vtxlist[i].end()) {
      fprintf(fd, "%ld ", long(*it));
      it++;
    }
    fprintf(fd, "\n");
  }

  fclose(fd);
}

void print_mesh_stats(const mt_meshdata & mesh)
{
  printf("\n");
  printf("-------------- mesh statistics --------------\n");
  if(mesh.e2n_cnt.size())
    printf("Number of elements:\t%ld\n", mesh.e2n_cnt.size());
  if(mesh.xyz.size())
    printf("Number of nodes:\t%ld\n", mesh.xyz.size() / 3);
  printf("\n");
}


void query_mode(int argc, char** argv)
{
  struct query_options opts;
  int ret = query_parse_options(argc, argv, opts);
  struct timeval t1, t2;

  if(ret != 0) return;

  mt_meshdata mesh;

  switch(opts.type) {
    case QRY_IDX: {
      // query index  ------------------------------------------------------------
      mt_vector<std::string> c;
      split_string(opts.coord, ',', c);

      if(c.size() != 3) {
        std::cerr << "Error: Number of provided coords != 3 ! Aborting!" << std::endl;
        break;
      }

      float c1 = atof(c[0].c_str());
      float c2 = atof(c[1].c_str());
      float c3 = atof(c[2].c_str());

      std::cout << "Reading points of mesh: " << opts.msh.base << std::endl;
      gettimeofday(&t1, NULL);
      read_mesh_selected(mesh, opts.msh.format, opts.msh.base, CRP_READ_PTS);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      mt_meshdata surfmesh;
      mt_vector<mt_int> nod;

      if(opts.surf.size() > 0) {
        std::cout << "Restricting index search to surface: " << opts.surf + SURF_EXT << std::endl;
        readElements(surfmesh, opts.surf + SURF_EXT);
        nod = surfmesh.e2n_con;
        binary_sort(nod); unique_resize(nod);
        nod.reallocate();
      }
      else {
        nod.resize(mesh.xyz.size() / 3);
        for(size_t i=0; i<nod.size(); i++) nod[i] = i;
      }

      if(opts.thr.size() > 0)
      {
        // here we print a list of vertices that are below "thr" away from the given coordinate
        mt_real rad = atof(opts.thr.c_str());
        printf("Searching for vertices in a %.2f radius to %.2f,%.2f,%.2f \n", rad, c1,c2,c3);

        rad *= rad;

        std::vector<mt_int> vtxlist;
#if 1
        gettimeofday(&t1, NULL);
        for(const mt_int & n : nod)
        {
          mt_real k1 = c1 - mesh.xyz[n*3+0];
          mt_real k2 = c2 - mesh.xyz[n*3+1];
          mt_real k3 = c3 - mesh.xyz[n*3+2];

          mt_real dist = (k1*k1 + k2*k2 + k3*k3);

          if(dist < rad) {
            vtxlist.push_back(n);
          }
        }
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

        if(vtxlist.size()) {
          std::sort(vtxlist.begin(), vtxlist.end());
          std::cout << "Vertex list:" << std::endl;
          for(size_t i=0; i<vtxlist.size() - 1; i++)
            std::cout << vtxlist[i] << ",";
          std::cout << vtxlist[vtxlist.size() - 1] << std::endl;
        }
        else {
          std::cerr << "No vertices found!" << std::endl;
        }
#else
        vtxlist.clear();

        mt_vector<vec3r> ptsvec;
        array_to_points(mesh.xyz, ptsvec);

        kdtree tree(10), cyl_tree(10);
        tree.build_vertex_tree(ptsvec, new cartesian_csys());

        for(vec3r & v : ptsvec) {
          vec3r cyl;
          cartesian_to_cylindrical(v, cyl); v = cyl;
        }
        cyl_tree.build_vertex_tree(ptsvec, new cylindrical_csys());

        vec3f inp_pt(c1, c2, c3);
        vec3f inp_cyl_pt; cartesian_to_cylindrical(inp_pt, inp_cyl_pt);

        {
          int min_per_leaf, max_per_leaf;
          float avrg_per_leaf;

          tree.balance_stats(min_per_leaf, max_per_leaf, avrg_per_leaf);
          printf("kdtree balance stats: min = %d, max = %d, avrg = %.2f\n",
                 min_per_leaf, max_per_leaf, avrg_per_leaf);

          cyl_tree.balance_stats(min_per_leaf, max_per_leaf, avrg_per_leaf);
          printf("cylindrical kdtree balance stats: min = %d, max = %d, avrg = %.2f\n",
                 min_per_leaf, max_per_leaf, avrg_per_leaf);

          int idx;
          float len2;
          vec3f cp, ccp;
          tree.closest_vertex(inp_pt, idx, cp, len2);
          printf("cartesion kdtree: closest vertex: %d, (%f, %f, %f) \n", idx, cp.x, cp.y, cp.z);

          cyl_tree.closest_vertex(inp_cyl_pt, idx, cp, len2);
          cylindrical_to_cartesian(cp, ccp);
          printf("cylindrical kdtree: closest vertex: %d, (%f, %f, %f) \n", idx, ccp.x, ccp.y, ccp.z);
        }

        mt_vector<mixed_tuple<float, int> > vertices;
        tree.vertices_in_sphere(inp_pt, sqrt(rad), vertices);

        for(auto it = vertices.begin(); it != vertices.end(); ++it)
          vtxlist.push_back(mt_int(it->v2));

        if(vtxlist.size()) {
          std::cout << "cartesian kdtree vertex list:" << std::endl;
          for(size_t i=0; i<vtxlist.size() - 1; i++)
            std::cout << vtxlist[i] << ",";
          std::cout << vtxlist[vtxlist.size() - 1] << std::endl;
        }
        else {
          std::cerr << "No vertices found!" << std::endl;
        }

        vtxlist.clear();

        cyl_tree.vertices_in_sphere(inp_cyl_pt, sqrt(rad), vertices);

        for(auto it = vertices.begin(); it != vertices.end(); ++it)
          vtxlist.push_back(mt_int(it->v2));

        if(vtxlist.size()) {
          std::cout << "cylindrical kdtree vertex list:" << std::endl;
          for(size_t i=0; i<vtxlist.size() - 1; i++)
            std::cout << vtxlist[i] << ",";
          std::cout << vtxlist[vtxlist.size() - 1] << std::endl;
        }
        else {
          std::cerr << "No vertices found!" << std::endl;
        }

        vtxlist.clear();
        tree.k_closest_vertices(10, inp_pt, vertices);

        for(auto it = vertices.begin(); it != vertices.end(); ++it)
          vtxlist.push_back(mt_int(it->v2));

        if(vtxlist.size()) {
          // std::sort(vtxlist.begin(), vtxlist.end());
          std::cout << "Vertex list:" << std::endl;
          for(size_t i=0; i<vtxlist.size() - 1; i++)
            std::cout << vtxlist[i] << ",";
          std::cout << vtxlist[vtxlist.size() - 1] << std::endl;
        }
        else {
          std::cerr << "No vertices found!" << std::endl;
        }
#endif
      }
      else
      {
        // here we print a single vertex, which is closest to the given coordinate
        printf("Searching for vertex closest to %.2f,%.2f,%.2f \n", c1,c2,c3);
        gettimeofday(&t1, NULL);

        mt_real k1 = c1 - mesh.xyz[0];
        mt_real k2 = c2 - mesh.xyz[1];
        mt_real k3 = c3 - mesh.xyz[2];
        mt_int minIdx = 0;
        mt_real mindist = (k1*k1 + k2*k2 + k3*k3);

        for(const mt_int & n : nod)
        {
          k1 = c1 - mesh.xyz[n*3+0];
          k2 = c2 - mesh.xyz[n*3+1];
          k3 = c3 - mesh.xyz[n*3+2];

          mt_real dist = (k1*k1 + k2*k2 + k3*k3);

          if(mindist > dist) {
            mindist = dist;
            minIdx = n;
          }
        }
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

        std::cout << "Vertex Idx: " << minIdx << " Distance: " << sqrt(mindist) << std::endl;
      }
      break;
    }

    case QRY_IDXLIST: {
      // query index list  ------------------------------------------------------------
      std::cout << "Reading points of mesh: " << opts.msh.base << std::endl;
      gettimeofday(&t1, NULL);
      read_mesh_selected(mesh, opts.msh.format, opts.msh.base, CRP_READ_PTS);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      mt_meshdata surfmesh;
      mt_vector<mt_int> nod;

      if(opts.surf.size() > 0) {
        std::cout << "Restricting index search to surface: " << opts.surf + SURF_EXT << std::endl;
        readElements(surfmesh, opts.surf + SURF_EXT);
        nod = surfmesh.e2n_con;
        binary_sort(nod); unique_resize(nod);
        nod.reallocate();
      }
      else {
        nod.resize(mesh.xyz.size() / 3);
        for(size_t i=0; i<nod.size(); i++) nod[i] = i;
      }

      mt_vector<mt_real> coords, thr;
      std::cout << "Reading coord file " << opts.coord << std::endl;
      read_coord_file(coords, thr, opts.coord);

      size_t numcoord = thr.size();
      mt_vector< std::list<mt_int> > vtxlist(numcoord);

      std::cout << "Searching for vertices .." << std::endl;
      gettimeofday(&t1, NULL);

      kdtree tree(10);
      tree.build_vertex_tree(mesh.xyz);

      for(size_t i=0; i<numcoord; i++)
      {
        if(thr[i]) {
          mt_vector<mixed_tuple<mt_real, int> > vertices;
          mt_real rad = thr[i];
          vec3r ctr(coords.data() + i*3);
          tree.vertices_in_sphere(ctr, rad, vertices);
          for(auto & v : vertices)
            vtxlist[i].push_back(v.v2);
        }
        else {
          int cls_idx;
          vec3r   cls_pos;
          mt_real dist;

          vec3f ref(coords.data() + i*3);
          tree.closest_vertex(ref, cls_idx, cls_pos, dist);
          vtxlist[i].push_back(cls_idx);
        }
      }

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::string outfile = opts.coord + ".out.txt";
      write_vtxlist_file(vtxlist, outfile);
      std::cout << "Wrote " << outfile << std::endl;
      break;
    }
    case QRY_IDXLIST_UVC: {
      // query index list  ------------------------------------------------------------
      mt_vector<mt_real> xyz_lv, xyz_rv, xyz_cart;
      mt_vector<mt_int> pos_lv, pos_rv;
      mt_vector<mixed_tuple<short, size_t> > map;
      gettimeofday(&t1, NULL);
      if (opts.uvc.empty())
        opts.uvc = opts.msh.base + UVC_PTS_EXT;
      readUVCPoints(xyz_lv, xyz_rv, pos_lv, pos_rv, map, opts.uvc);
      readPoints_general(xyz_cart, opts.msh.base);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;


      mt_vector<mt_real> coords_lv, coords_rv, thr_lv, thr_rv, eps_lv, eps_rv;
      mt_vector<mt_int> coords_pos_lv, coords_pos_rv;

      std::cout << "Reading uvc coord file " << opts.coord << std::endl;
      read_coord_file_uvc(coords_lv, coords_rv, coords_pos_lv, coords_pos_rv,
                          thr_lv, thr_rv, eps_lv, eps_rv, opts.coord);

      //Need two kd-trees for each side
      gettimeofday(&t1, NULL);
      kdtree tree_cart(10);
      mt_vector<vec3r> ptsvec_cart;
      array_to_points(xyz_cart, ptsvec_cart);
      tree_cart.build_vertex_tree(ptsvec_cart);

      const size_t numcoord_lv = thr_lv.size();
      const size_t numcoord_rv = thr_rv.size();
      mt_vector< std::list<mt_int> > vtxlist(numcoord_lv + numcoord_rv);

      //left ventricle
      if (xyz_lv.size() > 0) {
        std::cout << "Searching for vertices in LV .." << std::endl;
  
        kdtree tree_lv(10);
        mt_vector<vec3r> ptsvec_lv;
        array_to_points(xyz_lv, ptsvec_lv);
        tree_lv.build_vertex_tree(ptsvec_lv);

        for(size_t i=0; i<numcoord_lv; i++)
        {
          vec3f ctr(coords_lv.data() + i*3);
          vec3r   dummy;
          mt_real dummy_len;
          int idx = 0;
          tree_lv.closest_vertex(ctr, idx, dummy, dummy_len); //Closest vertex in uvc coordinate
          //Now locate all points in the given thresholds
          if(thr_lv[i] > 0.0) {
            mt_vector<mixed_tuple<mt_real,int> > vertices_cart;
            tree_cart.vertices_in_sphere(ptsvec_cart[pos_lv[idx]], thr_lv[i], vertices_cart);
            //Throw out anything that has a rho value outside of the epsilon-box
            for(auto &v : vertices_cart)
              if(map[v.v2].v1 == 1 && std::abs(ptsvec_lv[map[v.v2].v2].x - ptsvec_lv[idx].x) < eps_lv[i])
                vtxlist[coords_pos_lv[i]].push_back(v.v2);
          }
          else
            vtxlist[coords_pos_lv[i]].push_back(pos_lv[idx]);
        }
      }
      else {
        for(size_t i=0; i<numcoord_lv; i++)
          vtxlist[coords_pos_lv[i]].push_back(-1);
      }

      //right ventricle
      if (xyz_rv.size() > 0) {
        std::cout << "Searching for vertices in RV .." << std::endl;

        kdtree tree_rv(10);
        mt_vector<vec3r> ptsvec_rv;
        array_to_points(xyz_rv, ptsvec_rv);
        tree_rv.build_vertex_tree(ptsvec_rv);

        for(size_t i=0; i<numcoord_rv; i++)
        {
          vec3r ctr(coords_rv.data() + i*3);
          vec3r dummy;
          mt_real dummy_len;
          int idx = 0;
          tree_rv.closest_vertex(ctr, idx, dummy, dummy_len); //Closest vertex in uvc coordinate
          //Now locate all points in the given thresholds
          if(thr_rv[i] > 0.0) {
            mt_vector<mixed_tuple<mt_real,int> > vertices_cart;
            tree_cart.vertices_in_sphere(ptsvec_cart[pos_rv[idx]], thr_rv[i], vertices_cart);
            //Throw out anything that has a rho value outside of the epsilon-box
            for(auto &v : vertices_cart)
              if(map[v.v2].v1 == -1 && std::abs(ptsvec_rv[map[v.v2].v2].x - ptsvec_rv[idx].x) < eps_rv[i])
                vtxlist[coords_pos_rv[i]].push_back(v.v2);
          }
          else
            vtxlist[coords_pos_rv[i]].push_back(pos_rv[idx]);
        }
      }
      else {
        for(size_t i=0; i<numcoord_rv; i++)
          vtxlist[coords_pos_rv[i]].push_back(-1);
      }

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::string outfile = opts.coord + ".out.txt";
      write_vtxlist_file(vtxlist, outfile);
      std::cout << "Wrote " << outfile << std::endl;

      break;
    }
    case QRY_TAGS: {
      // query tags  ------------------------------------------------------------
      std::cout << "Reading mesh: " << opts.msh.base << std::endl;
      gettimeofday(&t1, NULL);
      read_mesh_selected(mesh, opts.msh.format, opts.msh.base, CRP_READ_ELEM | CRP_READ_LON);
      print_mesh_stats(mesh);

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      // compute and print all unique tags
      MT_USET<mt_int> alltags;
      alltags.insert(mesh.etags.begin(), mesh.etags.end());

      std::cout << "Extracting myocardium .." << std::endl;
      gettimeofday(&t1, NULL);
      // restrict mesh to myocard
      mt_vector<mt_int> eidx;
      short numFib = mesh.lon.size() == mesh.e2n_cnt.size()*6 ? 6 : 3;

      mt_vector<bool> keep(mesh.etags.size(), false);
      #ifdef OPENMP
      #pragma omp parallel for schedule(guided, 100)
      #endif
      for(size_t i=0; i<keep.size(); i++) {
        mt_real l1 = mesh.lon[i*numFib+0];
        mt_real l2 = mesh.lon[i*numFib+1];
        mt_real l3 = mesh.lon[i*numFib+2];

        if( !(l1 < 1e-6 && l2 < 1e-6 && l3 < 1e-6) )
          keep[i] = true;
      }
      restrict_meshdata(keep, mesh, eidx);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      MT_USET<mt_int> myotags;
      myotags.insert(mesh.etags.begin(), mesh.etags.end());

      MT_USET<mt_int> bathtags;
      for(auto t : alltags)
        if(myotags.count(t) == 0) bathtags.insert(t);

      myotags.sort();
      bathtags.sort();

      size_t pidx = 0;

      std::cout << "Myocardium tags: " << std::endl;
      for(auto t : myotags) {
        std::cout << t;
        if(pidx++ < myotags.size() - 1) std::cout << ",";
      }
      std::cout << std::endl;

      pidx = 0;
      std::cout << "Bath tags: " << std::endl;
      for(auto t : bathtags) {
        std::cout << t;
        if(pidx++ < bathtags.size() - 1) std::cout << ",";
      }
      std::cout << std::endl;
      break;
    }

    case QRY_BBOX:
    {
      // query bbox  ------------------------------------------------------------
      std::cout << "Reading points of mesh: " << opts.msh.base << std::endl;
      gettimeofday(&t1, NULL);
      read_mesh_selected(mesh, opts.msh.format, opts.msh.base, CRP_READ_PTS);
      print_mesh_stats(mesh);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Computing bounding box .. " << std::endl;
      gettimeofday(&t1, NULL);
      mt_real xmin = mesh.xyz[0], xmax = mesh.xyz[0],
              ymin = mesh.xyz[1], ymax = mesh.xyz[1],
              zmin = mesh.xyz[2], zmax = mesh.xyz[2];

      for(size_t i=0; i<mesh.xyz.size()/3; i++)
      {
        mt_real x = mesh.xyz[i*3+0];
        mt_real y = mesh.xyz[i*3+1];
        mt_real z = mesh.xyz[i*3+2];

        if(xmin > x) xmin = x;
        if(ymin > y) ymin = y;
        if(zmin > z) zmin = z;

        if(xmax < x) xmax = x;
        if(ymax < y) ymax = y;
        if(zmax < z) zmax = z;
      }
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      mt_real dx = xmax - xmin, dy = ymax - ymin, dz = zmax - zmin;

      printf("Bbox:\n");
      printf("x: ( %12.2f , %12.2f ), delta %.2f\n"
             "y: ( %12.2f , %12.2f ), delta %.2f\n"
             "z: ( %12.2f , %12.2f ), delta %.2f\n", xmin,xmax,dx,ymin,ymax,dy,zmin,zmax,dz);

      break;
    }

    case QRY_EDGE:
    {
      std::cout << "Reading mesh: " << opts.msh.base << std::endl;
      gettimeofday(&t1, NULL);
      read_mesh_selected(mesh, opts.msh.format, opts.msh.base);

      if(opts.tags.size()) {

        std::cout << "Restricting mesh to tags: " << opts.tags << std::endl;

        mt_vector<std::string> tags_str;
        mt_vector<bool> keep(mesh.e2n_cnt.size(), false);
        MT_USET<mt_int> selected_tags;

        split_string(opts.tags, ',', tags_str);

        for(size_t j=0; j<tags_str.size(); j++) {
          int t = atoi(tags_str[j].c_str());
          selected_tags.insert(t);
        }

        for(size_t i=0; i<mesh.etags.size(); i++)
          if(selected_tags.count(mesh.etags[i])) keep[i] = true;

        mt_vector<mt_int> nod, eidx;
        restrict_meshdata(keep, mesh, nod, eidx);
      }
      else
        reindex_nodes(mesh);

      print_mesh_stats(mesh);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Setting up n2e / n2n graphs .. " << std::endl;
      gettimeofday(&t1, NULL);
      compute_full_mesh_connectivity(mesh);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      // build up histogramms
      size_t num_his_int = 20, num_row = 15, num_col = 90;

      std::vector<size_t> edge_cnt(mesh.n2n_cnt.size());
      for(size_t i=0; i<edge_cnt.size(); i++)
        edge_cnt[i] = mesh.n2n_cnt[i];

      if(opts.odat.size()) {
        mt_vector<size_t> conv; conv.assign(edge_cnt.begin(), edge_cnt.end());
        std::string ofile = opts.odat + ".concnt" + DAT_EXT;
        write_vector_ascii(conv, ofile, 1);
      }

      MT_USET<tuple<mt_int> > edge_set;
      tuple<mt_int> e;

      for(size_t i=0, cidx=0; i<mesh.n2n_cnt.size(); i++)
      {
        for(mt_int j=0; j < mesh.n2n_cnt[i]; j++, cidx++)
        {
          mt_int nidx =  mesh.n2n_con[cidx];
          if (nidx != mt_int(i)) {
            sortTuple((mt_int) i, nidx, e.v1, e.v2);
            edge_set.insert(e);
          }
        }
      }

      std::vector<mt_real> edge_len(edge_set.size(), 0.0);
      size_t idx=0;
      for(auto it = edge_set.begin(); it != edge_set.end(); ++it)
      {
        mt_point<mt_real> p0(mesh.xyz.data() + it->v1*3),
                        p1(mesh.xyz.data() + it->v2*3);
        mt_point<mt_real> edge = p1 - p0;

        edge_len[idx++] = edge.length();
      }

      if(opts.odat.size()) {
        mt_vector<mt_real> conv; conv.assign(edge_len.begin(), edge_len.end());
        std::string ofile = opts.odat + ".len" + DAT_EXT;
        write_vector_ascii(conv, ofile, 1);
      }

      mt_vector<mt_real> elemvolumes;
      std::cout << "Computing element sizes .." << std::endl;
      gettimeofday(&t1, NULL);
      //mesh_min_angles(mesh, elemvolumes);
      mesh_volumes(mesh, mesh.xyz, elemvolumes);
      gettimeofday(&t2, NULL);

      std::vector<mt_real> vols;
      vols.assign(elemvolumes.begin(), elemvolumes.end());

      if(opts.odat.size()) {
        std::string ofile = opts.odat + ".vol" + DAT_EXT;
        write_vector_ascii(elemvolumes, ofile, 1);
      }

      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      histogramm_plotter<size_t> ecnt_histo(num_his_int, num_row, num_col);
      histogramm_plotter<mt_real> elen_histo(num_his_int, num_row, num_col);
      histogramm_plotter<mt_real> vol_histo(num_his_int, num_row, num_col);

      ecnt_histo.plot(edge_cnt, "\nNumber of connections");
      elen_histo.plot(edge_len, "\nEdge lengths");
      vol_histo.plot(vols, "\nElement volumes");
      break;
    }

    case QRY_QUAL:
    {
      std::cout << "Reading mesh: " << opts.msh.base << std::endl;
      gettimeofday(&t1, NULL);

      read_mesh_selected(mesh, opts.msh.format, opts.msh.base, CRP_READ_ELEM | CRP_READ_PTS);
      mesh.e2n_dsp.resize(mesh.e2n_cnt.size());
      bucket_sort_offset(mesh.e2n_cnt, mesh.e2n_dsp);

      print_mesh_stats(mesh);

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      int mode = 0;
      if(opts.mode.size()) mode = atoi(opts.mode.c_str());

      if(mode == 1) {
        std::cout << "Checking self intersection .." << std::endl;
        gettimeofday(&t1, NULL);

        mt_vector<tuple<mt_int>> intersec_edges;
        mt_vector<mt_int>        intersec_eidx;

        check_self_intersection(mesh, intersec_edges, intersec_eidx);

        if(intersec_edges.size() > 0) {
          mt_meshdata intersec;

          for(size_t i=0; i<intersec_edges.size(); i++) {
            const tuple<mt_int> & edge = intersec_edges[i];
            const mt_int & eidx = intersec_eidx[i];

            mt_int lncon[2] = {edge.v1, edge.v2};
            mt_int curtag = intersec.e2n_cnt.size();

            mesh_add_elem(intersec, Line, lncon, curtag);
            mesh_add_elem(intersec, mesh.etype[eidx], mesh.e2n_con.data() + mesh.e2n_dsp[eidx], curtag);
          }

          intersec.xyz = mesh.xyz;

          mt_vector<mt_int> nod;
          reindex_nodes(intersec, nod);

          const char* intersec_basename = "intersec.mesh";
          fprintf(stderr, "Self intersections detected. Outputting mesh intersections into mesh: %s.vtk \n", intersec_basename);
          write_mesh_selected(intersec, "vtk_bin", intersec_basename);
        }

        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      }

      mt_vector<mt_real> elemqual;
      std::cout << "Computing mesh quality .." << std::endl;
      gettimeofday(&t1, NULL);
      mesh_quality(mesh, elemqual);

      if(opts.odat.size()) {
        std::string ofile = opts.odat + ".qual" + DAT_EXT;
        write_vector_ascii(elemqual, ofile, 1);
      }

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;


      mt_real min = elemqual[0], max = 0;
      mt_real avrg = 0.0, stddev = 0.0;

      bool comp_thr = false;
      float thr;
      int above_thr = 0;

      if(opts.thr.size() > 0) {
        thr = atof(opts.thr.c_str());
        comp_thr = true;
      }

      // get basic min-max
      for(size_t i=0; i<elemqual.size(); i++)
      {
        mt_real c = elemqual[i];
        if(min > c) min = c;
        if(max < c) max = c;
        if( comp_thr && (c > thr) ) above_thr++;
        avrg += c;
      }
      avrg /= (double) elemqual.size();

      for(size_t i=0; i< elemqual.size(); i++) {
        mt_real c = elemqual[i];
        stddev += (c - avrg) * (c - avrg);
      }
      stddev /= ((double)(elemqual.size()) -1.0);
      stddev = std::sqrt(stddev);

      // build up histogramm
      size_t numint = 50;
      float delta = 1.0f / numint;
      std::vector<float> xval(numint+1), yval(numint+1, 0);

      asciiPlotter plotter(22, 90);

      for(size_t i=0; i<numint+1; i++) xval[i] = i*delta;
      for(size_t i=0; i<elemqual.size(); i++) {
        size_t idx = elemqual[i] / delta;
        if(idx < numint+1)
          yval[idx] += 1.0f;
      }
      for(size_t i=0; i<numint+1; i++)
        yval[i] /= (float)elemqual.size();

      plotter.set_xrange(0.0, 1.0);
      plotter.set_yrange(yval);
      plotter.add_graph(xval, yval, '*');
      std::cout << std::endl << "--------- element quality statistics ---------" << std::endl;
      std::cout << "Used method: " << QMETRIC_NAME << " with 0 == best, 1 == worst"<< std::endl;
      std::cout << "Min: " << min << " Max: " << max << " Mean: " << avrg << " Stddev: " << stddev << std::endl;
      if(comp_thr)
        std::cout << "Number of elements with quality above " << thr << " : " << above_thr << std::endl;

      std::cout << std::endl << "Histogramm of element quality:" << std::endl << std::endl;
      plotter.print();
      break;
    }

    case QRY_GRPH:
    {
      std::cout << "Reading mesh: " << opts.msh.base << std::endl;
      gettimeofday(&t1, NULL);

      read_mesh_selected(mesh, opts.msh.format, opts.msh.base, CRP_READ_ELEM );
      print_mesh_stats(mesh);
      compute_full_mesh_connectivity(mesh);

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      size_t num_his_int = 10, num_row = 15, num_col = 60;
      histogramm_plotter<mt_int> bw_histo(num_his_int, num_row, num_col);

      std::vector<mt_int> bwdth(mesh.n2n_cnt.size());
      for(size_t i=0, k=0; i<mesh.n2n_cnt.size(); i++)
      {
        mt_int max=0, min=mesh.n2n_cnt.size();
        for(mt_int j=0; j<mesh.n2n_cnt[i]; j++, k++)
        {
          mt_int c = mesh.n2n_con[k];
          if(min > c) min = c;
          if(max < c) max = c;
        }
        bwdth[i] = max - min;
      }

      bw_histo.plot(bwdth, "\n ------- matrix bandwidth ------- ");

      std::cout << "\n\n ------- mesh connectivity graph pattern -------" <<
      std::endl << std::endl;

      matrixgraph_plotter<mt_vector<mt_int> > plotter(30, 60);
      plotter.print(mesh.n2n_cnt, mesh.n2n_con, '*');
      std::cout << std::endl;

      // std::string outfile = opts.msh.base + ".grph.bin";
      // std::cout << "writing connectivity matrix " << outfile << std::endl;
      // write_graph(mesh.n2n_cnt, mesh.n2n_con, outfile.c_str());
      break;
    }

    case QRY_SMTH:
    {
      mt_meshdata surfmesh;

      std::cout << "Reading mesh: " << opts.msh.base << std::endl;
      gettimeofday(&t1, NULL);

      read_mesh_selected(mesh, opts.msh.format, opts.msh.base, CRP_READ_ELEM | CRP_READ_PTS);
      print_mesh_stats(mesh);
      compute_full_mesh_connectivity(mesh);

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Reading surface.. " << std::endl;
      gettimeofday(&t1, NULL);

      readElements(surfmesh, opts.surf + SURF_EXT);
      compute_full_mesh_connectivity(surfmesh);

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Computing smoothness .. " << std::endl;
      gettimeofday(&t1, NULL);

      size_t nnodes = mesh.n2n_cnt.size();
      mt_vector<double> smoothness(nnodes, 0.0);
      mt_vector<bool>   onManifold(nnodes, false);

      for(size_t i=0; i<surfmesh.e2n_con.size(); i++) onManifold[surfmesh.e2n_con[i]] = true;

      for(size_t i=0; i<nnodes; i++) {
        // smoothness[i] = distance_to_centroid(mesh, surfmesh, onManifold, i);

        if(onManifold[i])
          //smoothness[i] = nbhd_smoothness_functional(surfmesh, mesh.xyz, i, 0.9);
          smoothness[i] = minimal_surf_normal_correlation(surfmesh, mesh.xyz, i);
        else
          smoothness[i] = 0.0;
      }

      std::string outfile = opts.msh.base + ".smth" + DAT_EXT;
      write_vector_ascii(smoothness, outfile);

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      break;
    }
    case QRY_CURV:
    {
      mt_meshdata surfmesh;

      std::cout << "Reading mesh: " << opts.msh.base << std::endl;
      gettimeofday(&t1, NULL);

      read_mesh_selected(mesh, opts.msh.format, opts.msh.base, CRP_READ_ELEM | CRP_READ_PTS);
      print_mesh_stats(mesh);

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

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
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      }
      compute_full_mesh_connectivity(surfmesh);

      std::cout << "Computing curvature .. " << std::endl;
      gettimeofday(&t1, NULL);

      size_t nnodes = mesh.xyz.size() / 3;
      mt_vector<double> curvature(nnodes, 0.0);
      mt_vector<bool>   onManifold(nnodes, false);

      for(size_t i=0; i<surfmesh.e2n_con.size(); i++) onManifold[surfmesh.e2n_con[i]] = true;

      int rad = atoi(opts.rad.c_str());

      #ifdef OPENMP
      #pragma omp parallel for schedule(dynamic)
      #endif
      for(size_t nidx=0; nidx<nnodes; nidx++) {
        if(onManifold[nidx])
        {
          mt_vector<vec3r> vert;

          MT_USET<mt_int> nodes, elems;
          nodes.insert(nidx);

          for(int i=0; i<rad; i++) {
            nodeSet_to_elemSet(surfmesh, nodes, elems);
            elemSet_to_nodeSet(surfmesh, elems, nodes);
          }
          vert.resize(nodes.size());
          mt_int idx=0;
          for(const mt_int & n : nodes) vert[idx++] = vec3r(mesh.xyz.data() + n*3);

          mt_real R;
          vec3r pos;
          fit_sphere(vert, R, pos);

          #if 0
          mt_real Rb;
          bounding_sphere(vert, pos, Rb);
          R /= Rb;
          #endif

          curvature[nidx] = 1.0 / R;
        }
        else
          curvature[nidx] = 0.0;
      }

      std::string outfile = opts.msh.base + ".curv" + DAT_EXT;
      write_vector_ascii(curvature, outfile);

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      break;
    }

    case QRY_INSIDE:
    {
      std::cout << "Reading mesh: " << opts.msh.base << std::endl;
      gettimeofday(&t1, NULL);

      read_mesh_selected(mesh, opts.msh.format, opts.msh.base, CRP_READ_ELEM | CRP_READ_PTS);
      compute_full_mesh_connectivity(mesh);
      print_mesh_stats(mesh);

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;


      if(mesh.etype[0] == Tri)
      {
        float avg_edge_length = avrg_edgelength_estimate(mesh);
        kdtree tree(10);
        tree.build_tree(mesh);

        bool did_find = false;
        vec3f found_point(0,0,0);
        int iter = 0, max_iter = 10;

        do {
          mt_int eidx = drand48() * mesh.e2n_cnt.size();
          mt_int v1 = mesh.e2n_con[eidx*3+0];
          mt_int v2 = mesh.e2n_con[eidx*3+1];
          mt_int v3 = mesh.e2n_con[eidx*3+2];
          vec3f inside_normal = triangle_normal(v1, v2, v3, mesh.xyz) * mt_real(-1.0);
          vec3f centerpoint   = triangle_centerpoint(v1, v2, v3, mesh.xyz);

          vec3f insidepoint = centerpoint + inside_normal * avg_edge_length;
          if(inside_closed_surface(tree, insidepoint)) {
            did_find = true;
            found_point = insidepoint;
          }
          else {
            inside_normal *= -1.0f;
            insidepoint = centerpoint + inside_normal * avg_edge_length;
            if(inside_closed_surface(tree, insidepoint)) {
              did_find = true;
              found_point = insidepoint;
            }
          }
        } while(did_find == false && iter < max_iter);

        if(did_find) {
          printf("Point inside closed surface:\n%.2f,%.2f,%.2f\n",
                 found_point.x, found_point.y, found_point.z);
        }
        else {
          fprintf(stderr, "Could not find point inside closed surface!\n");
          exit(1);
        }
      }
      else {
        fprintf(stderr, "query inside error: This mode only operates on surface meshes! Aborting!\n");
        exit(1);
      }
      break;
    }

    case QRY_NRMLS:
    {
      mt_meshdata surfmesh;

      std::cout << "Reading mesh: " << opts.msh.base << std::endl;
      gettimeofday(&t1, NULL);

      read_mesh_selected(mesh, opts.msh.format, opts.msh.base, CRP_READ_ELEM | CRP_READ_PTS);
      print_mesh_stats(mesh);
      compute_full_mesh_connectivity(mesh);

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

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

      std::cout << "Computing nodal surface normals .. " << std::endl;
      mt_vector<mt_real> nrmls, snrmls;
      compute_nodal_surface_normals(surfmesh, mesh.xyz, snrmls);

      nrmls.resize(mesh.xyz.size(), 0.0);

      for(size_t i=0; i < snrmls.size() / 3; i++) {
        if(surfmesh.n2e_cnt[i] > 0) {
          nrmls[i*3+0] = snrmls[i*3+0];
          nrmls[i*3+1] = snrmls[i*3+1];
          nrmls[i*3+2] = snrmls[i*3+2];
        }
      }

      std::string filename = opts.msh.base + ".normals" + VEC_EXT;
      write_vector_ascii(nrmls, filename, 3);
      break;
    }
  }
}

