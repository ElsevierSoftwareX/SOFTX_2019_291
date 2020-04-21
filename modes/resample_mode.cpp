/**
* @file resample_mode.h
* @brief Refine mesh elements
* @author Aurel Neic
* @version
* @date 2017-02-13
*/

#include "mt_modes_base.h"

/**
* @brief Enum encoding the resample mode.
*/
enum RES_MODE {RES_PURK, RES_MESH, RES_SURF};

static const std::string pkje_in_par   = "-ips=";
static const std::string pkje_out_par  = "-ops=";
static const std::string surf_corr_par = "-surf_corr=";
static const std::string postsmth_par  = "-postsmth=";
static const std::string uni_par       = "-uniform=";
static const std::string conv_par      = "-conv=";

/**
* @brief The resample mode options struct.
*/
struct resample_options{
  enum RES_MODE mode;
  std::string pkje_in;
  std::string pkje_out;
  std::string msh_base;
  std::string outmsh_base;
  std::string ifmt;
  std::string ofmt;
  std::string size;
  std::string min;
  std::string max;
  std::string avrg;
  std::string tags;
  std::string surf_corr;
  std::string fix;
  std::string postsmth;
  std::string unif;
  std::string conv;
};

#define CORR_DFLT 0.95

/**
* @brief The resample purkinje help message.
*/
void print_resample_purkinje_help()
{
  fprintf(stderr, "resample purkinje: resample purkinje cables as close as possible to a segment size\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) Input Purkinje System file.\n", pkje_in_par.c_str());
  fprintf(stderr, "%s<float>\t (input) segment size to resample towards\n", size_par.c_str());
  fprintf(stderr, "%s<path>\t (output) Output Purkinje System file.\n", pkje_out_par.c_str());
  fprintf(stderr, "\n");
}
void print_resample_mesh_help()
{
  fprintf(stderr, "resample mesh: resample a tetrahedral mesh to fit a given edge size range\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t (input) path to basename of the input mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<float>\t\t (input) min edge size\n", min_par.c_str());
  fprintf(stderr, "%s<float>\t\t (input) max edge size\n", max_par.c_str());
  fprintf(stderr, "%s<float>\t\t (optional) average edge size. (min, max) are derived from it.\n", avrg_par.c_str());
  fprintf(stderr, "%s<path>\t\t (output) path to basename of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "%s<float>\t (optional) surface correlation parameter (default %.2f).\n"
                  "\t\t\t    Edges connecting surface nodes with surface normals correlating greater\n"
                  "\t\t\t    than the specified amount can be collapsed.\n", surf_corr_par.c_str(), CORR_DFLT);
  fprintf(stderr, "%s<int>\t\t (optional) Apply conservative smoothing between resampling iterations.\n"
                  "\t\t\t    0 = no, 1 = yes. Default is 1.\n", postsmth_par.c_str());
  fprintf(stderr, "%s<int>\t\t (optional) Edge-refinement is applied uniformly on selected tags.\n"
                  "\t\t\t    0 = no, 1 = yes. Default is 0.\n", uni_par.c_str());
  fprintf(stderr, "%s<int>\t\t (optional) Fix boundary of mesh. 0 = no, 1 = yes. Default is 0.\n",
                  fix_par.c_str());
  fprintf(stderr, "%s<int>\t\t (optional) Iterate until fully converged. 0 = no, 1 = yes. Default is 0.\n",
                  conv_par.c_str());
  fprintf(stderr, "%s<tag lists>\t (optional) element tag lists specifying the regions\n"
                  "\t\t\t    to perform resampling on\n", tags_par.c_str());
  fprintf(stderr, "%s<path>\t\t (optional) format of the input mesh\n", inp_format_par.c_str());
  fprintf(stderr, "%s<path>\t\t (optional) format of the output mesh\n\n", out_format_par.c_str());
  fprintf(stderr, "The supported input formats are:\n%s\n", input_formats.c_str());
  fprintf(stderr, "The supported output formats are:\n%s\n", output_formats.c_str());
  fprintf(stderr, "\n");
  fprintf(stderr, "The tag lists have the following syntax: multiple lists are seperated by a \"/\",\n"
                  "while the tags in one list are seperated by a \",\" character. The surface smoothness\n"
                  "of the submesh formed by one tag list will be preserved.\n");
  fprintf(stderr, "\n");
}
void print_resample_surfmesh_help()
{
  fprintf(stderr, "resample surfmesh: resample a triangle mesh to fit a given edge size range\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t (input) path to basename of the input mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<float>\t\t (input) min edge size\n", min_par.c_str());
  fprintf(stderr, "%s<float>\t\t (input) max edge size\n", max_par.c_str());
  fprintf(stderr, "%s<float>\t\t (optional) average edge size. (min, max) are derived from it.\n", avrg_par.c_str());
  fprintf(stderr, "%s<path>\t\t (output) path to basename of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "%s<float>\t (optional) surface correlation parameter (default %.2f).\n"
                  "\t\t\t    Edges connecting surface nodes with surface normals correlating greater\n"
                  "\t\t\t    than the specified amount can be collapsed.\n", surf_corr_par.c_str(), CORR_DFLT);
  fprintf(stderr, "%s<int>\t\t (optional) Fix boundary of non-closed surfaces. 0 = no, 1 = yes. Default is 1.\n",
                  fix_par.c_str());
  fprintf(stderr, "%s<int>\t\t (optional) Apply conservative post-resampling smoothing.\n"
                  "\t\t\t    0 = no, 1 = yes. Default is 1.\n", postsmth_par.c_str());
  fprintf(stderr, "%s<int>\t\t (optional) Edge-refinement is applied uniformly on selected tags.\n"
                  "\t\t\t    0 = no, 1 = yes. Default is 0.\n", uni_par.c_str());
  fprintf(stderr, "%s<format>\t\t (optional) format of the input mesh\n", inp_format_par.c_str());
  fprintf(stderr, "%s<format>\t\t (optional) format of the output mesh\n", out_format_par.c_str());
  fprintf(stderr, "%s<tag lists>\t (optional) element tag lists specifying the regions\n"
                  "\t\t\t    to perform resampling on\n\n", tags_par.c_str());
  fprintf(stderr, "The supported input formats are:\n%s\n", input_formats.c_str());
  fprintf(stderr, "The supported output formats are:\n%s\n", output_formats.c_str());
  fprintf(stderr, "\n");
  fprintf(stderr, "The tag lists have the following syntax: multiple lists are seperated by a \"/\",\n"
                  "while the tags in one list are seperated by a \",\" character. The surface smoothness\n"
                  "of the submesh formed by one tag list will be preserved.\n");
  fprintf(stderr, "\n");
}

/**
* @brief Refine mode options parsing
*/
int resample_parse_options(int argc, char** argv, struct resample_options & opts)
{
  if(argc < 3) {
    std::cerr << "Please choose one of the following resample modes: " << std::endl << std::endl;
    print_resample_purkinje_help();
    print_resample_mesh_help();
    print_resample_surfmesh_help();
    return 1;
  }

  std::string resample_type = argv[2];

  // parse all resample parameters -----------------------------------------------------------------
  for(int i=3; i<argc; i++)
  {
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, pkje_in_par, opts.pkje_in);
    if(!match) match = parse_param(param, pkje_out_par, opts.pkje_out);
    if(!match) match = parse_param(param, mesh_par, opts.msh_base);
    if(!match) match = parse_param(param, outmesh_par, opts.outmsh_base);
    if(!match) match = parse_param(param, inp_format_par, opts.ifmt);
    if(!match) match = parse_param(param, out_format_par, opts.ofmt);
    if(!match) match = parse_param(param, size_par, opts.size);
    if(!match) match = parse_param(param, min_par, opts.min);
    if(!match) match = parse_param(param, max_par, opts.max);
    if(!match) match = parse_param(param, avrg_par, opts.avrg);
    if(!match) match = parse_param(param, tags_par, opts.tags);
    if(!match) match = parse_param(param, surf_corr_par, opts.surf_corr);
    if(!match) match = parse_param(param, fix_par, opts.fix);
    if(!match) match = parse_param(param, postsmth_par, opts.postsmth);
    if(!match) match = parse_param(param, uni_par, opts.unif);
    if(!match) match = parse_param(param, conv_par, opts.conv);
    fixBasename(opts.msh_base);
    fixBasename(opts.outmsh_base);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 3;
    }
  }

  if(resample_type.compare("purkinje") == 0)
  {
    opts.mode = RES_PURK;
    // check if all relevant parameters have been set ---------------------------------------------------
    bool inok = opts.pkje_in.size() > 0, outok = opts.pkje_out.size() > 0,
         throk = opts.size.size() > 0;

    if( !(inok && outok && throk) )
    {
      std::cerr << "Error: Insufficient parameters provided." << std::endl;
      print_resample_purkinje_help();
      return 2;
    }
  }
  else if(resample_type.compare("mesh") == 0)
  {
    opts.mode = RES_MESH;

    bool inok = opts.msh_base.size() > 0, outok = opts.outmsh_base.size() > 0;

    if( !(inok && outok && (opts.min.size() + opts.max.size() + opts.avrg.size()) > 0) )
    {
      std::cerr << "Error: Insufficient parameters provided." << std::endl;
      print_resample_mesh_help();
      return 2;
    }
  }
  else if(resample_type.compare("surfmesh") == 0)
  {
    opts.mode = RES_SURF;

    bool inok = opts.msh_base.size() > 0, outok = opts.outmsh_base.size() > 0;

    if( !(inok && outok &&  (opts.min.size() + opts.max.size() + opts.avrg.size()) > 0) )
    {
      std::cerr << "Error: Insufficient parameters provided." << std::endl;
      print_resample_surfmesh_help();
      return 2;
    }
  }
  else {
    print_usage(argv[0]);
    return 4;
  }

  return 0;
}

/**
* @brief Check if the points of a cable are located along a line.
*
* @param [in] cab  The cable.
* @param [in] thr  Threshold for the scalar product testing alignment.
*                  Should be in (0,1), with 0 = orthogonality, and 1 = identical alignment
*
* @return True for line, else False.
*/
bool cable_is_line(const mt_pscable & cab, mt_real thr = 0.9)
{
  size_t num_inp_pts = cab.pts.size() / 3;
  mt_point<mt_real> s(cab.pts.data()), e(cab.pts.data() + (num_inp_pts-1)*3);
  mt_point<mt_real> d = e - s; d.normalize();

  for(size_t i=0; i<num_inp_pts; i++)
  {
    mt_point<mt_real> p = mt_point<mt_real>(cab.pts.data() + i*3) - s;
    p.normalize();
    mt_real sca = p.scaProd(d);
    if(sca < thr) return false;
  }

  return true;
}

/**
* @brief Descretize line cables of a PS to as many intervals as required for
*        an interval size matching given size parameter.
*
* @param [in, out] ps    Purkinje system
* @param [in]      size  Required interval size.
*/
void resample_purkinje(struct mt_psdata & ps, float size)
{
  for(size_t cidx = 0; cidx < ps.cables.size(); cidx++)
  {
    struct mt_pscable & ccab = ps.cables[cidx];

    // jump over cables that are not a geometric line
    bool isLine = cable_is_line(ccab, mt_real(0.9));
    if(!isLine) continue;

    size_t num_inp_pts = ccab.pts.size() / 3;
    mt_point<mt_real> s(ccab.pts.data()), e(ccab.pts.data() + (num_inp_pts-1)*3);
    mt_point<mt_real> d = e - s;

    // determine segment length and number of required intervals
    mt_real len = d.length();
    int numint = ((int)len + (int)size - 1) / (int)size;
    // scale segment vector
    d /= (mt_real) numint;

    int numpts = numint + 1;
    ccab.pts.resize( numpts*3);
    // start point
    ccab.pts[0] = s.x;
    ccab.pts[1] = s.y;
    ccab.pts[2] = s.z;
    // intermediate points
    for(int iidx = 1; iidx < numpts - 1; iidx++)
    {
      ccab.pts[iidx*3+0] = s.x + iidx*d.x;
      ccab.pts[iidx*3+1] = s.y + iidx*d.y;
      ccab.pts[iidx*3+2] = s.z + iidx*d.z;
    }
    // end point
    ccab.pts[(numpts-1)*3+0] = e.x;
    ccab.pts[(numpts-1)*3+1] = e.y;
    ccab.pts[(numpts-1)*3+2] = e.z;
  }
}

void preprocess_surf_tags(const mt_meshdata & mesh,
                          const std::string & oper_string,
                          mt_vector<MT_USET<mt_int> > & stags)
{
  // pre-process operations
  mt_vector<std::string> oper;
  if(oper_string.size() > 0)
    split_string(oper_string, '/' , oper);
  else {
    oper.resize(1);
    oper[0] = "*";
  }

  stags.resize(oper.size());
  // extract tags of all operations
  for(int oidx = 0; oidx < (int)oper.size(); oidx++)
  {
    mt_vector<std::string> tags;
    // either read in from user input or add all tags
    if(oper[oidx].compare("*") != 0) {
      split_string(oper[oidx], ',' , tags);
      for(const std::string & s : tags)
        stags[oidx].insert(atoi(s.c_str()));
    }
    else
      stags[oidx].insert(mesh.etags.begin(), mesh.etags.end());

    stags[oidx].sort();
  }
}


void preprocess_min_max(const std::string & minstr, const std::string & maxstr,
                        const std::string & avrgstr, const float avrg_est,
                        float & min, float & max)
{
  min = -1.0;
  max = -1.0;

  if(avrgstr.size()) {
    float avrg_req = atof(avrgstr.c_str());

    if(avrg_est * 0.5f > avrg_req) {
      // the requested edge length is significantly smaller than the current one
      // thus, we will be mainly refining
      min = 0.5f * avrg_req, max = 1.8f * avrg_req;
    }
    else if(avrg_est * 1.5f < avrg_req) {
      // coarsening will play a significant role. Since it has lesser impact on the
      // average, we move it closer to it.
      min = 0.75f * avrg_req, max = 1.8f * avrg_req;
    }
    else {
      min = 0.65f * avrg_req, max = 1.7f * avrg_req;
    }

    printf("Selected average edge length %.2f. Setting min = %.2f, max = %.2f.\n",
           avrg_req, min, max);
    min *= min;
    max *= max;
  }
  else {
    if(minstr.size()) {
      min = atof(minstr.c_str());
      min *= min;
    }
    if(maxstr.size()) {
      max = atof(maxstr.c_str());
      max *= max;
    }
  }
}


#define SMTH_BASE 0.14f
#define SMTH_CHANGE 0.95f
#define EDGE_ANG 30.0f
#define ITER_BASE 40
#define ITER_THRS 3

bool non_uniform_splitting(mt_meshdata & mesh,
                           const MT_USET<mt_int> & tags,
                           mt_edge_splitter & splitter,
                           mt_edge_manager & edge_manager,
                           const mt_vector<bool> & fixed_nod,
                           const mt_real max)
{
  timeval t1, t2;
  size_t iter = 0, nsel = 0;
  char msg[1024];

  bool having_fixed_nodes = fixed_nod.size();
  bool updated_mesh = false;

  mt_vector<edgeele> elem_edge_lengths;
  edge_len_asc ascending;

  // main splitting loop
  do {
    std::set<edgeele, edge_len_dsc> sel_edges;
    const MT_MAP<tuple<mt_int>,mt_int> & edges = edge_manager.edges();

    for(size_t eidx=0; eidx < mesh.e2n_cnt.size(); eidx++)
    {
      if(tags.count(mesh.etags[eidx])) {
        add_edges_lengths(mesh, eidx, elem_edge_lengths);

        edgeele longest = *std::max_element(elem_edge_lengths.begin(), elem_edge_lengths.end(), ascending);
        bool edge_is_fixed = false;
        if(having_fixed_nodes) {
          bool v1_fixed = fixed_nod[longest.v1],
               v2_fixed = fixed_nod[longest.v2];
          // we cannot touch edges going to the fixed surface. first, because sometimes vertices
          // switch and we might brake the logic, second for elem quality reasons
          edge_is_fixed = v1_fixed || v2_fixed;
        }

        if(!edge_is_fixed && longest.length > max) {
          longest.split_at = 0.5f;

          auto it = edges.find({longest.v1, longest.v2});
          longest.idx = it->second;

          sel_edges.insert(longest);
        }
      }
    }

    nsel = sel_edges.size();
    iter++;
    if(nsel) {
      sprintf(msg, "Iter %lu: splitting %lu edges, %.4f %% of total %lu ",
              iter, nsel, float(nsel) / edges.size() * 100.0f, edges.size());
      PROGRESS<size_t> prog(nsel, msg);
      gettimeofday(&t1, NULL);

      splitter(sel_edges, prog);

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      updated_mesh = true;
    }
  }
  while(nsel > 0);

  return updated_mesh;
}

bool uniform_splitting(mt_meshdata & mesh,
                       const MT_USET<mt_int> & tags,
                       mt_edge_manager & edge_manager,
                       const mt_real max)
{
  timeval t1, t2;
  size_t iter = 0, nsel = 0;
  bool updated_mesh = false;

  // main splitting loop
  do {
    std::set<edgeele, edge_len_dsc> sel_edges;
    // this is the set of edges for all tags
    const MT_MAP<tuple<mt_int>,mt_int> & edges = edge_manager.edges();
    mt_vector<edgeele> elem_edge_lengths;
    mt_real avrg_len = 0.0, num_added = 0.0;

    // first we determine if the average edge length of the specified tag region is
    // significantly larger than the specified max length
    for(size_t eidx=0; eidx < mesh.e2n_cnt.size(); eidx++)
    {
      if(tags.count(mesh.etags[eidx]))
        add_edges_lengths(mesh, eidx, elem_edge_lengths);
    }

    for(const edgeele & e : elem_edge_lengths) {
      avrg_len  += e.length;
      num_added += 1.0;
    }

    if(num_added)
      avrg_len /= num_added;

    if(avrg_len * 0.6 > max) {
      iter++;
      gettimeofday(&t1, NULL);
      printf("Iter %lu: Uniform refining .. \n", iter);

      refine_uniform(mesh, edges);
      bucket_sort_offset(mesh.e2n_cnt, mesh.e2n_dsp);
      edge_manager.refresh_edge_data(false);

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      updated_mesh = true;
    }
  }
  while(nsel > 0);

  return updated_mesh;
}

void update_surface_data(mt_meshdata & mesh,
                         const MT_USET<mt_int> & tags,
                         const bool fix_bnd,
                         mt_vector<bool> & surf_nod,
                         mt_vector<bool> & boundary_node,
                         mt_vector<mt_real> & surf_nrml)
{
  mt_meshdata surfmesh;
  compute_surface(mesh, tags, surfmesh, true);
  compute_full_mesh_connectivity(surfmesh);

  surf_nod.assign(mesh.n2e_cnt.size(), false);
  if(fix_bnd) boundary_node.assign(mesh.n2e_cnt.size(), false);

  for(auto c : surfmesh.e2n_con) {
    surf_nod[c] = true;
    if(fix_bnd) boundary_node[c] = true;
  }

  surfmesh.xyz.assign(mesh.xyz.size(), mesh.xyz.data(), false);
  compute_nodal_surface_normals(surfmesh, mesh.xyz, surf_nrml);
  surfmesh.xyz.assign(0, NULL, false);
}


void resample_purk_mode(const resample_options & opts)
{
  struct timeval t1, t2;
  struct mt_psdata ps;

  std::cout << "Reading purkinje file: " << opts.pkje_in << std::endl;
  gettimeofday(&t1, NULL);
  read_purkinje(ps, opts.pkje_in);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  float size = atof(opts.size.c_str());
  std::cout << "Refining PS .. " << std::endl;
  gettimeofday(&t1, NULL);
  resample_purkinje(ps, size);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  std::cout << "Writing purkinje file: " << opts.pkje_out << std::endl;
  gettimeofday(&t1, NULL);
  write_purkinje(ps, opts.pkje_out);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
}

void resample_mesh_mode(resample_options & opts)
{
  struct timeval t1, t2;

  bool fix_bnd   = opts.fix.size() ? atoi(opts.fix.c_str()) : 0;
  bool unif_splt = opts.unif.size() ? atoi(opts.unif.c_str()) : 0;
  bool do_conv   = opts.conv.size() ? atoi(opts.conv.c_str()) : 0;

  if(fix_bnd && unif_splt) {
    fprintf(stderr, "Warning: Boundary fixing and uniform refinement exclude each other!\n"
        "Turning off uniform refinement!\n");
    unif_splt = false;
  }

  if(unif_splt && opts.tags.size()) {
    fprintf(stderr, "Warning: Uniform refinement currently does not support local refinement! Ignoring tags!\n");
    opts.tags = "";
  }

  mt_meshdata mesh;
  mt_filename msh(opts.msh_base, opts.ifmt);
  std::cout << "Reading mesh: " << msh.base << std::endl;
  gettimeofday(&t1, NULL);
  read_mesh_selected(mesh, msh.format, msh.base);
  compute_full_mesh_connectivity(mesh, msh.base);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  float min, max;
  float avrg_est = avrg_edgelength_estimate(mesh);
  preprocess_min_max(opts.min, opts.max, opts.avrg, avrg_est, min, max);

  // make sure this is a tet mesh
  for(size_t i=0; i<mesh.etype.size(); i++)
    if(mesh.etype[i] != Tetra) {
      std::cerr << "Error: Mesh resampling implemented only for tetrahedral meshes!" << std::endl;
      exit(1);
    }

  // pre-process operations
  mt_vector<MT_USET<mt_int> > stags;
  preprocess_surf_tags(mesh, opts.tags, stags);

  bool dosmooth = opts.postsmth.size() ? atoi(opts.postsmth.c_str()) : true;

  // this flag controls the output of intermediate results
  bool output_iter = false;

  float smth_base = SMTH_BASE, smth_change = SMTH_CHANGE, edge_ang = EDGE_ANG;
  float iter_base = ITER_BASE;

  for(int oidx = 0; oidx < (int)stags.size(); oidx++)
  {
    // the set of tags that define the region we work on
    printf("Resampling op %d/%d: Processing tags: ", oidx+1, int(stags.size()));
    for(auto t : stags[oidx]) printf("%d ", int(t));
    printf("\n");

    mt_edge_manager edge_manager(mesh);
    mt_edge_splitter splitter(edge_manager);
    mt_edge_collapser collapser(edge_manager);

    // splitting run ========================================================================
    if(max > 0)
    {
      mt_vector<bool> boundary_node;
      if(fix_bnd) {
        if(unif_splt) {
          fprintf(stderr, "Warning: Boundary fixing and uniform refinement exclude each other!\n"
              "Turning off uniform refinement!\n");
          unif_splt = false;
        }

        mt_meshdata surfmesh;
        compute_surface(mesh, stags[oidx], surfmesh, false);

        // if fix_bnd is set, we insert the surface nodes of the geometric surface into
        // the boundary_node mask, which is tested in the splitting phase
        boundary_node.resize(mesh.n2e_cnt.size(), false);
        for(const mt_int & nidx : surfmesh.e2n_con) boundary_node[nidx] = true;
      }

      bool splt_update;

      if(unif_splt)
        splt_update = uniform_splitting(mesh, stags[oidx], edge_manager, max);
      else
        splt_update = non_uniform_splitting(mesh, stags[oidx], splitter, edge_manager, boundary_node, max);

      if(splt_update) {
        // finish splitting
        correct_insideOut(mesh);
        edge_manager.refresh_edge_data(true);
        compute_full_mesh_connectivity(mesh);
        if (dosmooth) {
          volumetric_smooth_from_tags(mesh, stags, iter_base, smth_base, edge_ang, RESAMPLE_QUAL_THR, false, true);
        }
      }
    }

    // collapsing run ========================================================================
    if(min > 0) {
      size_t iter = 0;
      size_t nsel = 0, nsel_old = 0, nsel_init = 0;
      char msg[1024];

      mt_real scorr = opts.surf_corr.size() > 0 ?
        atof(opts.surf_corr.c_str()) : CORR_DFLT;

      // compute geometric surface
      std::cout << "Computing surface and surface normals .." << std::endl;
      gettimeofday(&t1, NULL);
      mt_vector<bool> surf_nod;
      mt_vector<mt_real> surf_nrml;
      mt_vector<bool> fixed_nod,       // nodes that we have to collapse towards
                      boundary_node;   // nodes that cannot take part in collapsing

      // compute surface normals and fixed nods if fix_bnd == 1
      update_surface_data(mesh, stags[oidx], fix_bnd, surf_nod, boundary_node, surf_nrml);

      {
        mt_meshdata surfmesh, linemesh;
        unified_surface_from_tags(mesh, stags, surfmesh, NULL);
        compute_full_mesh_connectivity(surfmesh);
        compute_line_interfaces(mesh, surfmesh, 0.0, false, linemesh);
        fixed_nod.assign(mesh.n2e_cnt.size(), false);
        for(auto c : linemesh.e2n_con) fixed_nod[c] = true;
      }

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::set<edgeele, edge_len_asc> sel_edges;
      mt_vector<edgeele>              elem_edge_lengths;
      // in order to reduce access to the std::set we also use this boolean vector
      mt_vector<bool> sel_edges_map;

      bool collapse_update = false;

      // if we have do_conv == false, we want to limit the number of iterations and
      // the number of smoothing passes we want to allow to save time and not change
      // the mesh too much
      const size_t max_passes = 8;
      const size_t max_iter   = dosmooth ? max_passes * ITER_THRS : 10000;
      // we use this iteration variables in case we have do_conv == true. Then we
      // want to only stop after ITER_THRS iterations with zero improvement. the reasoning
      // is that we want to wait for a smoothing step to occur and then retry collapsing
      // before we call it a day.
      size_t max_zero_iter = dosmooth ? ITER_THRS : 1, zero_iter = 0;

      // we store here the ratio between selected edges and total edge count
      double diff_to_nsel;
      bool iter_control_continue = false;
      long int diff = 0;

      // main collapsing loop
      do {
        const MT_MAP<tuple<mt_int>, mt_int> & edges = edge_manager.edges();
        sel_edges.clear();
        sel_edges_map.assign(edges.size(), false);

        for(size_t eidx=0; eidx<mesh.e2n_cnt.size(); eidx++) {
          if(stags[oidx].count(mesh.etags[eidx])) {
            add_edges_lengths(mesh, eidx, elem_edge_lengths);

            for(edgeele & e : elem_edge_lengths) {
              bool edge_is_fixed = false;

              if(fix_bnd) {
                bool v1_fixed = boundary_node[e.v1],
                     v2_fixed = boundary_node[e.v2];
                edge_is_fixed = v1_fixed || v2_fixed;
              }

              // we select an edge if it is not going to a fixed boundary, if it is smaller
              // than the minimum allowed size and if it hasnt been selected yet
              if(!edge_is_fixed && e.length < min) {
                auto it = edges.find({e.v1, e.v2});
                e.idx = it->second;

                if(sel_edges_map[e.idx] == false) {
                  sel_edges.insert(e);
                  sel_edges_map[e.idx] = true;
                }
              }
            }
          }
        }

        nsel_old = nsel = sel_edges.size();
        iter++;
        diff_to_nsel = 0.0;

        if(nsel) {
          // store the number of initially selected edges
          if(nsel_init == 0) nsel_init = nsel;

          double nsel_to_total = double(nsel) / double(edges.size());
          sprintf(msg, "Iter %lu: Attempting to collapse %lu edges (%.4lf %%)", iter, nsel, nsel_to_total * 100.0);

          PROGRESS<size_t> prog(4, msg);
          gettimeofday(&t1, NULL);

          // call collapser to do collapsing work
          collapser(sel_edges, surf_nod, surf_nrml, fixed_nod, scorr, prog);
          nsel = sel_edges.size();
          diff = nsel_old - nsel;
          diff_to_nsel = double(diff) / double(nsel_init);

          gettimeofday(&t2, NULL);
          printf("Collapsed %lu (%.2f %%) in %.2f sec\n",
              diff, float(diff) / float(nsel_old) * 100.0f,
              (float)timediff_sec(t1, t2));

          // if we collapsed some edges we mark that work has been done
          if(diff > 0) collapse_update = true;

          // if post-smoothing is turned on, we do a smoothing step every ITER_THRS iterations
          // we reduce smoothing strength every time we smoothed to reduce volumetric
          // shrinking
          if(dosmooth && iter % ITER_THRS == 0) {
            compute_full_mesh_connectivity(mesh);
            volumetric_smooth_from_tags(mesh, stags, iter_base, smth_base, edge_ang, RESAMPLE_QUAL_THR, false, false);
            update_surface_data(mesh, stags[oidx], fix_bnd, surf_nod, boundary_node, surf_nrml);
            zero_iter = 0;

            if(do_conv == false) {
              smth_base *= smth_change;
              iter_base *= smth_change;
            }
          }

          if(output_iter) {
            correct_insideOut(mesh);
            mt_filename outfile(opts.outmsh_base, opts.ofmt);
            std::string filestr = outfile.base + ".iter" + std::to_string(iter);
            mt_filename file(filestr, outfile.format);
            std::cout << "Writing mesh: " << file.base << std::endl;
            write_mesh_selected(mesh, file.format, file.base);
          }
        }

        // in this flag we set all iteration control logic: if do_conv is true, we
        // iterate unter no more collapses could be made. if it is false, we stop earlier
        // to save iterations
        if(do_conv) {
          if(diff == 0) zero_iter++;
          iter_control_continue = nsel > 0 && zero_iter < max_zero_iter;
        }
        else
          iter_control_continue = nsel > 0 && diff_to_nsel > 1e-3 && iter < max_iter;
      }
      while(iter_control_continue);

      if(collapse_update) {
        // finish collapsing
        correct_insideOut(mesh);
        edge_manager.refresh_edge_data(true);
        compute_full_mesh_connectivity(mesh);
      }
    }
  }

  // do some final smoothing
  if(dosmooth)
    volumetric_smooth_from_tags(mesh, stags, iter_base, smth_base, edge_ang, RESAMPLE_QUAL_THR, false, true);

  // write final mesh
  mt_filename outmsh(opts.outmsh_base, opts.ofmt);
  std::cout << "Writing mesh: " << outmsh.base << std::endl;
  write_mesh_selected(mesh, outmsh.format, outmsh.base);
}


/// update the mesh with a collapsing iter
bool surf_collapse_update(mt_meshdata & mesh,
                          const MT_USET<mt_int> & tags,
                          mt_edge_manager & edge_manager,
                          mt_edge_collapser & collapser,
                          const mt_real min,
                          const mt_real scorr,
                          const bool fix_bnd,
                          const bool dosmooth)
{
  struct timeval t1, t2;
  bool did_update_mesh = false;
  size_t nsel = 0, nsel_old = 0;
  size_t iter = 0;
  char msg[1024];

  std::cout << "\n *** Collapsing step *** \n" << std::endl;

  mt_vector<mt_real> surf_nrml;
  compute_nodal_surface_normals(mesh, mesh.xyz, surf_nrml);

#if 0
  write_vector_ascii(surf_nrml, "normals.vec", 3);
#endif

  // we keep the edges connected to only one element as fixed, 
  // since they define the boundary of a non-closed surface
  mt_vector<bool> fixed_nod;
  fixed_nod.assign(mesh.n2e_cnt.size(), false);

  if(fix_bnd) {
    const mt_mapping<mt_int> & e2g = edge_manager.elem2edges();
    for(auto e : edge_manager.edges())
    {
      mt_int edge_idx = e.second;
      if(e2g.bwd_cnt[edge_idx] == 1) {
        fixed_nod[e.first.v1] = true;
        fixed_nod[e.first.v2] = true;
      }
    }
  }

  mt_vector<edgeele> elem_edge_lengths;

  // main collapsing loop
  do {
    std::set< edgeele, edge_len_asc> sel_edges;
    const MT_MAP<tuple<mt_int>, mt_int> & edges = edge_manager.edges();

    for(size_t eidx=0; eidx<mesh.e2n_cnt.size(); eidx++) {
      if(tags.count(mesh.etags[eidx])) {
        add_edges_lengths(mesh, eidx, elem_edge_lengths);

        for (edgeele & e : elem_edge_lengths) {
          if(e.length < min) {
            auto it = edges.find({e.v1, e.v2});
            e.idx = it->second;

            sel_edges.insert(e);
          }
        }
      }
    }

    nsel     = sel_edges.size();
    nsel_old = nsel;
    iter++;

    if(nsel) {
      double nsel_to_total = double(nsel) / double(edges.size()) * 100.0;
      sprintf(msg, "Iter %lu: Attempting to collapse %lu edges (%.4lf %%)", iter, nsel, nsel_to_total);
      PROGRESS<size_t> prog(4, msg);
      gettimeofday(&t1, NULL);

      collapser(sel_edges, surf_nrml, fixed_nod, scorr, prog);
      nsel = sel_edges.size();
      size_t diff = nsel_old - nsel;

      gettimeofday(&t2, NULL);
      printf("Collapsed %lu (%.2f %%) in %.2f sec\n",
          diff, float(diff) / float(nsel_old) * 100.0f,
          (float)timediff_sec(t1, t2));

      if(diff > 0) {
        did_update_mesh = true;
        compute_nodal_surface_normals(mesh, mesh.xyz, surf_nrml);
      }
    }
  } while(nsel > 0 && nsel_old > nsel);

  if(did_update_mesh) {
    // for checking if elements are insisde out we need the e2n_graph refreshed
    mesh.e2n_dsp.resize(mesh.e2n_cnt.size());
    bucket_sort_offset(mesh.e2n_cnt, mesh.e2n_dsp);
    correct_insideOut_tri(mesh, surf_nrml);

    // now we refresh edges and remap the points
    edge_manager.refresh_edge_data(true);

    // for simplicity we always recompte the connectivity as it might be needed in
    // splitting or smoothing runs
    compute_full_mesh_connectivity(mesh);

    if(dosmooth)
      surface_smooth(mesh, ITER_BASE, SMTH_BASE, EDGE_ANG, 0, false, false);
  }

  return did_update_mesh;
}

/// update the mesh with a splitting iter
bool surf_split_update(mt_meshdata & mesh,
                       const MT_USET<mt_int> & tags,
                       mt_edge_manager & edge_manager,
                       mt_edge_splitter & splitter,
                       const mt_real max,
                       const bool unif_splt,
                       const bool dosmooth)
{
  std::cout << "\n *** Splitting step *** \n" << std::endl;

  // to fix surface normal orientation issue, we will interpolate the surface normals
  // from the original mesh onto the refined mesh and use them to determine if the
  // element (triangle) definition is good or needs adjustment
  mt_meshdata origmesh = mesh;
  mt_vector<mt_real> surf_nrml, orig_surf_nrml;
  mt_vector<mt_point<mt_real>> surf_nrml_vec, orig_surf_nrml_vec;
  compute_nodal_surface_normals(origmesh, origmesh.xyz, orig_surf_nrml);

  // right now I dont see why we would need to fix the surface boundary during refinement, 
  // so the feature is off.
  mt_vector<bool> boundary_nodes;

  bool splt_update;

  if(unif_splt)
    splt_update = uniform_splitting(mesh, tags, edge_manager, max);
  else
    splt_update = non_uniform_splitting(mesh, tags, splitter, edge_manager, boundary_nodes, max);

  if(splt_update) {
    // finish splitting
    edge_manager.refresh_edge_data(true);
    compute_full_mesh_connectivity(mesh);

    mt_vector<mt_int> corr;
    mt_vector<mt_real> corr_dist;
    mt_manifold iman, oman;
    compute_correspondance(mesh, origmesh, corr, corr_dist);
    array_to_points(orig_surf_nrml, orig_surf_nrml_vec);
    nodal_interpolation(origmesh, mesh, iman, oman, corr, corr_dist, orig_surf_nrml_vec, surf_nrml_vec);
    points_to_array(surf_nrml_vec, surf_nrml);

    // TODO: smooth surface normal vecs prior to using them
    correct_insideOut_tri(mesh, surf_nrml);

    if(dosmooth)
      surface_smooth(mesh, ITER_BASE, SMTH_BASE, EDGE_ANG, 0, false, false);
  }

  return splt_update;
}


void resample_surf_mode(const resample_options & opts)
{
  struct timeval t1, t2;

  mt_meshdata mesh;
  mt_filename msh(opts.msh_base, opts.ifmt);
  std::cout << "Reading mesh: " << msh.base << std::endl;
  gettimeofday(&t1, NULL);
  read_mesh_selected(mesh, msh.format, msh.base);
  compute_full_mesh_connectivity(mesh);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  float min, max;
  float avrg_est = avrg_edgelength_estimate(mesh);
  preprocess_min_max(opts.min, opts.max, opts.avrg, avrg_est, min, max);

  // make sure this is a tri mesh
  for(size_t i=0; i<mesh.etype.size(); i++)
    if(mesh.etype[i] != Tri) {
      std::cerr << "Error: Surface resampling implemented only for triangular meshes!" << std::endl;
      exit(1);
    }

  // pre-process operations
  mt_vector<MT_USET<mt_int> > stags;
  preprocess_surf_tags(mesh, opts.tags, stags);

  bool dosmooth = opts.postsmth.size() ? atoi(opts.postsmth.c_str()) : true;

  for(int oidx = 0; oidx < (int)stags.size(); oidx++)
  {
    // the set of tags that define the region we work on
    printf("Resampling op %d/%d: Processing tags: ", oidx+1, int(stags.size()));
    for(auto t : stags[oidx]) printf("%d ", int(t));
    printf("\n");


    mt_edge_manager edge_manager(mesh);
    mt_edge_splitter splitter(edge_manager);
    mt_edge_collapser collapser(edge_manager);

    bool did_update_mesh;
    int  iter = 0, maxiter = 10;
    do {
      did_update_mesh = false;
      iter++;

      // collapsing run ========================================================================
      if(min > 0) {
        bool fix_bnd  = opts.fix.size() ? atoi(opts.fix.c_str()) : true;
        mt_real scorr = opts.surf_corr.size() > 0 ? atof(opts.surf_corr.c_str()) : CORR_DFLT;

        did_update_mesh = surf_collapse_update(mesh, stags[oidx], edge_manager, collapser, min, scorr, fix_bnd, dosmooth);
      }

      // splitting run ========================================================================
      if(max > 0) {
        bool unif_splt = opts.unif.size() ? atoi(opts.unif.c_str()) : 0;
        bool split_update = surf_split_update(mesh, stags[oidx], edge_manager, splitter, max, unif_splt, dosmooth);
        did_update_mesh = did_update_mesh || split_update;
      }
    } while(did_update_mesh && iter < maxiter);
  }

  if(dosmooth) {
    // do some final smoothing
    mt_meshdata & surfmesh = mesh;
    mt_meshdata linemesh;
    // extract surface and line manifolds
    MT_USET<mt_int> sm_vtx, ln_vtx;
    compute_line_interfaces(mesh, surfmesh, 0.25, true, linemesh);
    compute_full_mesh_connectivity(linemesh);

    sm_vtx.insert(mesh.e2n_con.begin(), mesh.e2n_con.end());
    ln_vtx.insert(linemesh.e2n_con.begin(), linemesh.e2n_con.end());
    // prepare data-struct for surface smoothing
    sm_vtx.sort(); ln_vtx.sort();
    for(auto n : ln_vtx) sm_vtx.erase(n);

    mt_vector<mt_int> sm_nod;
    mt_mask isMnfld(mesh.xyz.size() / 3);

    sm_nod.assign(sm_vtx.begin(), sm_vtx.end());
    isMnfld.insert(linemesh.e2n_con.begin(), linemesh.e2n_con.end());

    // smooth
    smooth_nodes(mesh, linemesh, isMnfld, sm_nod, 100, 0.25);
  }

#if 0
  {
    mt_vector<mt_real> surf_nrml;
    compute_nodal_surface_normals(mesh, mesh.xyz, surf_nrml);
    write_vector_ascii(surf_nrml, "normals.vec", 3);
  }
#endif

  // write final mesh
  mt_filename outmsh(opts.outmsh_base, opts.ofmt);
  std::cout << "Writing mesh: " << outmsh.base << std::endl;
  write_mesh_selected(mesh, outmsh.format, outmsh.base);
}



/// resample mode function
int resample_mode(int argc, char** argv)
{
  struct resample_options opts;

  int ret = resample_parse_options(argc, argv, opts);
  if (ret != 0) return 1;

  switch(opts.mode)
  {
    case RES_PURK:
      resample_purk_mode(opts);
      break;

    case RES_MESH:
      resample_mesh_mode(opts);
      break;

    case RES_SURF:
      resample_surf_mode(opts);
      break;

    default: break;
  }
  return 0;
}
