/**
* @file insert_mode.h
* @brief Meshtool insert mode.
* @author Aurel Neic
* @version
* @date 2017-02-13
*/

#include "mt_modes_base.h"

struct insert_options {
  enum op_type type;
  std::string msh_base;
  std::string submsh_base;
  std::string outmsh_base;
  std::string msh_dat_file;
  std::string submsh_dat_file;
  std::string out_dat_file;
  std::string idat;
  std::string ifmt;
  std::string ofmt;
  std::string trsp;
  std::string trsp_dist;
  std::string grad_thr;
  std::string con_thr;
  std::string oper;
  std::string mode;
};

static const std::string trsp_par="-trsp=";
static const std::string trsp_dist_par="-trsp_dist=";
static const std::string grad_par="-grad_thr=";
static const std::string con_par="-con_thr=";
#define GRAD_THR_DFLT 0.99
#define CON_THR_DFLT 1

void print_insert_submesh_help()
{
  fprintf(stderr, "insert submesh: a submesh is inserted back into a mesh and written to an output mesh\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the submesh to insert from\n", submesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh to insert into\n", mesh_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh output format. may be: %s\n", out_format_par.c_str(), output_formats.c_str());
  fprintf(stderr, "%s<path>\t (output) path to basename of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "\n");
  fprintf(stderr, "Note that the files defining the submesh must include a *.eidx and a *.nod file.\n"
                  "This files define how elements and nodes of the submesh map back into the original mesh.\n"
                  "The *.eidx and *.nod files are generated when using the \"extract mesh\" mode.\n");
  fprintf(stderr, "\n");
}
void print_insert_data_help()
{
  fprintf(stderr, "insert data: data defined on a submesh is inserted back into a mesh\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t (input) path to basename of the mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t\t (input) path to basename of the submesh\n", submesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to submesh data\n", submesh_data_par.c_str());
  fprintf(stderr, "%s<path>\t\t (output) path to output data\n", odat_par.c_str());
  fprintf(stderr, "%s<path>\t (optional) path to mesh data\n", mesh_data_par.c_str());
  fprintf(stderr, "%s<int>\t\t (optional) Data mode. 0 = nodal, 1 = element. Default is 0.\n", mode_par.c_str());
  fprintf(stderr, "\n");
  fprintf(stderr, "Note that the files defining the submesh must include a *.eidx and a *.nod file.\n"
                  "This files define how elements and nodes of the submesh map back into the original mesh.\n"
                  "The *.eidx and *.nod files are generated when using the \"extract mesh\" mode.\n");
  fprintf(stderr, "\n");
}
void print_insert_meshdata_help()
{
  fprintf(stderr, "insert meshdata: the fiber and tag data of a mesh is inserted into another mesh\n");
  fprintf(stderr, "                 based on vertex locations.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t (input) path to basename of the mesh to insert from\n",
          inp_mesh_par.c_str());
  fprintf(stderr, "%s<path>\t\t (input) path to basename of the mesh to insert into\n",
          mesh_par.c_str());
  fprintf(stderr, "%s<path>\t\t (input) Operation index: 0 = only tags, 1 = only fibers, 2 = both\n",
          oper_par.c_str());
  fprintf(stderr, "%s<format>\t\t (optional) mesh input format. may be: %s\n",
          inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "%s<format>\t\t (optional) mesh output format. may be: %s\n",
          out_format_par.c_str(), output_formats.c_str());
  fprintf(stderr, "%s<vec-file>\t (optional) element-based transport gradient\n", trsp_par.c_str());
  fprintf(stderr, "%s<float>\t (optional) transport distance threshold\n", trsp_dist_par.c_str());
  fprintf(stderr, "%s<float>\t (optional) gradient correlation threshold. Must be in (0, 1). Default is %g.\n",
          grad_par.c_str(), GRAD_THR_DFLT);
  fprintf(stderr, "%s<float>\t (optional) Connectivity threshold. Required number of nodes an element of mesh1\n"
                   "\t\t\t needs to share with an element of mesh2 in order for the data to be inserted. Default is %d.\n",
          con_par.c_str(), CON_THR_DFLT);
  fprintf(stderr, "%s<path>\t\t (output) path to basename of the output mesh\n",
          outmesh_par.c_str());
  fprintf(stderr, "\n");
}

void print_insert_tags_help()
{
  fprintf(stderr, "insert tags: insert an element data vector as new mesh element tags.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh.\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) .dat file the tags are inserted from.\n", idat_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to basename of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh input format.\n", inp_format_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh output format.\n", out_format_par.c_str());
  fprintf(stderr, "The supported input formats are:\n%s\n", input_formats.c_str());
  fprintf(stderr, "The supported output formats are:\n%s\n", output_formats.c_str());
  fprintf(stderr, "\n");
}




int insert_parse_options(int argc, char** argv, struct insert_options & opts)
{
  if(argc < 3) {
    //print_usage(argc, argv);
    std::cerr << "Error: Please choose one of the following insert modes: " << std::endl << std::endl;
    print_insert_submesh_help();
    print_insert_meshdata_help();
    print_insert_data_help();
    return 1;
  }

  std::string insert_type = argv[2];

  for(int i=3; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh_par, opts.msh_base);
    if(!match) match = parse_param(param, submesh_par, opts.submsh_base);
    if(!match) match = parse_param(param, inp_mesh_par, opts.submsh_base);
    if(!match) match = parse_param(param, outmesh_par, opts.outmsh_base);
    if(!match) match = parse_param(param, oper_par, opts.oper);
    if(!match) match = parse_param(param, mesh_data_par, opts.msh_dat_file);
    if(!match) match = parse_param(param, idat_par, opts.idat);
    if(!match) match = parse_param(param, submesh_data_par, opts.submsh_dat_file);
    if(!match) match = parse_param(param, odat_par, opts.out_dat_file);
    if(!match) match = parse_param(param, out_format_par, opts.ofmt);
    if(!match) match = parse_param(param, inp_format_par, opts.ifmt);
    if(!match) match = parse_param(param, trsp_par, opts.trsp);
    if(!match) match = parse_param(param, trsp_dist_par, opts.trsp_dist);
    if(!match) match = parse_param(param, grad_par, opts.grad_thr);
    if(!match) match = parse_param(param, con_par, opts.con_thr);
    if(!match) match = parse_param(param, mode_par, opts.mode);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 2;
    }
  }

  if(insert_type.compare("submesh") == 0)
  {
    opts.type = OP_MESH;

    if( !(opts.msh_base.size() > 0 && opts.submsh_base.size() > 0 && opts.outmsh_base.size() > 0) )
    {
      std::cerr << "Submesh insert error: Insufficient parameters provided." << std::endl;
      print_insert_submesh_help();
      return 4;
    }
  }
  else if(insert_type.compare("meshdata") == 0)
  {
    opts.type = OP_SURF;

    if(! (opts.msh_base.size() > 0 &&
          opts.submsh_base.size() > 0 &&
          opts.outmsh_base.size() > 0 &&
          opts.oper.size() > 0) )
    {
      std::cerr << "Meshdata insert error: Insufficient parameters provided." << std::endl;
      print_insert_meshdata_help();
      return 4;
    }
  }
  else if(insert_type.compare("data") == 0)
  {
    opts.type = OP_DATA;

    if( !(opts.msh_base.size() > 0 && opts.submsh_base.size() > 0 &&
          opts.submsh_dat_file.size() > 0 && opts.out_dat_file.size() > 0) )
    {
      std::cerr << "Data insert error: Insufficient parameters provided." << std::endl;
      print_insert_data_help();
      return 4;
    }
  }
  else if(insert_type.compare("tags") == 0)
  {
    opts.type = OP_TAGS;

    if( !(opts.msh_base.size() > 0 && opts.idat.size() > 0 && opts.outmsh_base.size() > 0) )
    {
      std::cerr << "Tags insert error: Insufficient parameters provided." << std::endl;
      print_insert_tags_help();
      return 4;
    }
  }
  else {
    print_usage(argv[0]);
    return 2;
  }

  return 0;
}


template<class T>
void insert_1dpn(const mt_vector<mt_int> & nod, mt_vector<T> & in, mt_vector<T> & out)
{
  #ifdef OPENMP
  #pragma omp parallel for schedule(dynamic, 256)
  #endif
  for(size_t i=0; i<nod.size(); i++)
    out[nod[i]] = in[i];
}
template<class T>
void insert_3dpn(const mt_vector<mt_int> & nod, mt_vector<T> & in, mt_vector<T> & out)
{
  #ifdef OPENMP
  #pragma omp parallel for schedule(dynamic, 128)
  #endif
  for(size_t i=0; i<nod.size(); i++) {
    out[nod[i]*3+0] = in[i*3+0];
    out[nod[i]*3+1] = in[i*3+1];
    out[nod[i]*3+2] = in[i*3+2];
  }
}



void insert_mode(int argc, char** argv)
{
  struct insert_options opts;
  int ret = insert_parse_options(argc, argv, opts);
  struct timeval t1, t2;

  if(ret != 0) return;

  struct mt_meshdata mesh;

  switch(opts.type) {
    case OP_MESH: {
      // insert mesh  ------------------------------------------------------------
      mt_filename mshfile   (opts.msh_base, opts.ifmt);
      mt_filename submshfile(opts.submsh_base, opts.ifmt);
      mt_filename outmshfile(opts.outmsh_base, opts.ofmt);

      std::cout << "Reading submesh: " << submshfile.base << std::endl;
      gettimeofday(&t1, NULL);

      mt_meshdata submsh;
      read_mesh_selected(submsh, submshfile.format, submshfile.base);
      mt_vector<mt_int> eidx(submsh.etags.size()), nod(submsh.xyz.size() / 3);
      binary_read(eidx, submshfile.base + EIDX_EXT);
      binary_read(nod,  submshfile.base + NOD_EXT);

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Reading mesh: " << mshfile.base << std::endl;
      gettimeofday(&t1, NULL);
      read_mesh_selected(mesh, mshfile.format, mshfile.base);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Inserting submesh ..." << std::endl;
      gettimeofday(&t1, NULL);
      insert_etags(mesh, submsh.etags, eidx);
      insert_fibers(mesh, submsh.lon, eidx);
      insert_points(mesh, submsh.xyz, nod);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Writing mesh: " << outmshfile.base << std::endl;
      gettimeofday(&t1, NULL);
      write_mesh_selected(mesh, outmshfile.format, outmshfile.base);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      break;
    }

    case OP_DATA: {
      bool nodal_data = true;
      if(opts.mode.size()) nodal_data = atoi(opts.mode.c_str()) != 1;

      // insert data  ------------------------------------------------------------
      mt_filename mshfile   (opts.msh_base, "");
      mt_filename submshfile(opts.submsh_base, "");

      std::string indices_file;
      if(nodal_data)
        indices_file = opts.submsh_base + NOD_EXT;
      else
        indices_file = opts.submsh_base + EIDX_EXT;

      size_t numidx = file_size(indices_file.c_str()) / sizeof(mt_int);
      if(numidx == 0) {
        fprintf(stderr, "%s error: Cannot determine number of indices. Aborting!\n", __func__);
        exit(EXIT_FAILURE);
      }

      // read mapping
      mt_vector<mt_int> indices(numidx);
      std::cout << "Reading " << indices_file << std::endl;
      binary_read(indices, indices_file);

      // deal with missing mesh data
      bool msh_data_exists = opts.msh_dat_file.size() > 0 && file_exists(opts.msh_dat_file);
      mt_int num_mesh_idx = 0;
      if(msh_data_exists == false) {
        read_mesh_selected(mesh, mshfile.format, mshfile.base, CRP_READ_ELEM | CRP_READ_PTS);
        num_mesh_idx = nodal_data ? mesh.xyz.size() / 3 : mesh.e2n_cnt.size();
      }

      // read and insert data
      short data_idx = -1, dpn = 1;
      mt_vector<mt_real> idat, odat;
      igb_header igb_msh, igb_submsh;

      setup_data_format(opts.submsh_dat_file, opts.msh_dat_file, data_idx, igb_submsh, igb_msh);
      if(data_idx == 2 || data_idx == 3) dpn = 3;

      gettimeofday(&t1, NULL);
      std::cout << "Inserting data .. " << std::endl;
      switch(data_idx) {
        case 0:
        case 2:

          std::cout << "Reading input data: " << opts.submsh_dat_file << std::endl;
          read_vector_ascii(idat, opts.submsh_dat_file, true);
          if(msh_data_exists) {
            std::cout << "Reading input data: " << opts.msh_dat_file << std::endl;
            read_vector_ascii(odat, opts.msh_dat_file, true);
          }
          else {
            odat.assign(num_mesh_idx*dpn, 0.0);
          }

          if(data_idx == 0) {
            insert_1dpn(indices, idat, odat);
            write_vector_ascii(odat, opts.out_dat_file, 1);
          }
          else {
            insert_3dpn(indices, idat, odat);
            write_vector_ascii(odat, opts.out_dat_file, 3);
          }

          std::cout << "Writing output data: " << opts.out_dat_file << std::endl;
          break;

        case 1:
        case 3:
        {
          if(igb_msh.v_t != igb_submsh.v_t) {
            fprintf(stderr, "Number of time-slices in mesh and submesh do not match. Aborting!\n");
            exit(1);
          }

          // set new spatial size and write output igb header
          igb_header igb_out = igb_msh;
          igb_out.filename = opts.out_dat_file;
          write_igb_header(igb_out);

          std::vector<mt_real> rbuff;
          printf("Inserting igb time-slices %s into %s, writing to %s. \n", opts.submsh_dat_file.c_str(),
                 opts.msh_dat_file.c_str(), opts.out_dat_file.c_str());

          for(int t=0; t<igb_submsh.v_t; t++) {
            printf("\rcurrent time-slice %d / %d .. ", t+1, int(igb_submsh.v_t));
            fflush(stdout);

            read_igb_slice(rbuff, igb_submsh);
            idat.assign(rbuff.begin(), rbuff.end());
            read_igb_slice(rbuff, igb_msh);
            odat.assign(rbuff.begin(), rbuff.end());

            if(msh_data_exists) {
              read_igb_slice(rbuff, igb_msh);
              odat.assign(rbuff.begin(), rbuff.end());
            }
            else {
              odat.assign(num_mesh_idx*dpn, 0.0);
            }

            if(data_idx == 1) insert_1dpn(indices, idat, odat);
            else              insert_3dpn(indices, idat, odat);

            rbuff.assign(odat.begin(), odat.end());
            write_igb_slice(rbuff, igb_out);
          }
          printf("\n");

          fclose(igb_msh.fileptr);
          fclose(igb_submsh.fileptr);
          fclose(igb_out.fileptr);
          break;
        }
        default: break;
      }
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      break;
    }

    /*  The module maps element data defined on a manifold (line, surface) of the
     *  mesh onto mesh elements.
     *
     *  The main idea is to create a manifoldEle-to-meshEle map. With this map,
     *  one can select the manifold element sharing the most nodes with each
     *  mesh element.
     *
     */
    case OP_SURF:
    {
      struct mt_meshdata manifold;
      size_t manifold_nnodes;

      mt_filename mshfile   (opts.msh_base, opts.ifmt);
      mt_filename submshfile(opts.submsh_base, opts.ifmt);
      mt_filename outmshfile(opts.outmsh_base, opts.ofmt);

      std::cout << "Reading mesh: " << mshfile.base << std::endl;
      gettimeofday(&t1, NULL);
      read_mesh_selected(mesh, mshfile.format, mshfile.base);
      reindex_nodes(mesh);
      compute_full_mesh_connectivity(mesh, mshfile.base);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      short numfib = mesh.e2n_cnt.size()*6 == mesh.lon.size() ? 2 : 1;

      std::cout << "Reading input mesh: " << submshfile.base << std::endl;
      gettimeofday(&t1, NULL);
      read_mesh_selected(manifold, submshfile.format, submshfile.base);
      reindex_nodes(manifold);
      manifold_nnodes = manifold.xyz.size() / 3;
      compute_full_mesh_connectivity(manifold);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      // consistency check
      short manifold_numfib = manifold.e2n_cnt.size()*6 == manifold.lon.size() ? 2 : 1;
      if (numfib != manifold_numfib) {
        std::cout << "Warning: imsh and msh must have the same fiber/sheet setup!\n";
        std::cout << "         Code is likely to crash...\n";
      }

      mt_vector<mt_int> corr;
      mt_vector<mt_real> corr_dist;
      std::cout << "Computing vertex correspondance .. " << std::endl;
      gettimeofday(&t1, NULL);
      compute_correspondance(manifold, mesh, corr, corr_dist);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Computing mapping between mesh and input-mesh.. " << std::endl;
      gettimeofday(&t1, NULL);

      MT_MAP<mt_int, mt_int> old2new, mesh2manifold;
      for(size_t i = 0; i<corr.size(); i++) mesh2manifold[corr[i]] = i;

      // select the initial elements we use in data insertion
      mt_meshdata res_mesh = mesh;
      mt_vector<mt_int> res_eidx;
      {
        MT_USET<mt_int> sel_elem, nodeSet;
        mt_vector<bool> keep(res_mesh.e2n_cnt.size(), false);

        nodeSet.insert(corr.begin(), corr.end());
        nodeSet_to_elemSet(mesh, nodeSet, sel_elem);

        // restrict mesh to selected elements
        for(const mt_int & e : sel_elem) keep[e] = true;
        restrict_meshdata(keep, res_mesh, res_eidx);
      }

      // map restricted mesh into a nodal indexing consistent with the manifold
      size_t vtxidx = manifold_nnodes;
      for(size_t i=0; i<res_mesh.e2n_con.size(); i++)
      {
        mt_int c = res_mesh.e2n_con[i], nc = -1;
        auto mc = mesh2manifold.find(c);
        if(mc != mesh2manifold.end()) {
          // set new connectivity from manifold correspondance
          nc = mc->second;
          old2new[c] = nc;
        }
        else {
          // set new connectivity from extended indexing
          mc = old2new.find(c);
          if(mc != old2new.end()) nc = mc->second;
          else {
            nc = vtxidx++;
            old2new[c] = nc;
          }
        }
        res_mesh.e2n_con[i] = nc;
      }
      compute_full_mesh_connectivity(res_mesh);

      mt_mapping<mt_int> mesh2manifld;
      mt_vector<mt_int> mesh_mult_cnt;
      /* We compute the mesh2manifld map by a matrix-matrix product.
       * The problem is that the manifold nodes dont span the whole
       * mesh node range. We must extend the manifold node range to
       * make the product work.
       */
      {
        mt_vector<mt_int> manifld_mult_cnt;
        size_t manifold_old_nsize = manifold.n2e_cnt.size();
        mt_vector<mt_int> aele(res_mesh.e2n_con.size(), 1);
        mt_vector<mt_int> bele(manifold.n2e_con.size(), 1);

        manifold.n2e_cnt.resize(vtxidx, 0);
        mat_mult_mat_crs(res_mesh.e2n_cnt, res_mesh.e2n_con, aele,
                         manifold.n2e_cnt, manifold.n2e_dsp, manifold.n2e_con, bele,
                         mesh2manifld.fwd_cnt, mesh2manifld.fwd_con, mesh_mult_cnt);
        manifold.n2e_cnt.resize(manifold_old_nsize);

        transpose_connectivity(mesh2manifld.fwd_cnt, mesh2manifld.fwd_con, mesh_mult_cnt,
                               mesh2manifld.bwd_cnt, mesh2manifld.bwd_con, manifld_mult_cnt);
        mesh2manifld.setup_dsp();
      }
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      // map the mesh data =====================================================
      std::cout << "Mapping element data .. " << std::endl;
      gettimeofday(&t1, NULL);

      mt_int mult_thr = opts.con_thr.size() > 0 ? atoi(opts.con_thr.c_str()) : CON_THR_DFLT;

      for(size_t eidx = 0; eidx < res_mesh.e2n_cnt.size(); eidx++)
      {
        mt_int maxidx = -1, maxcnt = 0;
        mt_int start = mesh2manifld.fwd_dsp[eidx], stop = start + mesh2manifld.fwd_cnt[eidx];
        // sweep over connected elements and save element with maximum connected nodes
        for(mt_int i=start; i<stop; i++)
        {
          mt_int midx = mesh2manifld.fwd_con[i];
          mt_int cc   = mesh_mult_cnt[i];
          if( cc >= mult_thr && cc > maxcnt) {
            maxcnt = mesh_mult_cnt[i];
            maxidx = midx;
          }
        }
        // set element data for saved element
        if(maxidx > -1) {
          // tags
          res_mesh.etags[eidx] = manifold.etags[maxidx];
          // fibers
          for(short i=0; i < numfib*3; i++)
            res_mesh.lon[eidx*3*numfib + i] = manifold.lon[maxidx*3*numfib + i];
        }
      }

      mesh2manifld.resize(0);
      short insert_op = atoi(opts.oper.c_str());

      // insert tags =====================================================
      if((insert_op == 0 || insert_op == 2))
      {
        if( (opts.trsp.size() > 0) &&  (opts.trsp_dist.size() > 0) )
        {
          // read in transport gradient
          mt_vector<mt_real> inp_grad;
          read_vector_ascii(inp_grad, opts.trsp);
          assert(inp_grad.size() / 3 == mesh.e2n_cnt.size());
          mt_vector<mt_point<mt_real> > grad(mesh.e2n_cnt.size());
          for(size_t i=0; i<grad.size(); i++) {
            grad[i].x = inp_grad[i*3+0];
            grad[i].y = inp_grad[i*3+1];
            grad[i].z = inp_grad[i*3+2];
            if(grad[i].length2() == 0)
              std::cerr << "Error: Zero gradient at " << i << std::endl;
            grad[i].normalize();
          }
          // read in transport distance
          mt_real trsp_dist = atof(opts.trsp_dist.c_str());
          // read gradient correlation threshold
          mt_real grad_thr  = opts.grad_thr.size() > 0 ? atof(opts.grad_thr.c_str()) :
                              GRAD_THR_DFLT;

          // set up e2e map
          mt_vector<mt_int> e2e_cnt, e2e_dsp, e2e_con;
          multiply_connectivities(mesh.e2n_cnt, mesh.e2n_con, mesh.n2e_cnt, mesh.n2e_con,
                                  e2e_cnt, e2e_con);
          e2e_dsp.resize(e2e_cnt.size()); bucket_sort_offset(e2e_cnt, e2e_dsp);

          // compute element centerpoints as average of node coords
          mt_vector<mt_point<mt_real> > vtx_nod, vtx_ele;
          array_to_points(mesh.xyz, vtx_nod);
          nodeData_to_elemData(mesh, vtx_nod, vtx_ele);

          // initialize element data transport
          mt_data_transport<mt_int, mt_real> transport(e2e_cnt, e2e_dsp, e2e_con, vtx_ele);

          std::set<mixed_tuple<mt_int,mt_int> > init;
          for(size_t i=0; i<res_eidx.size(); i++)
            init.insert( { res_eidx[i], res_mesh.etags[i] } );

          // apply data transport
          transport(grad, grad_thr, trsp_dist, mesh.etags, init);
        }
        else {
          insert_etags(mesh, res_mesh.etags, res_eidx);
        }
      }

      // insert fibers =========================================================
      if(insert_op == 1 || insert_op == 2)
        insert_fibers(mesh, res_mesh.lon, res_eidx);

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Writing mesh: " << outmshfile.base << std::endl;
      gettimeofday(&t1, NULL);
      write_mesh_selected(mesh, outmshfile.format, outmshfile.base);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      break;
    }

    case OP_TAGS: {
      mt_filename mshfile   (opts.msh_base, opts.ifmt);
      mt_filename outmshfile(opts.outmsh_base, opts.ofmt);

      std::cout << "Reading mesh: " << mshfile.base << std::endl;
      gettimeofday(&t1, NULL);
      read_mesh_selected(mesh, mshfile.format, mshfile.base);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      mt_vector<mt_int> newtags;
      read_vector_ascii(newtags, opts.idat, false);
      check_same_size(mesh.etags, newtags, __func__);
      mesh.etags = newtags;

      std::cout << "Writing mesh: " << outmshfile.base << std::endl;
      gettimeofday(&t1, NULL);
      write_mesh_selected(mesh, outmshfile.format, outmshfile.base);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      break;
    }

    default: break;
  }
}
