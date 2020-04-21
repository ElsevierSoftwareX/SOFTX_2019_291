/**
* @file interpolate_mode.cpp
* @brief Interpolate data between meshes.
* @author Aurel Neic
* @version
* @date 2017-10-27
*/

#include "mt_modes_base.h"
#include "dense_mat.hpp"

enum itp_type {ITP_NODES, ITP_ELEM, ITP_E2N, ITP_N2E, ITP_CLOUD};

struct interpolate_options {
  mt_filename imsh;
  mt_filename omsh;
  std::string pts;
  std::string idat;
  std::string odat;
  std::string thr;
  std::string pts_dynpt;
  std::string interp_mode;
  enum itp_type type;
};

static const std::string pts_par          = "-pts=";
static const std::string dynpts_par       = "-dynpts=";

void print_interpolate_nodedata_help()
{
  fprintf(stderr, "interpolate nodedata: interpolate nodal data from one mesh onto another\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh we interpolate to\n", out_mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh we interpolate from\n", inp_mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to input data.\n", idat_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to output data.\n", odat_par.c_str());
  fprintf(stderr, "%s<path>\t (optional) Dynpts describing point cloud movement.\n", dynpts_par.c_str());
  fprintf(stderr, "\n");
}
void print_interpolate_elemdata_help()
{
  fprintf(stderr, "interpolate elemdata: interpolate element data from one mesh onto another\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh we interpolate to\n", out_mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh we interpolate from\n", inp_mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to input data.\n", idat_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to output data.\n", odat_par.c_str());
  fprintf(stderr, "%s<path>\t (optional) Dynpts describing point cloud movement.\n", dynpts_par.c_str());
  fprintf(stderr, "\n");
}
void print_interpolate_clouddata_help()
{
  fprintf(stderr, "interpolate clouddata: interpolate data from a pointcloud onto a mesh using radial basis function interpolation\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh we interpolate to\n", out_mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to the coordinates of the point cloud\n", pts_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to input data.\n", idat_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to output data.\n", odat_par.c_str());
  fprintf(stderr, "%s<int>\t (optional) Choose between localized Shepard (=0), global Shepard (=1), and RBF interpolation (=2). Default is 2.\n", mode_par.c_str());
  fprintf(stderr, "%s<path>\t (optional) Dynpts describing point cloud movement.\n", dynpts_par.c_str());
  fprintf(stderr, "\n");
}
void print_interpolate_elem2node_help()
{
  fprintf(stderr, "interpolate elem2node: interpolate data from elements onto nodes\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh we interpolate to\n", out_mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to input data.\n", idat_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to output data.\n", odat_par.c_str());
  fprintf(stderr, "\n");
}
void print_interpolate_node2elem_help()
{
  fprintf(stderr, "interpolate node2elem: interpolate data from nodes onto elements\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh we interpolate to\n", out_mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to input data.\n", idat_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to output data.\n", odat_par.c_str());
  fprintf(stderr, "\n");
}


int interpolate_parse_options(int argc, char** argv, struct interpolate_options & opts)
{
  if(argc < 3) {
    std::cerr << "Please choose one of the following interpolate modes: " << std::endl << std::endl;
    print_interpolate_nodedata_help();
    print_interpolate_elemdata_help();
    print_interpolate_clouddata_help();
    print_interpolate_elem2node_help();
    print_interpolate_node2elem_help();
    return 1;
  }

  std::string interp_type = argv[2];

  std::string imsh_base;
  std::string omsh_base;

  for(int i=3; i<argc; i++) {
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, inp_mesh_par, imsh_base);
    if(!match) match = parse_param(param, out_mesh_par, omsh_base);
    if(!match) match = parse_param(param, idat_par, opts.idat, chk_fexists);
    if(!match) match = parse_param(param, odat_par, opts.odat);
    if(!match) match = parse_param(param, thr_par, opts.thr);
    if(!match) match = parse_param(param, pts_par, opts.pts);
    if(!match) match = parse_param(param, dynpts_par, opts.pts_dynpt);
    if(!match) match = parse_param(param, mode_par, opts.interp_mode);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 2;
    }
  }

  if(interp_type.compare("nodedata") == 0)
  {
    opts.imsh.assign(imsh_base, "");
    opts.omsh.assign(omsh_base, "");

    opts.type = ITP_NODES;

    if(! (opts.imsh.isSet() && opts.omsh.isSet() &&
          opts.idat.size() && opts.odat.size()) ) {
      std::cerr << "Interpolate nodedata error: Insufficient parameters provided." << std::endl;
      print_interpolate_nodedata_help();
      return 3;
    }
  }
  else if(interp_type.compare("elemdata") == 0)
  {
    opts.imsh.assign(imsh_base, "");
    opts.omsh.assign(omsh_base, "");

    opts.type = ITP_ELEM;

    if(! (opts.imsh.isSet() && opts.omsh.isSet() &&
          opts.idat.size() && opts.odat.size()) ) {
      std::cerr << "Interpolate elemdata error: Insufficient parameters provided." << std::endl;
      print_interpolate_elemdata_help();
      return 3;
    }
  }
  else if(interp_type.compare("clouddata") == 0)
  {
    opts.omsh.assign(omsh_base, "");

    opts.type = ITP_CLOUD;

    if(! (opts.omsh.isSet() && opts.pts.size() &&
          opts.idat.size() && opts.odat.size()) ) {
      std::cerr << "Interpolate clouddata error: Insufficient parameters provided." << std::endl;
      print_interpolate_clouddata_help();
      return 3;
    }
  }
  else if(interp_type.compare("elem2node") == 0)
  {
    opts.omsh.assign(omsh_base, "");

    opts.type = ITP_E2N;

    if(! (opts.omsh.isSet() && opts.idat.size() && opts.odat.size()) ) {
      std::cerr << "Interpolate elem2node error: Insufficient parameters provided." << std::endl;
      print_interpolate_elem2node_help();
      return 3;
    }
  }
  else if(interp_type.compare("node2elem") == 0)
  {
    opts.omsh.assign(omsh_base, "");

    opts.type = ITP_N2E;

    if(! (opts.omsh.isSet() && opts.idat.size() && opts.odat.size()) ) {
      std::cerr << "Interpolate node2elem error: Insufficient parameters provided." << std::endl;
      print_interpolate_node2elem_help();
      return 3;
    }
  }
  else {
    print_usage(argv[0]);
    return 2;
  }

  return 0;
}


void interpolate_mode(int argc, char** argv)
{
  struct interpolate_options opts;
  int ret = interpolate_parse_options(argc, argv, opts);
  if(ret != 0) return;

  struct timeval t1, t2;
  struct mt_meshdata imesh, omesh;
  mt_vector<mt_real> ipts;

  // data_idx may be:
  // 0 : scalar dat file
  // 1 : scalar igb file
  // 2 : vec file
  // 3 : vec3 igb file
  // 4 : vec9 igb file
  short data_idx = -1;
  mt_vector<mt_real>            idat,     odat;
  mt_vector<mt_point<mt_real> > idat_vec, odat_vec;
  mt_vector<dmat<mt_real> >     idat_ten, odat_ten;
  igb_header igb, igb_out;

  setup_data_format(opts.idat, opts.odat, data_idx, igb, igb_out);

  // first read mesh
  gettimeofday(&t1, NULL);
  std::cout << "Reading output mesh: " << opts.omsh.base << std::endl;
  read_mesh_selected(omesh, opts.omsh.format, opts.omsh.base, CRP_READ_ELEM | CRP_READ_PTS);
  compute_full_mesh_connectivity(omesh, opts.omsh.base);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  switch(opts.type) {
    case ITP_NODES:
    {
      std::cout << "Reading input mesh: " << opts.imsh.base << std::endl;
      gettimeofday(&t1, NULL);
      read_mesh_selected(imesh, opts.imsh.format, opts.imsh.base, CRP_READ_ELEM | CRP_READ_PTS);
      compute_full_mesh_connectivity(imesh, opts.imsh.base);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Compute correspondance .." << std::endl;
      gettimeofday(&t1, NULL);
      mt_vector<mt_int> corr;
      mt_vector<mt_real> corr_dist;
      compute_correspondance(omesh, imesh, corr, corr_dist);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      mt_manifold imnfld, omnfld;
      imnfld.on_mnfld.resize(imesh.xyz.size() / 3),
      omnfld.on_mnfld.resize(imesh.xyz.size() / 3);

      if(imesh.etype[0] == Tetra) {
        compute_surface(imesh.etype, imesh.e2n_cnt, imesh.e2n_con, imnfld.mesh, true);
        compute_surface(omesh.etype, omesh.e2n_cnt, omesh.e2n_con, omnfld.mesh, true);
        compute_full_mesh_connectivity(imnfld.mesh);
        compute_full_mesh_connectivity(omnfld.mesh);
        imnfld.mesh.xyz = imesh.xyz; omnfld.mesh.xyz = omesh.xyz;
        imnfld.set_mask(), omnfld.set_mask();
      }

      std::cout << "Compute interpolation .." << std::endl;
      gettimeofday(&t1, NULL);
      switch(data_idx) {
        case 0:
        case 2:
          read_vector_ascii(idat, opts.idat, true);

          if(data_idx == 0) {
            check_size(idat, imesh.xyz.size() / 3, "node interpolation");
            nodal_interpolation(imesh, omesh, imnfld, omnfld, corr, corr_dist, idat, odat);
            write_vector_ascii(odat, opts.odat, 1);
          }
          else {
            check_size(idat, imesh.xyz.size(), "node interpolation");
            array_to_points(idat, idat_vec);
            nodal_interpolation(imesh, omesh, imnfld, omnfld, corr, corr_dist, idat_vec, odat_vec);
            write_vector_ascii(odat_vec, opts.odat);
          }
          break;

        case 1:
        case 3:
        case 4:
        {
          const bool have_dynpts = opts.pts_dynpt.size() > 0;
          igb_header * igb_pts_dynpt = have_dynpts ? new igb_header
            : nullptr;
          if(have_dynpts) {
            init_igb_header(opts.pts_dynpt, *igb_pts_dynpt);
            read_igb_header(*igb_pts_dynpt);
          }
          check_size(size_t(igb.v_x), imesh.xyz.size() / 3, "element interpolation");

          // set new spatial size and write output igb header
          igb_out.v_x = corr.size();
          write_igb_header(igb_out);

          std::vector<mt_real> rbuff;
          printf("Processing igb time-slices: %s to %s: \n", opts.idat.c_str(), opts.odat.c_str());

          for(int t=0; t<igb.v_t; t++) {
            printf("\rcurrent time-slice %d / %d .. ", t+1, int(igb.v_t));
            fflush(stdout);

            read_igb_slice(rbuff, igb);
            idat.assign(rbuff.begin(), rbuff.end());
            if(have_dynpts) {
              read_igb_slice(rbuff, *igb_pts_dynpt);
              imesh.xyz.assign(rbuff.begin(), rbuff.end());
              compute_correspondance(omesh, imesh, corr, corr_dist);
              if(imesh.etype[0] == Tetra)
                imnfld.mesh.xyz = imesh.xyz;
            }
            if(data_idx == 1) {
              nodal_interpolation(imesh, omesh, imnfld, omnfld, corr, corr_dist, idat, odat);
            }
            else if(data_idx == 3) {
              array_to_points(idat, idat_vec);
              nodal_interpolation(imesh, omesh, imnfld, omnfld, corr, corr_dist, idat_vec, odat_vec);
              points_to_array(odat_vec, odat);
            }
            else {
              array_to_tensors(idat, idat_ten);
              nodal_interpolation(imesh, omesh, imnfld, omnfld, corr, corr_dist, idat_ten, odat_ten);
              tensors_to_array(odat_ten, odat);
            }

            rbuff.assign(odat.begin(), odat.end());
            write_igb_slice(rbuff, igb_out);
          }
          printf("\n");

          fclose(igb.fileptr);
          fclose(igb_out.fileptr);
          if(have_dynpts) {
            fclose(igb_pts_dynpt->fileptr);
            delete igb_pts_dynpt;
          }
          break;
        }
        default: break;
      }

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      break;
    }

    case ITP_ELEM:
    {
      std::cout << "Reading input mesh: " << opts.imsh.base << std::endl;
      gettimeofday(&t1, NULL);
      read_mesh_selected(imesh, opts.imsh.format, opts.imsh.base, CRP_READ_ELEM | CRP_READ_PTS);
      compute_full_mesh_connectivity(imesh, opts.imsh.base);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Setting up kdtree for input mesh vertices .." << std::endl;
      gettimeofday(&t1, NULL);

      kdtree vert_tree(10);
      vert_tree.build_vertex_tree(imesh.xyz);

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Compute interpolation .." << std::endl;
      gettimeofday(&t1, NULL);
      switch(data_idx) {
        case 0:
        case 2:
          read_vector_ascii(idat, opts.idat, true);

          if(data_idx == 0) {
            check_size(idat.size(), imesh.e2n_cnt.size(), "element interpolation");
            elementwise_interpolation(imesh, omesh, vert_tree, idat, odat);
            write_vector_ascii(odat, opts.odat, 1);
          }
          else {
            check_size(idat.size(), imesh.e2n_cnt.size()*3, "element interpolation");
            array_to_points(idat, idat_vec);
            elementwise_interpolation(imesh, omesh, vert_tree, idat_vec, odat_vec);
            write_vector_ascii(odat_vec, opts.odat);
          }
          break;

        case 1:
        case 3:
        case 4:
        {
          const bool have_dynpts = opts.pts_dynpt.size() > 0;
          igb_header * igb_pts_dynpt = have_dynpts ? new igb_header
            : nullptr;
          if(have_dynpts) {
            init_igb_header(opts.pts_dynpt, *igb_pts_dynpt);
            read_igb_header(*igb_pts_dynpt);
          }
          kdtree * act_tree = nullptr;
          check_size(size_t(igb.v_x), imesh.e2n_cnt.size(), "element interpolation");

          // set new spatial size and write output igb header
          igb_out.v_x = omesh.e2n_cnt.size();
          write_igb_header(igb_out);

          std::vector<mt_real> rbuff;
          printf("Processing igb time-slices: %s to %s: \n", opts.idat.c_str(), opts.odat.c_str());

          for(int t=0; t<igb.v_t; t++) {
            printf("\rcurrent time-slice %d / %d .. ", t+1, int(igb.v_t));
            fflush(stdout);

            read_igb_slice(rbuff, igb);
            idat.assign(rbuff.begin(), rbuff.end());
            if(have_dynpts) {
              read_igb_slice(rbuff, *igb_pts_dynpt);
              imesh.xyz.assign(rbuff.begin(), rbuff.end());
              act_tree = new kdtree(10);
              act_tree->build_vertex_tree(imesh.xyz);
            } else {
              act_tree = &vert_tree;
            }
            if(data_idx == 1) {
              elementwise_interpolation(imesh, omesh, *act_tree, idat, odat);
            }
            else if(data_idx == 3) {
              array_to_points(idat, idat_vec);
              elementwise_interpolation(imesh, omesh, *act_tree, idat_vec, odat_vec);
              points_to_array(odat_vec, odat);
            }
            else if(data_idx == 4) {
              array_to_tensors(idat, idat_ten);
              elementwise_interpolation(imesh, omesh, *act_tree, idat_ten, odat_ten);
              tensors_to_array(odat_ten, odat);
            }

            rbuff.assign(odat.begin(), odat.end());
            write_igb_slice(rbuff, igb_out);
            if(have_dynpts)
              delete act_tree;
          }
          printf("\n");

          fclose(igb.fileptr);
          fclose(igb_out.fileptr);
          if(have_dynpts) {
            fclose(igb_pts_dynpt->fileptr);
            delete igb_pts_dynpt;
          }
          break;
        }
        default: break;
      }

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      break;
    }

    case ITP_CLOUD:
    {
      bool shepard = false;
      bool shepard_global = false;
      if(opts.interp_mode.size() > 0)
      {
        shepard = (atoi(opts.interp_mode.c_str()) < 2);
        shepard_global = (atoi(opts.interp_mode.c_str()) == 1);
      }

      std::cout << "Reading input points: " << opts.pts << std::endl;
      readPoints(ipts, opts.pts);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Compute interpolation .." << std::endl;
      gettimeofday(&t1, NULL);
      switch(data_idx) {
        case 0:
        case 2:
        {
          read_vector_ascii(idat, opts.idat, true);
          //Setup pointcloud interpolator
          PointCloudInterpolator cloud_interp(ipts);
          const int dpn = data_idx == 0 ? 1 : 3;
          const int num_neigbours = 7;
          cloud_interp.setup_roi(num_neigbours);
          if(shepard)
            cloud_interp.interpolate_shepard(omesh.xyz, idat, dpn, shepard_global, odat);
          else
          {
            cloud_interp.setup_rbf_coeffs(idat, dpn);
            cloud_interp.interpolate_data(omesh.xyz, odat);
          }
          write_vector_ascii(odat, opts.odat, dpn);
          break;
        }

        case 1:
        case 3:
        case 4:
        {
          // set new spatial size and write output igb header
          const bool have_dynpts = opts.pts_dynpt.size() > 0;
          igb_header igb_pts_dynpt;
          if(have_dynpts) {
            init_igb_header(opts.pts_dynpt, igb_pts_dynpt);
            read_igb_header(igb_pts_dynpt);
          }

          igb_out.v_x = omesh.xyz.size() / 3;
          write_igb_header(igb_out);

          int dpn = 0;
          const int num_neigbours = 7;
          PointCloudInterpolator cloud_interp;
          if(!have_dynpts) {
            cloud_interp.set_pts(ipts);
            cloud_interp.setup_roi(num_neigbours);
          }

          if (data_idx == 1)       dpn = 1;
          else if (data_idx == 3)  dpn = 3;
          else                     dpn = 9;

          std::vector<std::vector<mt_real> > rbuff;
          printf("Processing igb time-slices: %s to %s: \n", opts.idat.c_str(), opts.odat.c_str());

          for(int t=0; t<igb.v_t; t++) {
            printf("\rcurrent time-slice %d / %d .. ", t+1, int(igb.v_t));
            fflush(stdout);
            read_igb_block(rbuff, 1, igb);
            idat.assign(rbuff[0].begin(), rbuff[0].end());
            // apply deformation on input point set
            if(have_dynpts) {
              read_igb_block(rbuff, 1, igb_pts_dynpt);
              ipts.assign(rbuff[0].begin(), rbuff[0].end());
              cloud_interp.set_pts(ipts);
              cloud_interp.setup_roi(num_neigbours);
              if(!shepard)
                cloud_interp.setup_rbf_coeffs(idat, dpn);
            } else {
              if(!shepard && t == 0)
                cloud_interp.setup_rbf_coeffs(idat, dpn);
              else if(!shepard && t > 0)
                cloud_interp.recompute_rbf_coeffs(idat);
            }

            if(shepard)
              cloud_interp.interpolate_shepard(omesh.xyz, idat, dpn, shepard_global, odat);
            else
              cloud_interp.interpolate_data(omesh.xyz, odat);

            rbuff[0].assign(odat.begin(), odat.end());
            write_igb_block(rbuff, igb_out);
          }
          printf("\n");

          fclose(igb.fileptr);
          fclose(igb_out.fileptr);

          if(have_dynpts) {
            fclose(igb_pts_dynpt.fileptr);
          }
          break;
        }
        default: break;
      }
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      break;
    }

    case ITP_E2N:
    {
      std::cout << "Interpolating element data onto nodes .." << std::endl;
      gettimeofday(&t1, NULL);
      switch(data_idx) {
        case 0:
        case 2:
          read_vector_ascii(idat, opts.idat, true);

          if(data_idx == 0) {
            elemData_to_nodeData(omesh, idat, odat);
            write_vector_ascii(odat, opts.odat, 1);
          }
          else {
            array_to_points(idat, idat_vec);
            elemData_to_nodeData(omesh, idat_vec, odat_vec);
            points_to_array(odat_vec, odat);
            write_vector_ascii(odat, opts.odat, 3);
          }
          break;

        case 1:
        case 3:
        case 4:
        {
          // set new spatial size and write output igb header
          igb_out.v_x = omesh.n2e_cnt.size();
          write_igb_header(igb_out);

          std::vector<std::vector<mt_real> > rbuff;
          printf("Processing igb time-slices: %s to %s: \n", opts.idat.c_str(), opts.odat.c_str());

          for(int t=0; t<igb.v_t; t++) {
            printf("\rcurrent time-slice %d / %d .. ", t+1, int(igb.v_t));
            fflush(stdout);

            read_igb_block(rbuff, 1, igb);
            idat.assign(rbuff[0].begin(), rbuff[0].end());

            if(data_idx == 1) {
              elemData_to_nodeData(omesh, idat, odat);
            }
            else if(data_idx == 3) {
              array_to_points(idat, idat_vec);
              elemData_to_nodeData(omesh, idat_vec, odat_vec);
              points_to_array(odat_vec, odat);
            }
            else {
              array_to_tensors(idat, idat_ten);
              elemData_to_nodeData(omesh, idat_ten, odat_ten);
              tensors_to_array(odat_ten, odat);
            }

            rbuff[0].assign(odat.begin(), odat.end());
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
    case ITP_N2E:
    {
      std::cout << "Interpolating node data onto elements .." << std::endl;
      gettimeofday(&t1, NULL);
      switch(data_idx) {
        case 0:
        case 2:
          read_vector_ascii(idat, opts.idat, true);

          if(data_idx == 0) {
            nodeData_to_elemData(omesh, idat, odat);
            write_vector_ascii(odat, opts.odat, 1);
          }
          else {
            array_to_points(idat, idat_vec);
            nodeData_to_elemData(omesh, idat_vec, odat_vec);
            points_to_array(odat_vec, odat);
            write_vector_ascii(odat, opts.odat, 3);
          }
          break;

        case 1:
        case 3:
        case 4:
        {
          // set new spatial size and write output igb header
          igb_out.v_x = omesh.n2e_cnt.size();
          write_igb_header(igb_out);

          std::vector<std::vector<mt_real> > rbuff;
          printf("Processing igb time-slices: %s to %s: \n", opts.idat.c_str(), opts.odat.c_str());

          for(int t=0; t<igb.v_t; t++) {
            printf("\rcurrent time-slice %d / %d .. ", t+1, int(igb.v_t));
            fflush(stdout);

            read_igb_block(rbuff, 1, igb);
            idat.assign(rbuff[0].begin(), rbuff[0].end());

            if(data_idx == 1) {
              nodeData_to_elemData(omesh, idat, odat);
            }
            else if(data_idx == 3) {
              array_to_points(idat, idat_vec);
              nodeData_to_elemData(omesh, idat_vec, odat_vec);
              points_to_array(odat_vec, odat);
            }
            else {
              array_to_tensors(idat, idat_ten);
              nodeData_to_elemData(omesh, idat_ten, odat_ten);
              tensors_to_array(odat_ten, odat);
            }

            rbuff[0].assign(odat.begin(), odat.end());
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

    default: break;
  }


}

