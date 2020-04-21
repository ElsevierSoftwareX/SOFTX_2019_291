#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include <string>

#include "mt_modes_base.h"


struct fibmap_options{
  std::string msh_base;
  std::string la_surf;
  std::string ra_surf;
  std::string la_epi_fib;
  std::string la_endo_fib;
  std::string ra_endo_fib;
  std::string la_sol;
  std::string la_grad;
  std::string ra_grad;
  std::string numfib;
};

static const std::string la_surf_par     = "-la_surf=";
static const std::string ra_surf_par     = "-ra_surf=";
static const std::string la_epi_fib_par  = "-la_epi_fib=";
static const std::string la_endo_fib_par = "-la_endo_fib=";
static const std::string ra_endo_fib_par = "-ra_endo_fib=";
static const std::string la_sol_par      = "-la_sol=";
static const std::string la_grad_par     = "-la_grad=";
static const std::string ra_grad_par     = "-ra_grad=";
static const std::string numfib_par     = "-numfib=";


void print_fibmap_help()
{
  fprintf(stderr, "fibmap: Map fiber directions defined on the atrial endocard onto the whole atrial volume.");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to the basename of the mesh.\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to the LA endocard surface file.\n", la_surf_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to the RA endocard surface file.\n", ra_surf_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to the LA epicard fiber file.\n",  la_epi_fib_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to the LA endocard fiber file.\n", la_endo_fib_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to the RA endocard fiber file.\n", ra_endo_fib_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to the LA laplace solution data.\n", la_sol_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to the LA laplace solution gradient data.\n", la_grad_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to the RA laplace solution gradient data.\n", ra_grad_par.c_str());
  fprintf(stderr, "%s[1|2]\t (input) number of fiber directions in output. Must be either 1 or 2.\n", numfib_par.c_str());
  fprintf(stderr, "\n");
}

int fibmap_parse_options(int argc, char** argv, struct fibmap_options & opts)
{
  if(argc < 2) {
    print_fibmap_help();
    return 1;
  }

  // parse all parameters -----------------------------------------------------------------
  for(int i=1; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh_par, opts.msh_base);
    if(!match) match = parse_param(param, la_surf_par, opts.la_surf);
    if(!match) match = parse_param(param, ra_surf_par, opts.ra_surf);
    if(!match) match = parse_param(param, la_epi_fib_par,  opts.la_epi_fib);
    if(!match) match = parse_param(param, la_endo_fib_par, opts.la_endo_fib);
    if(!match) match = parse_param(param, ra_endo_fib_par, opts.ra_endo_fib);
    if(!match) match = parse_param(param, la_sol_par, opts.la_sol);
    if(!match) match = parse_param(param, la_grad_par, opts.la_grad);
    if(!match) match = parse_param(param, ra_grad_par, opts.ra_grad);
    if(!match) match = parse_param(param, numfib_par, opts.numfib);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 3;
    }
  }
  fixBasename(opts.msh_base);

  // check if all relevant parameters have been set ---------------------------------------------------
  bool msh_ok  = opts.msh_base.size() > 0,
       surf_ok = opts.la_surf.size() > 0 && opts.ra_surf.size() > 0,
       fib_ok  = opts.la_epi_fib.size() > 0 && opts.la_endo_fib.size() > 0 && opts.ra_endo_fib.size() > 0,
       sol_ok  = opts.la_sol.size() > 0,
       grad_ok = opts.la_grad.size() > 0 && opts.ra_grad.size() > 0,
       numfib_ok = opts.numfib.size() > 0 && ( opts.numfib.compare("1") == 0 || opts.numfib.compare("2") == 0);

  if( !(msh_ok && surf_ok && fib_ok && sol_ok && grad_ok && numfib_ok) )
  {
    std::cerr << "Error: Insufficient parameters provided." << std::endl;
    print_fibmap_help();
    return 2;
  }

  return 0;
}

void expand_active_nodes(const mt_meshdata & mesh,
                         const mt_vector<mt_real> & la_sol,
                         mt_vector<bool> & traversed,
                         mt_int nidx,
                         std::set<mt_int> & active_nodes)
{
  bool la_base = la_sol[nidx] <= 1.0;

  mt_int start = mesh.n2n_dsp[nidx], stop = start + mesh.n2n_cnt[nidx];
  for(mt_int i=start; i<stop; i++)
  {
    mt_int idx = mesh.n2n_con[i];
    if(!traversed[idx]) {
      bool both_on_la = la_sol[idx] <= 1.0 && la_base;
      bool both_on_ra = la_sol[idx] > 1.0 && (!la_base);

      if(both_on_la || both_on_ra)
        active_nodes.insert(idx);
    }
  }
  active_nodes.erase(nidx);
}

mt_int max_grad_neighbour(const mt_meshdata & mesh,
                     const mt_vector<bool> & traversed,
                     const mt_vector<bool> & onSubset,
                     const mt_point<mt_real> & grad,
                     mt_int nidx)
{
  mt_int max_idx = -1;
  mt_real max_sca = -1.0;
  mt_point<mt_real> rp(mesh.xyz.data() + nidx*3);

  mt_int start = mesh.n2n_dsp[nidx], stop = start + mesh.n2n_cnt[nidx];
  for(mt_int i=start; i<stop; i++)
  {
    mt_int idx = mesh.n2n_con[i];
    if(onSubset[idx] && traversed[idx]) {
      mt_point<mt_real> p = rp - mt_point<mt_real>(mesh.xyz.data() + idx*3);
      p.normalize();

      mt_real sca = p.scaProd(grad);

      if(max_sca < sca) {
        max_sca = sca;
        max_idx = idx;
      }
    }
  }
  return max_idx;
}

mt_int max_grad_neighbour(const mt_meshdata & mesh,
                     const mt_vector<bool> & traversed,
                     const mt_point<mt_real> & grad,
                     mt_int nidx)
{
  mt_int max_idx = -1;
  mt_real max_sca = -1.0;
  mt_point<mt_real> rp(mesh.xyz.data() + nidx*3);

  mt_int start = mesh.n2n_dsp[nidx], stop = start + mesh.n2n_cnt[nidx];
  for(mt_int i=start; i<stop; i++)
  {
    mt_int idx = mesh.n2n_con[i];
    if(traversed[idx]) {
      mt_point<mt_real> p = rp - mt_point<mt_real>(mesh.xyz.data() + idx*3);
      p.normalize();

      mt_real sca = p.scaProd(grad);

      if(max_sca < sca) {
        max_sca = sca;
        max_idx = idx;
      }
    }
  }
  return max_idx;
}


int main(int argc, char** argv)
{
  struct timeval t1, t2;
  struct fibmap_options opts;
  struct mt_meshdata mesh, la_surf, ra_surf;
  mt_vector<mt_real> la_epi_lon, la_sol, la_grad, ra_grad;

  int ret = fibmap_parse_options(argc, argv, opts);
  if (ret != 0) return 1;

  // read input data =====================================================

  gettimeofday(&t1, NULL);
  std::cout << "Reading mesh: " << opts.msh_base << std::endl;
  readElements_general(mesh, opts.msh_base);
  readPoints_general(mesh.xyz, opts.msh_base);
  compute_full_mesh_connectivity(mesh);

  std::cout << "\nReading LA endo surface and fibers: " << opts.la_surf  << ", "
            << opts.la_endo_fib << std::endl;
  readElements(la_surf, opts.la_surf);
  compute_full_mesh_connectivity(la_surf);
  readFibers(la_surf.lon, la_surf.e2n_cnt.size(), opts.la_endo_fib);
  std::cout << "Reading LA epi fibers: " << opts.la_epi_fib << std::endl;
  readFibers(la_epi_lon, la_surf.e2n_cnt.size(), opts.la_epi_fib);

  std::cout << "\nReading RA endo surface and fibers: " << opts.ra_surf << ", "
            << opts.ra_endo_fib << std::endl;
  readElements(ra_surf, opts.ra_surf);
  compute_full_mesh_connectivity(ra_surf);
  readFibers(ra_surf.lon, ra_surf.e2n_cnt.size(), opts.ra_endo_fib);

  std::cout << "\nReading LA solution data: " << opts.la_sol << std::endl;
  read_vector_ascii(la_sol, opts.la_sol);
  std::cout << "Reading LA gradient data: " << opts.la_grad << std::endl;
  readPoints(la_grad, opts.la_grad);
  std::cout << "Reading RA gradient data: " << opts.ra_grad << std::endl;
  readPoints(ra_grad, opts.ra_grad);

  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  gettimeofday(&t1, NULL);
  std::cout << "Convert element based fiber data to node based data .." << std::endl;

  // element based fiber data
  mt_vector<mt_point<mt_real> > la_endo_fib_ele(la_surf.e2n_cnt.size()),
                                la_epi_fib_ele (la_surf.e2n_cnt.size()),
                                ra_endo_fib_ele(ra_surf.e2n_cnt.size());

  array_to_points(la_surf.lon, la_endo_fib_ele);
  array_to_points(la_epi_lon,  la_epi_fib_ele);
  array_to_points(ra_surf.lon, ra_endo_fib_ele);

  // node based fiber data
  mt_vector<mt_point<mt_real> > la_endo_fib_nod(la_surf.n2n_cnt.size()),
                                la_epi_fib_nod (la_surf.n2n_cnt.size()),
                                ra_endo_fib_nod(ra_surf.n2n_cnt.size());

  la_surf.xyz.assign(mesh.xyz.size(), mesh.xyz.data(), false);
  ra_surf.xyz.assign(mesh.xyz.size(), mesh.xyz.data(), false);

  elemData_to_nodeData(la_surf, la_epi_fib_ele, la_epi_fib_nod);
  elemData_to_nodeData(la_surf, la_endo_fib_ele, la_endo_fib_nod);
  elemData_to_nodeData(ra_surf, ra_endo_fib_ele, ra_endo_fib_nod);

  la_surf.xyz.assign(0, NULL, false);
  ra_surf.xyz.assign(0, NULL, false);

  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  size_t nnodes = mesh.n2n_cnt.size();

  // prepare gradient ================================================================
  mt_vector<mt_point<mt_real> > grad(nnodes);
  mt_vector<bool> la_flag(nnodes, false), ra_flag(nnodes, false);
  {
    std::set<mt_int> la_nodes, ra_nodes;

    for(size_t i=0; i<nnodes; i++)
    {
      if(la_sol[i] <= 1.0) {
        grad[i].x = la_grad[i*3+0];
        grad[i].y = la_grad[i*3+1];
        grad[i].z = la_grad[i*3+2];
        la_flag[i] = true;
        la_nodes.insert(i);
      }
      else {
        grad[i].x = ra_grad[i*3+0];
        grad[i].y = ra_grad[i*3+1];
        grad[i].z = ra_grad[i*3+2];
        ra_flag[i] = true;
        ra_nodes.insert(i);
      }
    }

    mt_vector<mt_int> la_nod_vec, ra_nod_vec;
    la_nod_vec.assign(la_nodes.begin(), la_nodes.end());
    ra_nod_vec.assign(ra_nodes.begin(), ra_nodes.end());

    smooth_data(mesh, la_nod_vec, la_flag, 10, mt_real(0.2), grad, true);
    smooth_data(mesh, ra_nod_vec, ra_flag, 10, mt_real(0.2), grad, true);
    for(size_t i=0; i<nnodes; i++)
      grad[i].normalize();
    smooth_data(mesh, la_nod_vec, la_flag, 50, mt_real(0.2), grad, true);
    smooth_data(mesh, ra_nod_vec, ra_flag, 50, mt_real(0.2), grad, true);
  }

  // the fibers and sheets we compute
  mt_vector<mt_point<mt_real> > fib_ele(mesh.e2n_cnt.size()), she_ele(mesh.e2n_cnt.size());
  // we can already compute the sheets, since they correspond to the gradient
  nodeData_to_elemData(mesh, grad, she_ele);

  #if 0
  {
    mt_vector<float> fib_out(she_ele.size()*3);
    for(size_t i=0; i<she_ele.size(); i++)
    {
      fib_out[i*3+0] = she_ele[i].x;
      fib_out[i*3+1] = she_ele[i].y;
      fib_out[i*3+2] = she_ele[i].z;
    }
    std::string fibfile = opts.msh_base + "used_grad.vec";
    std::cout << "Writing vectors " << fibfile << std::endl;
    write_vector_ascii(fib_out, fibfile, 3);
    std::cout << "Done .." << std::endl;
  }
  #endif

  // start with mapping algorithm =====================================================
  gettimeofday(&t1, NULL);
  std::cout << "Computing fibers .." << std::endl;

  std::set<mt_int> active_list;
  mt_vector<mt_int> active_list_vec(nnodes);

  size_t numfin = 0;
  short  numfib = opts.numfib.compare("1") == 0 ? 1 : 2;

  mt_vector<bool> fin_map(nnodes, false);
  // fibers in nodal representation
  mt_vector<mt_point<mt_real> > fib_nod (nnodes, mt_point<mt_real>(0, 0, 0));
  // auxiliary fibers holding the epicardium fibers
  mt_vector<mt_point<mt_real> > auxfib_nod (nnodes, mt_point<mt_real>(0, 0, 0));

  // initialization: first set all surface nodes as finished
  for(size_t i=0; i<la_surf.n2n_cnt.size(); i++)
  {
    // test if node is represented in the surface
    if(la_surf.n2n_cnt[i]) {
      fib_nod[i]     = la_endo_fib_nod[i];
      auxfib_nod[i]  = la_epi_fib_nod[i];
      fin_map[i]  = true;
      numfin++;
    }
  }
  for(size_t i=0; i<ra_surf.n2n_cnt.size(); i++)
  {
    // test if node is represented in the surface
    if(ra_surf.n2n_cnt[i]) {
      fib_nod[i]  = ra_endo_fib_nod[i];
      fin_map[i]  = true;
      numfin++;
    }
  }
  // then add all connected nodes into the active list
  for(size_t i=0; i<nnodes; i++)
  {
    if(fin_map[i])
    {
      int start = mesh.n2n_dsp[i], stop = start + mesh.n2n_cnt[i];
      for(int j=start; j<stop; j++)
      {
        int idx = mesh.n2n_con[j];
        if(!fin_map[idx])
          active_list.insert(idx);
      }
    }
  }

  // write start data
  #if 0
  mt_vector<float> fib_out(nnodes*numfib*3);
  {
    for(size_t i=0; i<fib_nod.size(); i++)
    {
      fib_out[i*3+0] = fib_nod[i].x;
      fib_out[i*3+1] = fib_nod[i].y;
      fib_out[i*3+2] = fib_nod[i].z;
    }
    std::string fibfile = opts.msh_base + "fibersStart.vec";
    std::cout << "Writing fibers " << fibfile << std::endl;
    write_vector_ascii(fib_out, fibfile, 3);
    std::cout << "Done .." << std::endl;
  }
  #endif

  mt_vector<short> upd(nnodes, 0);
  const short upd_thr = 10;

  // computation loop
  while(active_list.size() > 0)
  {
    active_list_vec.assign(active_list.begin(), active_list.end());

    std::cout << "Active list size: " << active_list_vec.size() << std::endl;

    for(size_t i=0; i<active_list_vec.size(); i++)
    {
      mt_int nidx = active_list_vec[i];
      mt_point<mt_real> & cg = grad[nidx];
      mt_int vidx = -1;
      upd[nidx]++;

      if(la_sol[nidx] <= 1.0)
      {
        // we are in left atrium
        if(upd[nidx] < upd_thr)
          vidx = max_grad_neighbour(mesh, fin_map, la_flag, cg, nidx);
        else
          vidx = max_grad_neighbour(mesh, fin_map, cg, nidx);

        if(vidx > -1) {
          // always copy the epicard fibers
          auxfib_nod[nidx] = auxfib_nod[vidx];

          // if we are past the mid-section the fibers are taken from the epicard
          if(la_sol[nidx] > 0.5)
            fib_nod[nidx] = auxfib_nod[vidx];
          else
            fib_nod[nidx] = fib_nod[vidx];

          fin_map[nidx] = true;
        }
      }
      else
      {
        // we are in right atrium
        if(upd[nidx] < upd_thr)
          vidx = max_grad_neighbour(mesh, fin_map, ra_flag, cg, nidx);
        else
          vidx = max_grad_neighbour(mesh, fin_map, cg, nidx);

        if(vidx > -1)
        {
          fib_nod[nidx] = fib_nod[vidx];
          fin_map[nidx] = true;
        }
      }

      if(fin_map[nidx])
      {
        // adapt fiber directionn w.r.t. gradient
        mt_point<mt_real> n = fib_nod[nidx].crossProd(cg);
        fib_nod[nidx] = n.crossProd(cg) * mt_real(-1.0);
        fib_nod[nidx].normalize();

        // check if value is ok
        bool xok = ! ( std::isnan(fib_nod[nidx].x) || std::isinf(fib_nod[nidx].x) );
        bool yok = ! ( std::isnan(fib_nod[nidx].y) || std::isinf(fib_nod[nidx].y) );
        bool zok = ! ( std::isnan(fib_nod[nidx].z) || std::isinf(fib_nod[nidx].z) );
        if( ! (xok && yok && zok) )
          std::cerr << "Warning: Bad fiber at node " << nidx << std::endl;

        expand_active_nodes(mesh, la_sol, fin_map, nidx, active_list);
        numfin++;
      }
    }

    if(active_list.size() == 0 && numfin < nnodes)
    {
      std::cerr << "Warning: not all nodes reached! " << nnodes - numfin <<
                   " left unreached!" << std::endl;
      std::cerr << "Inserting all unreached nodes into active list." << std::endl;
      for(size_t i=0; i<nnodes; i++)
        if(fin_map[i] == false) active_list.insert(i);
    }
  }

  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;


  #if 0
  {
    for(size_t i=0; i<fib_nod.size(); i++)
    {
      fib_out[i*3+0] = fib_nod[i].x;
      fib_out[i*3+1] = fib_nod[i].y;
      fib_out[i*3+2] = fib_nod[i].z;
    }
    std::string fibfile = opts.msh_base + "fibersEnd.vec";
    std::cout << "Writing fibers " << fibfile << std::endl;
    write_vector_ascii(fib_out, fibfile, 3);
    std::cout << "Done .." << std::endl;
  }
  #endif

  nodeData_to_elemData(mesh, fib_nod, fib_ele);
  mesh.lon.resize(fib_ele.size()*numfib*3);
  for(size_t i=0; i<fib_ele.size(); i++)
  {
    mesh.lon[i*3*numfib+0] = fib_ele[i].x;
    mesh.lon[i*3*numfib+1] = fib_ele[i].y;
    mesh.lon[i*3*numfib+2] = fib_ele[i].z;
    if(numfib == 2) {
      mesh.lon[i*6+3] = she_ele[i].x;
      mesh.lon[i*6+4] = she_ele[i].y;
      mesh.lon[i*6+5] = she_ele[i].z;
    }
  }

  std::string lonfile = opts.msh_base + CARPBIN_LON_EXT;
  if(file_exists(lonfile))
  {
    std::cout << "Writing fibers " << lonfile << std::endl;
    writeFibersBinary(mesh.lon, mesh.e2n_cnt.size(), lonfile);
    std::cout << "Done .." << std::endl;
  }
  else {
    lonfile = opts.msh_base + CARPTXT_LON_EXT;
    std::cout << "Writing fibers " << lonfile << std::endl;
    writeFibers(mesh.lon, mesh.e2n_cnt.size(), lonfile);
    std::cout << "Done .." << std::endl;
  }

  return 0;
}
