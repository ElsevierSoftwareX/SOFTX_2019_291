/**
* @file mesh_smoothing.cpp
* @brief Mesh smoothing algorithms and classes.
* @author Aurel Neic
* @version
* @date 2017-08-16
*/

#include "mt_utils_base.h"
#include "mesh_smoothing.h"
#include "topology_utils.h"
#include "mesh_quality.h"
#include "mesh_utils.h"
#include "lookup_table.hpp"

bool centroid_smoothing_kernel(const mt_meshdata & m,
                               const mt_real * xyz,
                               const mt_real* vols,
                               const mt_int nidx,
                               const mt_real sca,
                               vec3r & np)
{
  // full mesh connectivity needs to be set up
  assert(m.e2n_dsp.size() > 0);
  assert(m.n2n_cnt.size() > 0 && m.n2e_cnt.size() > 0);
  mt_real totvol = 0.0;
  //get the point to smooth
  np.get(xyz + 3*nidx);
  vec3r avrg(0,0,0);

  const mt_int start = m.n2e_dsp[nidx], stop = start + m.n2e_cnt[nidx];
  for(mt_int i=start; i < stop; i++) {
    const mt_int eidx = m.n2e_con[i];
    const mt_int* con = m.e2n_con.data() + m.e2n_dsp[eidx];
    vec3r bary = barycenter(m.e2n_cnt[eidx], con, xyz);

    const mt_real vol = vols[eidx];
    avrg += (bary * vol);
    totvol += vol;
  }
  if(totvol > 0.0) {
    avrg /= totvol;
    np += (np - avrg) * sca;
  }
  return (totvol > 0.0);
}

bool circumsphere_smoothing_kernel(const mt_meshdata & m,
                                   const mt_real * xyz,
                                   const mt_real* vols,
                                   const mt_int nidx,
                                   const mt_real sca,
                                   vec3r & np)
{
  // full mesh connectivity needs to be set up
  assert(m.e2n_dsp.size() > 0);
  assert(m.n2n_cnt.size() > 0 && m.n2e_cnt.size() > 0);
  mt_real totvol = 0.0;
  //get the point to smooth
  np.get(xyz + 3*nidx);
  vec3r avrg(0,0,0);

  const mt_int start = m.n2e_dsp[nidx], stop = start + m.n2e_cnt[nidx];
  for(mt_int i=start; i < stop; i++) {
    const mt_int eidx = m.n2e_con[i];
    const mt_int* con = m.e2n_con.data() + m.e2n_dsp[eidx];
    vec3r midpnt;
    circumsphere(m.etype[eidx], con, xyz, midpnt);
    const mt_real vol = vols[eidx];
    avrg.x += vol * midpnt.x;
    avrg.y += vol * midpnt.y;
    avrg.z += vol * midpnt.z;
    totvol += vol;
  }
  if(totvol > 0.0) {
    avrg /= totvol;
    np += (np - avrg) * sca;
  }
  return (totvol > 0.0);
}

bool laplace_smoothing_kernel(const mt_meshdata & m, const mt_real * xyz, const mt_int nidx, const mt_real sca, vec3r & np)
{
  const mt_int start = m.n2n_dsp[nidx], stop = start + m.n2n_cnt[nidx];
  vec3r avrg(0,0,0);
  mt_int numadd = 0;

  //get the point to smooth
  np.get(xyz + 3*nidx);
  vec3r op = np;

  mt_real x,y,z;
  for(mt_int i=start; i<stop; i++) {
    const mt_int cidx = m.n2n_con[i];
    #ifdef OPENMP
    #pragma omp atomic read
    x = xyz[cidx*3 + 0];
    #pragma omp atomic read
    y = xyz[cidx*3 + 1];
    #pragma omp atomic read
    z = xyz[cidx*3 + 2];
    #else
    x = xyz[cidx*3 + 0], y = xyz[cidx*3 + 1], z = xyz[cidx*3 + 2];
    #endif
    avrg.x += x, avrg.y += y, avrg.z += z;
    numadd++;
  }
  // instead of testing cidx != nidx, we remove nidx after the loop
  numadd--;
  avrg -= op;

  if(numadd) {
    avrg /= static_cast<mt_real>(numadd);
    np += (avrg - np) * sca;
  }
  return (numadd > 0);
}


void smooth_iter(mt_meshdata & mesh,
                 const mt_vector<mt_real> & old_xyz,
                 const mt_vector<mt_real> & vols,
                 const mt_vector<mt_real> & vols_mnfld,
                 const mt_meshdata & manifold,
                 const mt_mask & mnfld_mask,
                 const mt_vector<mt_int> & nodes,
                 mt_real sca,
                 smoothing_t type)
{
  size_t nnodes = nodes.size();
  const mt_real *__restrict xyz = old_xyz.data();

  #ifdef OPENMP
  #pragma omp parallel for schedule(dynamic, 100)
  #endif
  for(size_t k=0; k < nnodes; k++)
  {
    const mt_int nidx = nodes[k];
    mt_point<mt_real> np;

    const mt_meshdata &        msh = mnfld_mask.count(nidx) ? manifold : mesh;
    const mt_vector<mt_real> & vol = mnfld_mask.count(nidx) ? vols_mnfld : vols;

    bool did_smth = false;
    if(type == LAPLACE_SMOOTHING)
      did_smth = laplace_smoothing_kernel(msh, xyz, nidx, sca, np);
    else if(type == CIRCUMSPHERE_SMOOTHING)
      did_smth = circumsphere_smoothing_kernel(msh, xyz, vol.data(), nidx, sca, np);
    else if(type == CENTROID_SMOOTHING)
      did_smth = centroid_smoothing_kernel(msh, xyz, vol.data(), nidx, sca, np);

    if(did_smth) {
      mesh.xyz[nidx*3+0] = np.x;
      mesh.xyz[nidx*3+1] = np.y;
      mesh.xyz[nidx*3+2] = np.z;
    }
  }
}

void quality_aware_smoother::init_qual_thr(const mt_meshdata & checkmesh,
                   const mt_vector<mt_int> & nodes,
                   const mt_real inp_thr)
{
  const size_t nnodes = nodes.size();
  _qual_thr.resize(nnodes);

  #ifdef OPENMP
  #pragma omp parallel for schedule(dynamic, 100)
  #endif
  for(size_t k=0; k < nnodes; k++)
  {
    const mt_int nidx = nodes[k];
    mt_real qual; mt_int qual_eidx;

    nbhd_quality_max(checkmesh, nidx, qual, qual_eidx);
    mt_meshdata manifold;  // empty manifold mesh to satisfy API
    mt_mask isSurf;        // an empty mask to satisfy API

    if(qual != qual) {
      _qual_thr[k] = inp_thr;
    }
    else if(qual < QUAL_THR_CHK_SCA * inp_thr) {
      _qual_thr[k] = 0.0;
    }
    else {
      _qual_thr[k] = qual < inp_thr ? inp_thr : qual;
    }
  }
}

void quality_aware_smoother::qual_smooth_iter_fwd(mt_meshdata & checkmesh,
                          mt_meshdata & smoothmesh,
                          mt_meshdata & manifold,
                          const mt_mask & mnfld_mask,
                          const mt_vector<mt_int> & nodes,
                          const mt_real sca,
                          const mt_real inp_thr)
{
  const size_t nnodes = nodes.size();
  mt_real *__restrict xyz = checkmesh.xyz.data();

  #ifdef OPENMP
  #pragma omp parallel for schedule(guided)
  #endif
  for(size_t k=0; k < nnodes; k++)
  {
    const mt_int nidx = nodes[k];
    mt_point<mt_real> np, op(xyz + nidx*3);

    const mt_meshdata &        msh = mnfld_mask.count(nidx) ? manifold : smoothmesh;
    const mt_vector<mt_real> & vol = mnfld_mask.count(nidx) ? _vols_mnfld : _vols_msh;

    bool did_smth = false;
    if(smoothing_type == LAPLACE_SMOOTHING)
      did_smth = laplace_smoothing_kernel(msh, xyz, nidx, sca, np);
    else if(smoothing_type == CIRCUMSPHERE_SMOOTHING)
      did_smth = circumsphere_smoothing_kernel(msh, xyz, vol.data(), nidx, sca, np);
    else if(smoothing_type == CENTROID_SMOOTHING)
      did_smth = centroid_smoothing_kernel(msh, xyz, vol.data(), nidx, sca, np);

    if(did_smth) {
      #ifdef OPENMP
      # pragma omp atomic write
      xyz[nidx*3+0] = np.x;
      # pragma omp atomic write
      xyz[nidx*3+1] = np.y;
      # pragma omp atomic write
      xyz[nidx*3+2] = np.z;
      #else
      xyz[nidx*3+0] = np.x, xyz[nidx*3+1] = np.y, xyz[nidx*3+2] = np.z;
      #endif
    }
    _applied[k] = true;
    const double qual_thr = _qual_thr[k];

    if(qual_thr > 0.0) {
      const double qual = nbhd_quality_max(checkmesh, nidx);

      if(qual < qual_thr) {
        if(qual_thr > inp_thr)
          _qual_thr[k] = qual;
      }
      else {
        #ifdef OPENMP
        # pragma omp atomic write
        xyz[nidx*3+0] = op.x;
        # pragma omp atomic write
        xyz[nidx*3+1] = op.y;
        # pragma omp atomic write
        xyz[nidx*3+2] = op.z;
        #else
        xyz[nidx*3+0] = op.x, xyz[nidx*3+1] = op.y, xyz[nidx*3+2] = op.z;
        #endif
        _applied[k] = false;
      }
    }
  }
}

void quality_aware_smoother::qual_smooth_iter_fwd(const mt_vector<mt_meshdata*> & manifolds,
                                                  const mt_vector<int> & mnfld_idx,
                                                  mt_vector<mt_real> & inp_xyz,
                                                  const mt_vector<mt_int> & nodes,
                                                  const mt_real sca,
                                                  const mt_real inp_thr)
{
  const size_t nnodes = nodes.size();
  mt_real *__restrict xyz = inp_xyz.data();

  // we check the mesh quality w.r.t. the first mesh in manifolds.
  const mt_meshdata & checkmesh = *manifolds[0];

  #ifdef OPENMP
  #pragma omp parallel for schedule(guided)
  #endif
  for(size_t k=0; k < nnodes; k++)
  {
    const mt_int nidx = nodes[k];
    mt_point<mt_real> avrg(0, 0, 0), np(xyz + nidx*3);
    mt_point<mt_real> op = np;

    // each node is associated to a manifold graph
    const mt_meshdata & msh = *manifolds[mnfld_idx[nidx]];
    bool did_smth = laplace_smoothing_kernel(msh, xyz, nidx, sca, np);

    if(did_smth) {
      #ifdef OPENMP
      # pragma omp atomic write
      xyz[nidx*3+0] = np.x;
      # pragma omp atomic write
      xyz[nidx*3+1] = np.y;
      # pragma omp atomic write
      xyz[nidx*3+2] = np.z;
      #else
      xyz[nidx*3+0] = np.x, xyz[nidx*3+1] = np.y, xyz[nidx*3+2] = np.z;
      #endif
    }
    _applied[k] = true;
    const double qual_thr = _qual_thr[k];

    if(qual_thr > 0.0) {
      const double qual = nbhd_quality_max(checkmesh, nidx);

      if(qual < qual_thr) {
        if(qual_thr > inp_thr)
          _qual_thr[k] = qual;
      }
      else {
        #ifdef OPENMP
        # pragma omp atomic write
        xyz[nidx*3+0] = op.x;
        # pragma omp atomic write
        xyz[nidx*3+1] = op.y;
        # pragma omp atomic write
        xyz[nidx*3+2] = op.z;
        #else
        xyz[nidx*3+0] = op.x, xyz[nidx*3+1] = op.y, xyz[nidx*3+2] = op.z;
        #endif
        _applied[k] = false;
      }
    }
  }
}




void quality_aware_smoother::qual_smooth_iter_bwd(mt_meshdata & checkmesh,
                          mt_meshdata & smoothmesh,
                          mt_meshdata & manifold,
                          const mt_mask & mnfld_mask,
                          const mt_vector<mt_int> & nodes,
                          mt_real sca,
                          mt_real inp_thr)
{
  size_t nnodes = nodes.size();
  mt_real *__restrict xyz = checkmesh.xyz.data();

  #ifdef OPENMP
  #pragma omp parallel for schedule(guided)
  #endif
  for(size_t k=0; k < nnodes; k++)
  {
    const mt_int nidx = nodes[k];

    if(_applied[k]) {
      mt_point<mt_real> avrg(0, 0, 0), np(xyz + nidx*3);
      mt_point<mt_real> op = np;

      const mt_meshdata &        msh = mnfld_mask.count(nidx) ? manifold : smoothmesh;
      const mt_vector<mt_real> & vol = mnfld_mask.count(nidx) ? _vols_mnfld : _vols_msh;

      bool did_smth = false;
      if(smoothing_type == LAPLACE_SMOOTHING)
        did_smth = laplace_smoothing_kernel(msh, xyz, nidx, sca, np);
      else if(smoothing_type == CIRCUMSPHERE_SMOOTHING)
        did_smth = circumsphere_smoothing_kernel(msh, xyz, vol.data(), nidx, sca, np);
      else if(smoothing_type == CENTROID_SMOOTHING)
        did_smth = centroid_smoothing_kernel(msh, xyz, vol.data(), nidx, sca, np);

      if(did_smth) {
        #ifdef OPENMP
        # pragma omp atomic write
        xyz[nidx*3+0] = np.x;
        # pragma omp atomic write
        xyz[nidx*3+1] = np.y;
        # pragma omp atomic write
        xyz[nidx*3+2] = np.z;
        #else
        xyz[nidx*3+0] = np.x, xyz[nidx*3+1] = np.y, xyz[nidx*3+2] = np.z;
        #endif
      }
      _applied[k] = false;
      const double qual_thr = _qual_thr[k];

      if(qual_thr > 0.0) {
        const double qual = nbhd_quality_max(checkmesh, nidx);

        if(qual < qual_thr) {
          if(qual_thr > inp_thr)
            _qual_thr[k] = qual;
        }
        else {
          #ifdef OPENMP
          # pragma omp atomic write
          xyz[nidx*3+0] = op.x;
          # pragma omp atomic write
          xyz[nidx*3+1] = op.y;
          # pragma omp atomic write
          xyz[nidx*3+2] = op.z;
          #else
          xyz[nidx*3+0] = op.x, xyz[nidx*3+1] = op.y, xyz[nidx*3+2] = op.z;
          #endif
        }
      }
    }
  }
}

void quality_aware_smoother::qual_smooth_iter_bwd(const mt_vector<mt_meshdata*> & manifolds,
                                                  const mt_vector<int> & mnfld_idx,
                                                  mt_vector<mt_real> & inp_xyz,
                                                  const mt_vector<mt_int> & nodes,
                                                  const mt_real sca,
                                                  const mt_real inp_thr)
{
  size_t nnodes = nodes.size();
  mt_real *__restrict xyz = inp_xyz.data();

  // we check the mesh quality w.r.t. the first mesh in manifolds.
  const mt_meshdata & checkmesh = *manifolds[0];

  #ifdef OPENMP
  #pragma omp parallel for schedule(guided)
  #endif
  for(size_t k=0; k < nnodes; k++)
  {
    const mt_int nidx = nodes[k];

    if(_applied[k]) {
      mt_point<mt_real> avrg(0, 0, 0), np(xyz + nidx*3);
      mt_point<mt_real> op = np;

      // each node is associated to a manifold graph
      const mt_meshdata & msh = *manifolds[mnfld_idx[nidx]];
      bool did_smth = laplace_smoothing_kernel(msh, xyz, nidx, sca, np);

      if(did_smth) {
        #ifdef OPENMP
        # pragma omp atomic write
        xyz[nidx*3+0] = np.x;
        # pragma omp atomic write
        xyz[nidx*3+1] = np.y;
        # pragma omp atomic write
        xyz[nidx*3+2] = np.z;
        #else
        xyz[nidx*3+0] = np.x, xyz[nidx*3+1] = np.y, xyz[nidx*3+2] = np.z;
        #endif
      }
      _applied[k] = false;
      const double qual_thr = _qual_thr[k];

      if(qual_thr > 0.0) {
        const double qual = nbhd_quality_max(checkmesh, nidx);

        if(qual < qual_thr) {
          if(qual_thr > inp_thr)
            _qual_thr[k] = qual;
        }
        else {
          #ifdef OPENMP
          # pragma omp atomic write
          xyz[nidx*3+0] = op.x;
          # pragma omp atomic write
          xyz[nidx*3+1] = op.y;
          # pragma omp atomic write
          xyz[nidx*3+2] = op.z;
          #else
          xyz[nidx*3+0] = op.x, xyz[nidx*3+1] = op.y, xyz[nidx*3+2] = op.z;
          #endif
        }
      }
    }
  }
}


void quality_aware_smoother::operator()(mt_meshdata & checkmesh,
                mt_meshdata & smoothmesh,
                mt_meshdata & manifold,
                const mt_mask & mnfld_mask,
                const mt_vector<mt_int> & nodes,
                const size_t nsmooth,
                const mt_real sca,
                const mt_real inp_thr)
{
  _applied.assign(nodes.size(), false);
  init_qual_thr(checkmesh, nodes, inp_thr);

  bool need_vols = false;
  if(smoothing_type != LAPLACE_SMOOTHING)
    need_vols = true;

  if(need_vols) {
    mesh_volumes(smoothmesh, smoothmesh.xyz, _vols_msh);
    mesh_volumes(manifold,   smoothmesh.xyz, _vols_mnfld);
  }

  std::string msg = get_smoothing_type_string(smoothing_type);
  msg += " smoothing progress: ";

  PROGRESS<size_t> progress(nsmooth, msg.c_str());

  for(size_t sidx = 0; sidx < nsmooth; sidx++) {
    progress.next();

    qual_smooth_iter_fwd(checkmesh, smoothmesh, manifold, mnfld_mask, nodes, sca, inp_thr);
    qual_smooth_iter_bwd(checkmesh, smoothmesh, manifold, mnfld_mask, nodes, -(sca*SMOOTH_FREQ_SCA), inp_thr);

    if((sidx % QUAL_THR_UPD) == 0)
      init_qual_thr(checkmesh, nodes, inp_thr);

    if(need_vols && (sidx % VOL_UPD) == 0) {
      mesh_volumes(smoothmesh, smoothmesh.xyz, _vols_msh);
      mesh_volumes(manifold,   smoothmesh.xyz, _vols_mnfld);
    }
  }
  progress.finish();
}

void quality_aware_smoother::operator()(const mt_vector<mt_meshdata*> & manifolds,
                                                  const mt_vector<int> & mnfld_idx,
                                                  mt_vector<mt_real> & inp_xyz,
                                                  const mt_vector<mt_int> & nodes,
                                                  const size_t  nsmooth,
                                                  const mt_real sca,
                                                  const mt_real inp_thr)
{
  // we check the mesh quality w.r.t. the first mesh in manifolds.
  mt_meshdata & checkmesh    = *manifolds[0];
  bool  checkmesh_had_no_xyz = false;

  if(checkmesh.xyz.size() != inp_xyz.size()) {
    // if the checkmesh has no coords, we do a shallow copy of the input xyz
    if(checkmesh.xyz.size()) {
      fprintf(stderr, "%s error: Coords of base manifold do not match input coords! Aborting!\n",
              __func__);
      exit(EXIT_FAILURE);
    }

    checkmesh.xyz.assign(inp_xyz.size(), inp_xyz.data(), false);
    checkmesh_had_no_xyz = true;
  }

  _applied.assign(nodes.size(), false);
  init_qual_thr(checkmesh, nodes, inp_thr);

  PROGRESS<size_t> progress(nsmooth, "Smoothing progress: ");

  for(size_t sidx = 0; sidx < nsmooth; sidx++) {
    progress.next();

    qual_smooth_iter_fwd(manifolds, mnfld_idx, inp_xyz, nodes, sca, inp_thr);
    qual_smooth_iter_bwd(manifolds, mnfld_idx, inp_xyz, nodes, -(sca*SMOOTH_FREQ_SCA), inp_thr);

    if((sidx % QUAL_THR_UPD) == 0)
      init_qual_thr(checkmesh, nodes, inp_thr);
  }
  progress.finish();

  if(checkmesh_had_no_xyz)
    checkmesh.xyz.assign(size_t(0), NULL, false);
}



void smooth_nodes(mt_meshdata & mesh,
                  mt_meshdata & manifold,
                  const mt_mask & mnfld_mask,
                  const mt_vector<mt_int> & nodes,
                  size_t nsmooth,
                  mt_real sca,
                  smoothing_t type)
{
  mt_vector<mt_real> old_xyz, vols, vols_mnfld;

  std::string msg = get_smoothing_type_string(type);
  msg += " smoothing progress: ";

  PROGRESS<size_t> progress(nsmooth, msg.c_str());

  bool need_vols = false;
  if(type != LAPLACE_SMOOTHING)
    need_vols = true;

  if(need_vols) {
    mesh_volumes(mesh,     mesh.xyz, vols);
    mesh_volumes(manifold, mesh.xyz, vols_mnfld);
  }

  for(size_t sidx = 0; sidx < nsmooth; sidx++)
  {
    progress.next();

    old_xyz.assign(mesh.xyz.begin(), mesh.xyz.end());
    smooth_iter(mesh, old_xyz, vols, vols_mnfld, manifold, mnfld_mask, nodes, sca, type);
    old_xyz.assign(mesh.xyz.begin(), mesh.xyz.end());
    smooth_iter(mesh, old_xyz, vols, vols_mnfld, manifold, mnfld_mask, nodes, -(sca*SMOOTH_FREQ_SCA), type);

    if(need_vols && (sidx % VOL_UPD) == 0) {
      mesh_volumes(mesh,     mesh.xyz, vols);
      mesh_volumes(manifold, mesh.xyz, vols_mnfld);
    }
  }
  progress.finish();
}

void volumetric_smooth_from_tags(mt_meshdata & mesh,
                                 const mt_vector<MT_USET<mt_int> > & stags,
                                 const int iter,
                                 const mt_real smth,
                                 const mt_real edge_ang,
                                 const mt_real max_qual,
                                 const bool skip_lines,
                                 const bool verbose)
{
  mt_meshdata surfmesh, linemesh;
  MT_USET<mt_int> sm_vtx, ln_vtx, bad_surf_vtx, bad_ln_vtx;
  int viter = iter*0.5, siter = iter*0.5;

  // extract surface and line manifolds
  if(verbose)
    std::cout << "Processing tag surfaces into unified surface .." << std::endl;

  unified_surface_from_tags(mesh, stags, surfmesh, &sm_vtx);
  compute_full_mesh_connectivity(surfmesh);

  // check if we are working on an open surface
  identify_surface_border_nodes(surfmesh, bad_surf_vtx);
  if(bad_surf_vtx.size()) {
    printf("%s warning: The derived surface to smooth is open! "
           "Skipping surface border smoothing.\n", __func__);
  }

  if(verbose)
    std::cout << "Computing line interfaces between surfaces .." << std::endl;

  compute_line_interfaces(mesh, surfmesh, edge_ang, false, linemesh);
  compute_full_mesh_connectivity(linemesh);

  // line vertices connected to more than two neighbours are corners and not proper
  // line manifold vertices -> smoothing does not work for them. So we remove from the
  // smoothing vertex set later on. Note that we check n2n_cnt > 3 since one of the
  // edges in the n2n graph is the node itself.
  for(size_t i=0; i<linemesh.n2n_cnt.size(); i++)
    if(linemesh.n2n_cnt[i] > 3) bad_ln_vtx.insert(i);

  ln_vtx.insert(linemesh.e2n_con.begin(), linemesh.e2n_con.end());

  // remove lines from smoothing
  sm_vtx.sort(); ln_vtx.sort();
  for(auto n : ln_vtx) sm_vtx.erase(n);

  // remove open surface border from smoothing
  for(auto n : bad_surf_vtx) sm_vtx.erase(n);

  mt_vector<mt_int> sm_nod;
  mt_mask isMnfld(mesh.xyz.size() / 3);
  sm_nod.assign(sm_vtx.begin(), sm_vtx.end());
  isMnfld.insert(surfmesh.e2n_con.begin(), surfmesh.e2n_con.end());

  if(verbose)
    std::cout << "Smoothing .." << std::endl;

  // smooth
  quality_aware_smoother smoother;

  if(max_qual) {
    smoother.smoothing_type = LAPLACE_SMOOTHING;
    smoother(mesh, mesh, surfmesh, isMnfld, sm_nod, viter, smth, max_qual);
  }
  else {
    smooth_nodes(mesh, surfmesh, isMnfld, sm_nod, viter, smth, LAPLACE_SMOOTHING);
  }

  // prepare data-structs for line smoothing
  sm_vtx.clear(); isMnfld.clear();
  sm_vtx.insert(surfmesh.e2n_con.begin(), surfmesh.e2n_con.end());
  for(auto n : bad_ln_vtx)   sm_vtx.erase(n);
  for(auto n : bad_surf_vtx) sm_vtx.erase(n);

  if(skip_lines) {
    for(const mt_int & n : ln_vtx) sm_vtx.erase(n);
  }

  sm_nod.assign(sm_vtx.begin(), sm_vtx.end());
  isMnfld.insert(ln_vtx.begin(), ln_vtx.end());

  // a shallow copy of coords into surfmesh
  surfmesh.xyz.assign(mesh.xyz.size(), mesh.xyz.data(), false);

  if(max_qual) {
    smoother.smoothing_type = LAPLACE_SMOOTHING;
    smoother(mesh, surfmesh, linemesh, isMnfld, sm_nod, siter, smth, max_qual);
  }
  else {
    smooth_nodes(surfmesh, linemesh, isMnfld, sm_nod, siter, smth, LAPLACE_SMOOTHING);
  }

  surfmesh.xyz.assign(0, NULL, false);
}


void surface_smooth(mt_meshdata & surfmesh,
                    const int iter,
                    const mt_real smth,
                    const mt_real edge_ang,
                    const mt_real max_qual,
                    const bool skip_lines,
                    const bool verbose)
{
  mt_meshdata linemesh;
  MT_USET<mt_int> sm_vtx, ln_vtx, bad_ln_vtx;

  if(verbose)
    std::cout << "Computing line interfaces between surfaces .." << std::endl;

  compute_line_interfaces(surfmesh, surfmesh, edge_ang, false, linemesh);
  compute_full_mesh_connectivity(linemesh);

  // line vertices connected to more than two neighbours are corners and not proper
  // line manifold vertices -> smoothing does not work for them. So we remove from the
  // smoothing vertex set later on. Note that we check n2n_cnt > 3 since one of the
  // edges in the n2n graph is the node itself.
  for(size_t i=0; i<linemesh.n2n_cnt.size(); i++)
    if(linemesh.n2n_cnt[i] > 3) bad_ln_vtx.insert(i);

  ln_vtx.insert(linemesh.e2n_con.begin(), linemesh.e2n_con.end());
  ln_vtx.sort();

  mt_vector<mt_int> sm_nod;
  mt_mask isMnfld(surfmesh.xyz.size() / 3);

  if(verbose) std::cout << "Smoothing .." << std::endl;
  quality_aware_smoother smoother;

  // prepare data-structs for line smoothing
  sm_vtx.insert(surfmesh.e2n_con.begin(), surfmesh.e2n_con.end());
  for(auto n : bad_ln_vtx) sm_vtx.erase(n);

  if(skip_lines) {
    for(const mt_int & n : ln_vtx) sm_vtx.erase(n);
  }

  sm_nod.assign(sm_vtx.begin(), sm_vtx.end());
  isMnfld.insert(ln_vtx.begin(), ln_vtx.end());

  if(max_qual) smoother(surfmesh, surfmesh, linemesh, isMnfld, sm_nod, iter, smth, max_qual);
  else         smooth_nodes(surfmesh, linemesh, isMnfld, sm_nod, iter, smth);
}


void directional_smoothing(mt_meshdata & mesh,
                           const mt_vector<mt_int> & nod,
                           const mt_vector<mt_real> & scale_dir,
                           mt_real scale,
                           mt_real ortho_scale,
                           const bool only_positive,
                           const mt_real smth,
                           const int iter)
{
  size_t nnodes = nod.size();

  for(int sidx = 0; sidx < iter; sidx++) {
    #ifdef OPENMP
    #pragma omp parallel for schedule(dynamic, 100)
    #endif
    for(size_t k=0; k < nnodes; k++)
    {
      const mt_int nidx = nod[k];
      vec3r avrg(mesh.xyz.data() + nidx*3);
      vec3r np(mesh.xyz.data() + nidx*3);
      //This way np contains np = np + 1.0 * (avrg - np) = avrg
      bool did_smth = laplace_smoothing_kernel(mesh, mesh.xyz.data(), nidx, 1.0, avrg);
      if(did_smth) {
        vec3r dir = avrg - np;       // the direction towards the average
        mt_real len = dir.length();  // the step length towards the average
        dir /= len;

        // we (at least partially) orthogonalize w.r.t. to scale_dir
        vec3r sdir(scale_dir.data() + nidx*3);

        // we compute the vectors spanning the orthogonal plane
        vec3r odir1 = sdir.x < 0.95 ? vec3r(1,0,0) : vec3r(0,1,0);
        odir1 -= sdir * (sdir.scaProd(odir1)); odir1.normalize();
        vec3r odir2 = sdir.crossProd(odir1);

        mt_real sdir_proj = sdir.scaProd(dir);

        // if only_positive == true, we allow only positive movement, thus
        // for sdir_proj < 0, we set the scale to 1.0 regardless of user
        // setting, in order to remove that component from dir
        if(only_positive && sdir_proj < 0.0)
          dir -= sdir * sdir_proj;
        else
          dir -= sdir * (sdir_proj * scale);

        // now we (partially) remove the orthogonal directions
        dir -= odir1 * (odir1.scaProd(dir) * ortho_scale);
        dir -= odir2 * (odir2.scaProd(dir) * ortho_scale);

        // we reapply the original step length
        mt_real newlen = dir.length();
        if(newlen > 1e-10) {
          dir = dir * len / newlen;
          np += dir * smth;

          mesh.xyz[nidx*3+0] = np.x;
          mesh.xyz[nidx*3+1] = np.y;
          mesh.xyz[nidx*3+2] = np.z;
        }
      }
    }
  }
}


void directional_smoothing_elem(mt_meshdata & mesh,
                           const mt_vector<mt_int> & nod,
                           const mt_vector<mt_real> & ortho_dir,
                           mt_real ortho_scale,
                           const bool only_positive,
                           const mt_real smth,
                           const int iter)
{
  size_t nnodes = nod.size();

  for(int sidx = 0; sidx < iter; sidx++) {
    #ifdef OPENMP
    #pragma omp parallel for schedule(dynamic, 100)
    #endif
    for(size_t k=0; k < nnodes; k++)
    {
      const mt_int nidx = nod[k];
      vec3r avrg(mesh.xyz.data() + nidx*3);
      vec3r np(mesh.xyz.data() + nidx*3);
      bool did_smth = laplace_smoothing_kernel(mesh, mesh.xyz.data(), nidx, 1.0, avrg);
      if(did_smth) {
        vec3r dir = avrg - np;       // the direction towards the average
        mt_real len = dir.length();  // the step length towards the average
        dir /= len;

        for(mt_int i=0; i<mesh.n2e_cnt[nidx]; i++) {
          mt_int eidx = mesh.n2e_con[mesh.n2e_dsp[nidx]+i];

          // we (at least partially) orthogonalize w.r.t. to ortho_dir
          vec3r ortho(ortho_dir.data() + eidx*3);
          mt_real ortho_proj = ortho.scaProd(dir);

          // if only_positive == true, we allow only positive movement, thus
          // for ortho_proj < 0, we set the ortho_scale to 1.0 regardless of user
          // setting, in order to remove that component from dir
          if(only_positive && ortho_proj < 0.0)
            dir -= ortho * ortho_proj;
          else
            dir -= ortho * (ortho_proj * ortho_scale);

        }

        // we reapply the original step length
        mt_real newlen = dir.length();

        if(newlen > 1e-10) {
          dir = dir * len / newlen;
          np += dir * smth;

          mesh.xyz[nidx*3+0] = np.x;
          mesh.xyz[nidx*3+1] = np.y;
          mesh.xyz[nidx*3+2] = np.z;
        }
      }
    }
  }
}

