/**
* @file mesh_refining.cpp
* @brief Mesh refining classes.
* @author Aurel Neic
* @version
* @date 2017-08-07
*/


#include "mt_utils_base.h"
#include "topology_utils.h"
#include "mesh_utils.h"
#include "mesh_refining.h"
#include "mesh_quality.h"

void mt_edge_splitter::split_line(const mt_int eidx,
    const mt_int ele_start,
    const mt_int con_start,
    const tuple<mt_int> & edge,
    const mt_int nv)
{
  mt_meshdata & mesh = _edge_manager._mesh;
  short numfib = mesh.lon.size() == mesh.e2n_cnt.size()*6 ? 6 : 3;

  mt_int v1 = mesh.e2n_con[mesh.e2n_dsp[eidx] + 0],
  v2 = mesh.e2n_con[mesh.e2n_dsp[eidx] + 1];

  // generate connectivities
  mesh.e2n_con[con_start + 0] = v1;
  mesh.e2n_con[con_start + 1] = nv;
  mesh.e2n_con[con_start + 2] = nv;
  mesh.e2n_con[con_start + 3] = v2;

  // copy element data
  mesh.e2n_cnt[ele_start + 0] = 2;
  mesh.e2n_cnt[ele_start + 1] = 2;
  mesh.etype[ele_start + 0] = Line;
  mesh.etype[ele_start + 1] = Line;
  mesh.etags[ele_start + 0] = mesh.etags[eidx];
  mesh.etags[ele_start + 1] = mesh.etags[eidx];

  memcpy(mesh.lon.data() + ele_start*numfib,     mesh.lon.data() + eidx*numfib, numfib*sizeof(mt_real));
  memcpy(mesh.lon.data() + (ele_start+1)*numfib, mesh.lon.data() + eidx*numfib, numfib*sizeof(mt_real));
}

// Split a triangle into 2 triangles along a given edge.
void mt_edge_splitter::split_tri(const mt_int eidx,
    const mt_int ele_start,
    const mt_int con_start,
    const tuple<mt_int> & edge,
    const mt_int nv)
{
  mt_meshdata & mesh = _edge_manager._mesh;
  short numfib = mesh.lon.size() == mesh.e2n_cnt.size()*6 ? 6 : 3;

  mt_int b1 = edge.v1, b2 = edge.v2, v;
  mt_int v1 = mesh.e2n_con[mesh.e2n_dsp[eidx] + 0],
  v2 = mesh.e2n_con[mesh.e2n_dsp[eidx] + 1],
  v3 = mesh.e2n_con[mesh.e2n_dsp[eidx] + 2];

  if(v1 != b1 && v1 != b2)      v = v1;
  else if(v2 != b1 && v2 != b2) v = v2;
  else                          v = v3;

  /*
   * the new triangles are (b1, nv, v) and (b2, v, nv)
   *     v
   *    .    .
   *   .       .
   *  .           .
   * b1 . .  *nv . .b2
   */
  // generate connectivities
  mesh.e2n_con[con_start + 0] = b1;
  mesh.e2n_con[con_start + 1] = nv;
  mesh.e2n_con[con_start + 2] = v;
  mesh.e2n_con[con_start + 3] = b2;
  mesh.e2n_con[con_start + 4] = v;
  mesh.e2n_con[con_start + 5] = nv;

  // copy element data
  mesh.e2n_cnt[ele_start + 0] = 3;
  mesh.e2n_cnt[ele_start + 1] = 3;
  mesh.etype[ele_start + 0] = Tri;
  mesh.etype[ele_start + 1] = Tri;
  mesh.etags[ele_start + 0] = mesh.etags[eidx];
  mesh.etags[ele_start + 1] = mesh.etags[eidx];

  memcpy(mesh.lon.data() + ele_start*numfib,     mesh.lon.data() + eidx*numfib, numfib*sizeof(mt_real));
  memcpy(mesh.lon.data() + (ele_start+1)*numfib, mesh.lon.data() + eidx*numfib, numfib*sizeof(mt_real));
}


void mt_edge_splitter::split_tet(const mt_int eidx,
    const mt_int ele_start,
    const mt_int con_start,
    const tuple<mt_int> & edge,
    const mt_int nv)
{
  mt_meshdata & mesh = _edge_manager._mesh;
  short numfib = mesh.lon.size() == mesh.e2n_cnt.size()*6 ? 6 : 3;

  /*
   * the new tetrahedra are (va, b1, vb, nv) and (va, nv, vb, b2)
   *     b2
   *    .  .  .
   *   .    .    .
   *  .     *nv    .
   * va . -  . -  - .vb
   *      .   .   .
   *         . b1
   */

  mt_int b1 = edge.v1, b2 = edge.v2, va, vb;
  mt_int v1 = mesh.e2n_con[mesh.e2n_dsp[eidx] + 0],
  v2 = mesh.e2n_con[mesh.e2n_dsp[eidx] + 1],
  v3 = mesh.e2n_con[mesh.e2n_dsp[eidx] + 2],
  v4 = mesh.e2n_con[mesh.e2n_dsp[eidx] + 3];

  if     (v1 != b1 && v1 != b2) va = v1;
  else if(v2 != b1 && v2 != b2) va = v2;
  else if(v3 != b1 && v3 != b2) va = v3;
  else                          va = v4;

  if     (v1 != b1 && v1 != b2 && v1 != va) vb = v1;
  else if(v2 != b1 && v2 != b2 && v2 != va) vb = v2;
  else if(v3 != b1 && v3 != b2 && v3 != va) vb = v3;
  else                                      vb = v4;

  // generate connectivities
  mesh.e2n_con[con_start + 0] = va;
  mesh.e2n_con[con_start + 1] = b1;
  mesh.e2n_con[con_start + 2] = vb;
  mesh.e2n_con[con_start + 3] = nv;
  mesh.e2n_con[con_start + 4] = va;
  mesh.e2n_con[con_start + 5] = nv;
  mesh.e2n_con[con_start + 6] = vb;
  mesh.e2n_con[con_start + 7] = b2;

  // copy element data
  mesh.e2n_cnt[ele_start + 0] = 4;
  mesh.e2n_cnt[ele_start + 1] = 4;
  mesh.etype[ele_start + 0] = Tetra;
  mesh.etype[ele_start + 1] = Tetra;
  mesh.etags[ele_start + 0] = mesh.etags[eidx];
  mesh.etags[ele_start + 1] = mesh.etags[eidx];

  memcpy(mesh.lon.data() + ele_start*numfib,     mesh.lon.data() + eidx*numfib, numfib*sizeof(mt_real));
  memcpy(mesh.lon.data() + (ele_start+1)*numfib, mesh.lon.data() + eidx*numfib, numfib*sizeof(mt_real));
}

void mt_edge_splitter::split_iter(MT_MAP<mt_int, mixed_triple<mt_int,mt_int,float> > & split_edge,
    MT_USET<mt_int> & split_elem)
{
  mt_meshdata                   & mesh          = _edge_manager._mesh;
  mt_mapping<mt_int>            & e2g           = _edge_manager._ele2edge;
  MT_MAP<tuple<mt_int>, mt_int> & edge_map      = _edge_manager._edge_map;

  ++_inner_iter;
  short nlon = (mesh.lon.size() / 6) == mesh.e2n_cnt.size() ? 6 : 3;

  // prepare datastructs =============================================
  mt_int mesh_xyz_idx      = mesh.xyz.size() / 3;
  mt_int mesh_e2n_cnt_idx  = mesh.e2n_cnt.size();
  mt_int mesh_e2n_con_idx  = mesh.e2n_con.size();

  mt_int ele2edge_fwd_con_idx  = e2g.fwd_con.size();

  size_t new_nelem = 0, new_ncon = 0, new_nedge = 0;

  for(auto it = split_elem.begin(); it != split_elem.end(); ++it)
  {
    switch(mesh.etype[*it])
    {
      case Line:  new_nelem += 2; new_ncon += 4; new_nedge += 4;  break;
      case Tri:   new_nelem += 2; new_ncon += 6; new_nedge += 6;  break;
      case Tetra: new_nelem += 2; new_ncon += 8; new_nedge += 12; break;
      default:
                  std::cerr << "split_edges_iter: Error, element type not yet supported. Aborting!"
                    << std::endl;
                  exit(1);
    }
  }

  mesh_resize_elemdata(mesh, mesh.e2n_cnt.size() + new_nelem,
      mesh.e2n_con.size() + new_ncon);
  mesh.xyz.resize(mesh.xyz.size() + split_edge.size()*3);

  e2g.fwd_cnt.resize(mesh.e2n_cnt.size(), 0);
  e2g.fwd_con.resize(e2g.fwd_con.size() + new_nedge, 0);

  // start splitting =================================================
  for(auto it = split_edge.begin(); it != split_edge.end(); ++it)
  {
    // first split edge ----------------------------------------------
    tuple<mt_int> tupbuff;
    mt_int v1 = it->second.v1, v2 = it->second.v2, new_v = mesh_xyz_idx++;
    mt_int old_e = it->first;

    // add new vertex
    mt_point<mt_real> p1(mesh.xyz.data() + v1*3), p2(mesh.xyz.data() + v2*3);
    mt_point<mt_real> np = p1 + (p2 - p1) * mt_real(it->second.v3);
    mesh.xyz[new_v*3+0] = np.x;
    mesh.xyz[new_v*3+1] = np.y;
    mesh.xyz[new_v*3+2] = np.z;

    // old edge
    sortTuple(v1, v2, tupbuff.v1, tupbuff.v2);

    // now split elements connected to the edge ----------------------
    mt_int start = e2g.bwd_dsp[old_e], stop = start + e2g.bwd_cnt[old_e];
    for(mt_int i=start; i<stop; i++)
    {
      mt_int eidx = e2g.bwd_con[i];
      switch(mesh.etype[eidx])
      {
        case Line:
          split_line(eidx, mesh_e2n_cnt_idx, mesh_e2n_con_idx, tupbuff, new_v);

          e2g.fwd_cnt[mesh_e2n_cnt_idx+0] = 2;
          e2g.fwd_cnt[mesh_e2n_cnt_idx+1] = 2;

          add_edges_line(mesh.e2n_con.data() + mesh_e2n_con_idx + 0,
              e2g.fwd_con.data() + ele2edge_fwd_con_idx + 0,
              _edgeidx,
              edge_map);
          add_edges_line(mesh.e2n_con.data() + mesh_e2n_con_idx + 2,
              e2g.fwd_con.data() + ele2edge_fwd_con_idx + 2,
              _edgeidx,
              edge_map);

          mesh_e2n_cnt_idx += 2;
          mesh_e2n_con_idx += 4;
          ele2edge_fwd_con_idx += 4;
          break;

        case Tri:
          split_tri(eidx, mesh_e2n_cnt_idx, mesh_e2n_con_idx, tupbuff, new_v);

          e2g.fwd_cnt[mesh_e2n_cnt_idx+0] = 3;
          e2g.fwd_cnt[mesh_e2n_cnt_idx+1] = 3;

          add_edges_tri(mesh.e2n_con.data() + mesh_e2n_con_idx + 0,
              e2g.fwd_con.data() + ele2edge_fwd_con_idx + 0,
              _edgeidx,
              edge_map);
          add_edges_tri(mesh.e2n_con.data() + mesh_e2n_con_idx + 3,
              e2g.fwd_con.data() + ele2edge_fwd_con_idx + 3,
              _edgeidx,
              edge_map);

          mesh_e2n_cnt_idx += 2;
          mesh_e2n_con_idx += 6;
          ele2edge_fwd_con_idx += 6;
          break;

        case Tetra:
          split_tet(eidx, mesh_e2n_cnt_idx, mesh_e2n_con_idx, tupbuff, new_v);

          e2g.fwd_cnt[mesh_e2n_cnt_idx+0] = 6;
          e2g.fwd_cnt[mesh_e2n_cnt_idx+1] = 6;

          add_edges_tet(mesh.e2n_con.data() + mesh_e2n_con_idx + 0,
              e2g.fwd_con.data() + ele2edge_fwd_con_idx + 0,
              _edgeidx,
              edge_map);
          add_edges_tet(mesh.e2n_con.data() + mesh_e2n_con_idx + 4,
              e2g.fwd_con.data() + ele2edge_fwd_con_idx + 6,
              _edgeidx,
              edge_map);

          mesh_e2n_cnt_idx += 2;
          mesh_e2n_con_idx += 8;
          ele2edge_fwd_con_idx += 12;
          break;
        default:
          std::cerr << "split_edges_iter: Error, element type not yet supported. Aborting!"
            << std::endl;
          exit(1);
      }
    }
    edge_map.erase(tupbuff);

    _prg->next();
  }
  // remove split elements from mesh and ele2edge
  size_t widx=0;
  size_t e_ridx=0, e_widx=0, e2g_ridx=0, e2g_widx=0;
  size_t numelems = mesh.e2n_cnt.size();
  for(size_t i=0; i<numelems; i++)
  {
    if( split_elem.count(i) == 0 )
    {
      // mesh element data
      mesh.e2n_cnt[widx] = mesh.e2n_cnt[i];
      mesh.etags[widx] = mesh.etags[i];
      mesh.etype[widx] = mesh.etype[i];
      for(int j=0; j<nlon; j++)
        mesh.lon[widx*nlon + j] = mesh.lon[i*nlon + j];

      // ele2edge element data
      e2g.fwd_cnt[widx] = e2g.fwd_cnt[i];

      widx++;

      // mesh connectivity
      for(int j=0; j<mesh.e2n_cnt[i]; j++)
        mesh.e2n_con[e_widx++] = mesh.e2n_con[e_ridx++];

      // ele2edge connectivity
      for(int j=0; j<e2g.fwd_cnt[i]; j++)
        e2g.fwd_con[e2g_widx++] = e2g.fwd_con[e2g_ridx++];
    }
    else {
      e_ridx   += mesh.e2n_cnt[i];
      e2g_ridx += e2g.fwd_cnt[i];
    }
  }
  mesh.e2n_cnt.resize(widx);
  mesh.etags.resize(widx);
  mesh.etype.resize(widx);
  mesh.lon.resize(widx*nlon);
  mesh.e2n_con.resize(e_widx);

  e2g.fwd_cnt.resize(widx);
  e2g.fwd_con.resize(e2g_widx);
  e2g.transpose();
  e2g.setup_dsp();

  mesh.e2n_dsp.resize(mesh.e2n_cnt.size());
  bucket_sort_offset(mesh.e2n_cnt, mesh.e2n_dsp);
}

void mt_edge_splitter::operator()(std::set<edgeele, edge_len_dsc> & edge_set,
                                  PROGRESS<size_t> & inp_prg)
{
  _prg = & inp_prg;
  mt_mapping<mt_int> & e2g = _edge_manager._ele2edge;

  MT_MAP<mt_int, mixed_triple<mt_int,mt_int,float> > sel_edge;
  MT_USET<mt_int> split_elem;

  while(edge_set.size() > 0)
  {
    auto eit = edge_set.begin();
    split_elem.clear();
    sel_edge.clear();

    // iterate over edges in descending length order.
    while(eit != edge_set.end())
    {
      mt_int edge_idx = eit->idx;
      bool can_split = true;
      mt_int start = e2g.bwd_dsp[edge_idx], stop = start + e2g.bwd_cnt[edge_idx];
      for(mt_int i=start; i<stop; i++) {
        if(split_elem.count(e2g.bwd_con[i]) > 0)
        {
          can_split = false;
          break;
        }
      }

      if(can_split) {
        sel_edge[edge_idx] = {eit->v1, eit->v2,eit->split_at};
        for(mt_int i=start; i<stop; i++)
          split_elem.insert(e2g.bwd_con[i]);

        // remove edge from set since it was selected
        auto rem_it = eit++;
        edge_set.erase(rem_it);
      }
      else ++eit;
    }

    split_iter(sel_edge, split_elem);
  }
  _prg->finish();

  if( _inner_iter > RECOMP_THR ) {
    std::cout << "Recomputing node and edge indices.. " << std::endl;
    _inner_iter = 0;
    _edge_manager.refresh_edge_data(true);
    _edgeidx = _edge_manager.edges().size();
  }
}

void mt_edge_collapser::collapse_iter(MT_MAP<mt_int, mixed_tuple<mt_int, mt_point<mt_real> > > & nodmap,
    MT_USET<mt_int> & rem_elem)
{
  mt_meshdata & mesh = _edge_manager._mesh;

  short nlon = (mesh.lon.size() / 6) == mesh.e2n_cnt.size() ? 6 : 3;

  // loop over elements, and update connectivites
  for(size_t i=0; i<mesh.e2n_con.size(); i++) {
    mt_int c = mesh.e2n_con[i];
    if(nodmap.count(c)) {
      mt_int idx = nodmap[c].v1;
      mt_point<mt_real> p = nodmap[c].v2;

      mesh.e2n_con[i] = idx;
      mesh.xyz[idx*3+0] = p.x;
      mesh.xyz[idx*3+1] = p.y;
      mesh.xyz[idx*3+2] = p.z;
    }
  }

  // remove elements from mesh and e2g
  size_t widx=0;
  size_t e_ridx=0, e_widx=0;
  size_t numelems = mesh.e2n_cnt.size();
  for(size_t i=0; i<numelems; i++)
  {
    if( rem_elem.count(i) == 0 )
    {
      // mesh element data
      mesh.e2n_cnt[widx] = mesh.e2n_cnt[i];
      mesh.etags[widx] = mesh.etags[i];
      mesh.etype[widx] = mesh.etype[i];
      for(int j=0; j<nlon; j++)
        mesh.lon[widx*nlon + j] = mesh.lon[i*nlon + j];

      widx++;

      // mesh connectivity
      for(int j=0; j<mesh.e2n_cnt[i]; j++)
        mesh.e2n_con[e_widx++] = mesh.e2n_con[e_ridx++];

    }
    else {
      e_ridx   += mesh.e2n_cnt[i];
    }
  }
  mesh.e2n_cnt.resize(widx);
  mesh.etags.resize(widx);
  mesh.etype.resize(widx);
  mesh.lon.resize(widx*nlon);
  mesh.e2n_con.resize(e_widx);
  transpose_connectivity(mesh.e2n_cnt, mesh.e2n_con, mesh.n2e_cnt, mesh.n2e_con);

  mesh.e2n_dsp.resize(mesh.e2n_cnt.size());
  bucket_sort_offset(mesh.e2n_cnt, mesh.e2n_dsp);
  mesh.n2e_dsp.resize(mesh.n2e_cnt.size());
  bucket_sort_offset(mesh.n2e_cnt, mesh.n2e_dsp);
}

bool mt_edge_collapser::check_bad_tet(const MT_USET<mt_int> & edgeelems,
          const mt_vector<bool> & nodmask,
          const MT_MAP<mt_int, mixed_tuple<mt_int, mt_point<mt_real> > > & nodmap)
{
  mt_meshdata & mesh = _edge_manager._mesh;
  mt_int  con[4] = {0, 1, 2, 3};
  mt_real xyz[12];

  bool bad = false;

  for(const mt_int & eidx : edgeelems) {
    // we create an auxiliary element holding the remapped connectivity
    mt_int estart = mesh.e2n_dsp[eidx];
    for(size_t i=0; i<4; i++) {
      mt_int c = mesh.e2n_con[estart+i];
      if(nodmask[c]) {
        const mt_point<mt_real> & p = nodmap.find(c)->second.v2;
        xyz[i*3+0] = p.x;
        xyz[i*3+1] = p.y;
        xyz[i*3+2] = p.z;
      }
      else {
        xyz[i*3+0] = mesh.xyz[c*3+0];
        xyz[i*3+1] = mesh.xyz[c*3+1];
        xyz[i*3+2] = mesh.xyz[c*3+2];
      }
    }

    double qual = QMETRIC(con, xyz);
    if(has_bad_vol(Tetra, con, xyz, mt_real(0.2)) || qual > RESAMPLE_QUAL_THR) {
      bad = true;
      break;
    }
  }
  return bad;
}

bool mt_edge_collapser::check_bad_tri(const MT_USET<mt_int> & edgeelems,
          MT_MAP<mt_int, mixed_tuple<mt_int, mt_point<mt_real> > > & nodmap)
{
  mt_meshdata & mesh = _edge_manager._mesh;
  mt_int  con[3] = {0, 1, 2};
  mt_real xyz[9];
  bool bad = false;

  for(const mt_int & eidx : edgeelems) {
    // we create an auxiliary element holding the remapped connectivity
    mt_int estart = mesh.e2n_dsp[eidx];
    for(size_t i=0; i<3; i++) {
      mt_int c = mesh.e2n_con[estart+i];
      if(nodmap.count(c)) {
        const mt_point<mt_real> & p = nodmap[c].v2;
        xyz[i*3+0] = p.x;
        xyz[i*3+1] = p.y;
        xyz[i*3+2] = p.z;
      }
      else {
        xyz[i*3+0] = mesh.xyz[c*3+0];
        xyz[i*3+1] = mesh.xyz[c*3+1];
        xyz[i*3+2] = mesh.xyz[c*3+2];
      }
    }

    double qual = tri_qmetric_surface(con, xyz);
    if(has_bad_vol(Tri, con, xyz, mt_real(0.2)) || qual > RESAMPLE_QUAL_THR) {
      bad = true;
      break;
    }
  }

  return bad;
}

void mt_edge_collapser::operator()(std::set<edgeele, edge_len_asc> & edge_set,
                                   const mt_vector<bool> & surf_nodes,
                                   const mt_vector<mt_real> & surf_nrml,
                                   const mt_vector<bool> & fixed_surf_nodes,
                                   const mt_real surf_corr,
                                   PROGRESS<size_t> & inp_prg)
{
  _prg = & inp_prg;
  mt_mapping<mt_int> &  e2g = _edge_manager._ele2edge;
  mt_meshdata        & mesh = _edge_manager._mesh;

  MT_USET<tuple<mt_int> > sel_edge;
  MT_USET<mt_int> rem_elem;
  MT_MAP<mt_int, mixed_tuple<mt_int, mt_point<mt_real> > > nodmap;
  mt_vector<bool> nodmask(mesh.xyz.size() / 3, false);

  auto eit = edge_set.begin();
  size_t clp_cnt = 0;

  // iterate over edges in descending length order.
  while(eit != edge_set.end())
  {
    const mt_int v1 = eit->v1, v2 = eit->v2;

    bool mod_v1 = nodmask[v1], mod_v2 = nodmask[v2];
    bool sel = false, move_vtx = false, abort = false;

    mt_int rem_nod, keep_nod;

    if(!(mod_v1 || mod_v2))
    {
      bool v1_surf  = surf_nodes[v1],  v2_surf  = surf_nodes[v2];
      bool v1_fixed = fixed_surf_nodes[v1], v2_fixed = fixed_surf_nodes[v2];
      rem_nod = v1, keep_nod = v2;

      if(v1_surf) {
        if(v2_surf) {
          // both nodes are surface nodes
          mt_point<mt_real> n1(surf_nrml.data()+v1*3), n2(surf_nrml.data()+v2*3);
          if( (surf_corr == 0.0) || (n1.scaProd(n2) > surf_corr) ) {
            // both node normals correlate sufficiently
            if(v1_fixed || v2_fixed) {
              if(!v1_fixed) // v1 can be removed -> no vertex moving
                sel = true;
              if(!v2_fixed) { // v2 can be removed -> swap and no moving
                rem_nod = v2, keep_nod = v1;
                sel = true;
              }
            }
            else { // no vtx is fixed -> no problem
              sel = true;
              move_vtx = true;
            }
          }
        }
        else {
          // v2 is not on a surface -> we swap rem and keep nodes
          rem_nod = v2, keep_nod = v1;
          sel = true;
        }
      }
      else {
        sel = true;
        if(!v2_surf) move_vtx = true;
      }
    }

    if(sel) {
      mt_int edge_idx = eit->idx;
      MT_USET<mt_int> edgeelems, edge_rem;

      mixed_tuple<mt_int, mt_point<mt_real> > tb;
      tb.v1 = keep_nod;
      tb.v2 = mt_point<mt_real>(mesh.xyz.data() + keep_nod*3);
      nodmap[rem_nod] = tb;
      nodmask[rem_nod] = true;
      nodmap[keep_nod] = tb;
      nodmask[keep_nod] = true;

      mt_point<mt_real> & rem_pt = nodmap[rem_nod].v2;
      mt_point<mt_real> & keep_pt = nodmap[keep_nod].v2;

      // add elements connected to collapsed edge for removal
      mt_int start = e2g.bwd_dsp[edge_idx], stop = start + e2g.bwd_cnt[edge_idx];
      for(mt_int i=start; i<stop; i++)
        edge_rem.insert(e2g.bwd_con[i]);

      start = mesh.n2e_dsp[rem_nod], stop = start + mesh.n2e_cnt[rem_nod];
      for(mt_int i = start; i<stop; i++) {
        mt_int e = mesh.n2e_con[i];
        if(edge_rem.count(e) == 0) edgeelems.insert(e);
      }
      start = mesh.n2e_dsp[keep_nod], stop = start + mesh.n2e_cnt[keep_nod];
      for(mt_int i = start; i<stop; i++) {
        mt_int e = mesh.n2e_con[i];
        if(edge_rem.count(e) == 0) edgeelems.insert(e);
      }

      // TODO: factor this into function
      if(move_vtx) {
        // set keep_nod coord to center between rem_nod and keep_nod
        mt_point<mt_real> p1(mesh.xyz.data()+rem_nod*3), p2(mesh.xyz.data()+keep_nod*3);
        mt_point<mt_real> c = (p1 + p2) * mt_real(0.5);
        keep_pt = c, rem_pt = c;
        // check if collapsing to center yields bad elems
        bool bad = check_bad_tet(edgeelems, nodmask, nodmap);

        if(bad) {
          c = p1 + ((p2 - p1)*mt_real(0.25));
          keep_pt = c, rem_pt = c;
          bad = check_bad_tet(edgeelems, nodmask, nodmap);
        }
        if(bad) {
          c = p1 + ((p2 - p1)*mt_real(0.75));
          keep_pt = c, rem_pt = c;
          bad = check_bad_tet(edgeelems, nodmask, nodmap);
        }
        if(bad) abort = true;
      }
      else
        abort = check_bad_tet(edgeelems, nodmask, nodmap);

      if(abort) {
        nodmap.erase(keep_nod);
        nodmap.erase(rem_nod);
        nodmask[rem_nod] = false;
        nodmask[keep_nod] = false;
      }
      else {
        // the selected edge can be removed
        for(auto e : edge_rem) rem_elem.insert(e);
        auto rem_it = eit++;
        edge_set.erase(rem_it);
        clp_cnt++;
      }
    }

    if(abort || (!sel)) ++eit;
  }
  _prg->next();

  // apply node remapping to mesh
  collapse_iter(nodmap, rem_elem);
  _prg->next();

  // recompute edges
  _edge_manager.refresh_edge_data(false);
  _prg->next();

  // recompute surface normals on the fly
  #if 0
  for(size_t nidx=0; nidx < mesh.n2e_cnt.size(); nidx++) {
    if(mesh.n2e_cnt[nidx] > 0 && surf_nodes[nidx]) {
      MT_USET<mt_int> nset, eset, dset;
      nset.insert(nidx);
      nodeSet_to_elemSet(mesh, nset, eset);
      elemSet_to_nodeSet(mesh, eset, nset);
      nset.erase(nidx);

      for(auto n : nset)
        if(!surf_nodes[n]) dset.insert(n);
      for(auto n : dset) nset.erase(n);
    }
  }
  #endif
  _prg->finish();
}

void mt_edge_collapser::operator()(std::set<edgeele, edge_len_asc> & edge_set,
                                   const mt_vector<mt_real> & surf_nrml,
                                   const mt_vector<bool> & fixed_surf_nodes,
                                   const mt_real surf_corr,
                                   PROGRESS<size_t> & inp_prg)
{
  _prg = & inp_prg;
  mt_mapping<mt_int> &  e2g = _edge_manager._ele2edge;
  mt_meshdata        & mesh = _edge_manager._mesh;

  MT_USET<tuple<mt_int> > sel_edge;
  MT_USET<mt_int> rem_elem;
  MT_MAP<mt_int, mixed_tuple<mt_int, mt_point<mt_real> > > nodmap;

  auto eit = edge_set.begin();
  size_t clp_cnt = 0;

  // iterate over edges in descending length order.
  while(eit != edge_set.end())
  {
    const mt_int v1 = eit->v1, v2 = eit->v2;

    bool mod_v1  = nodmap.count(v1)  > 0, mod_v2  = nodmap.count(v2) > 0;
    bool sel = false, move_vtx = false, abort = false;

    mt_int rem_nod, keep_nod;

    if(!(mod_v1 || mod_v2))
    {
      bool v1_fixed = fixed_surf_nodes[v1], v2_fixed = fixed_surf_nodes[v2];
      rem_nod = v1, keep_nod = v2;

      // both nodes are surface nodes
      mt_point<mt_real> n1(surf_nrml.data()+v1*3), n2(surf_nrml.data()+v2*3);
      if( (surf_corr == 0.0) || (n1.scaProd(n2) > surf_corr) ) {
        // both node normals correlate sufficiently
        if(v1_fixed || v2_fixed) {
          if(!v1_fixed) // v1 can be removed -> no vertex moving
            sel = true;
          if(!v2_fixed) { // v2 can be removed -> swap and no moving
            rem_nod = v2, keep_nod = v1;
            sel = true;
          }
        }
        else { // no vtx is fixed -> no problem
          sel = true;
          move_vtx = true;
        }
      }
    }

    if(sel) {
      mt_int edge_idx = eit->idx;
      MT_USET<mt_int> edgeelems, edge_rem;

      mixed_tuple<mt_int, mt_point<mt_real> > tb;
      tb.v1 = keep_nod;
      tb.v2 = mt_point<mt_real>(mesh.xyz.data() + keep_nod*3);
      nodmap[rem_nod] = tb;
      nodmap[keep_nod] = tb;

      mt_point<mt_real> & rem_pt = nodmap[rem_nod].v2;
      mt_point<mt_real> & keep_pt = nodmap[keep_nod].v2;

      // add elements connected to collapsed edge for removal
      mt_int start = e2g.bwd_dsp[edge_idx], stop = start + e2g.bwd_cnt[edge_idx];
      for(mt_int i=start; i<stop; i++)
        edge_rem.insert(e2g.bwd_con[i]);

      start = mesh.n2e_dsp[rem_nod], stop = start + mesh.n2e_cnt[rem_nod];
      for(mt_int i = start; i<stop; i++) {
        mt_int e = mesh.n2e_con[i];
        if(edge_rem.count(e) == 0) edgeelems.insert(e);
      }
      start = mesh.n2e_dsp[keep_nod], stop = start + mesh.n2e_cnt[keep_nod];
      for(mt_int i = start; i<stop; i++) {
        mt_int e = mesh.n2e_con[i];
        if(edge_rem.count(e) == 0) edgeelems.insert(e);
      }

      // TODO: factor this into function
      if(move_vtx) {
        // set keep_nod coord to center between rem_nod and keep_nod
        mt_point<mt_real> p1(mesh.xyz.data()+rem_nod*3), p2(mesh.xyz.data()+keep_nod*3);
        mt_point<mt_real> c = (p1 + p2) * mt_real(0.5);
        keep_pt = c, rem_pt = c;
        // check if collapsing to center yields bad elems
        bool bad = check_bad_tri(edgeelems, nodmap);

        if(bad) {
          c = p1 + ((p2 - p1)*mt_real(0.25));
          keep_pt = c, rem_pt = c;
          bad = check_bad_tri(edgeelems, nodmap);
        }
        if(bad) {
          c = p1 + ((p2 - p1)*mt_real(0.75));
          keep_pt = c, rem_pt = c;
          bad = check_bad_tri(edgeelems, nodmap);
        }
        if(bad) {
          c = p1;
          keep_pt = c, rem_pt = c;
          bad = check_bad_tri(edgeelems, nodmap);
        }
        if(bad) {
          c = p2;
          keep_pt = c, rem_pt = c;
          bad = check_bad_tri(edgeelems, nodmap);
        }
        if(bad) abort = true;
      }
      else
        abort = check_bad_tri(edgeelems, nodmap);

      if(abort) {
        nodmap.erase(keep_nod);
        nodmap.erase(rem_nod);
      }
      else {
        // the selected edge can be removed
        for(auto e : edge_rem) rem_elem.insert(e);
        auto rem_it = eit++;
        edge_set.erase(rem_it);
        clp_cnt++;
      }
    }

    if(abort || (!sel)) ++eit;
  }
  _prg->next();

  // apply node remapping to mesh
  collapse_iter(nodmap, rem_elem);
  _prg->next();

  // recompute edges
  _edge_manager.refresh_edge_data(false);
  _prg->next(); _prg->finish();
}

