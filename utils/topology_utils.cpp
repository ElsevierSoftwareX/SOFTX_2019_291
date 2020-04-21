/**
* @file topology_utils.cpp
* @brief Topology manipulation utility functions.
* @author Aurel Neic
* @version
* @date 2017-08-16
*/

#include "mt_utils_base.h"
#include "io_utils.h"
#include "mesh_utils.h"
#include "mesh_smoothing.h"
#include "mesh_quality.h"
#include "topology_utils.h"
#include "dense_mat.hpp"


void elemSet_to_nodeSet(const mt_meshdata & mesh,
                        const MT_USET<mt_int> & elemSet,
                        MT_USET<mt_int> & nodeSet)
{
  MT_USET<mt_int>::const_iterator eit = elemSet.begin();
  while(eit != elemSet.end())
  {
    mt_int estart = mesh.e2n_dsp[*eit], estop = estart + mesh.e2n_cnt[*eit];
    for(mt_int i = estart; i<estop; i++)
      nodeSet.insert(mesh.e2n_con[i]);
    ++eit;
  }
}

void nodeSet_to_elemSet(const mt_meshdata & mesh,
                        const MT_USET<mt_int> & nodeSet,
                        MT_USET<mt_int> & elemSet)
{
  MT_USET<mt_int>::const_iterator it = nodeSet.begin();
  while(it != nodeSet.end())
  {
    mt_int start = mesh.n2e_dsp[*it], stop = start + mesh.n2e_cnt[*it];
    for(mt_int i = start; i<stop; i++)
      elemSet.insert(mesh.n2e_con[i]);
    ++it;
  }
}

// TODO: refactor to be able to specify what connectivity data to compute
void compute_full_mesh_connectivity(mt_meshdata & mesh, bool verbose)
{
  PROGRESS<short>* prog = NULL;

  if(verbose)
    prog = new PROGRESS<short>(3, "Computing connectivity graphs: ");

  mesh.e2n_dsp.resize(mesh.e2n_cnt.size());
  bucket_sort_offset(mesh.e2n_cnt, mesh.e2n_dsp);
  if(prog) prog->next();

  transpose_connectivity(mesh.e2n_cnt, mesh.e2n_con, mesh.n2e_cnt, mesh.n2e_con);
  if(prog) prog->next();
  multiply_connectivities(mesh.n2e_cnt, mesh.n2e_con, mesh.e2n_cnt, mesh.e2n_con, mesh.n2n_cnt, mesh.n2n_con);
  if(prog) prog->next();

  mesh.n2e_dsp.resize(mesh.n2e_cnt.size());
  mesh.n2n_dsp.resize(mesh.n2n_cnt.size());
  bucket_sort_offset(mesh.n2e_cnt, mesh.n2e_dsp);
  bucket_sort_offset(mesh.n2n_cnt, mesh.n2n_dsp);

  if(prog) {
    prog->finish();
    delete prog;
  }
}

void compute_full_mesh_connectivity(mt_meshdata & mesh, std::string basename, bool verbose)
{
  std::string filename = basename + FULL_CON_EXT;
  bool have_read_con = false;

  if(file_exists(filename)) {
    read_full_mesh_connectivity(mesh, filename);
    have_read_con = true;

    // do checks
    if(mesh.n2e_cnt.size() != (mesh.xyz.size() / 3))
      have_read_con = false;

    if(mesh.n2n_cnt.size() != (mesh.xyz.size() / 3))
      have_read_con = false;

    mt_int e2n_edges = bucket_sort_size(mesh.e2n_cnt);
    mt_int n2e_edges = bucket_sort_size(mesh.n2e_cnt);

    if(e2n_edges != n2e_edges)
      have_read_con = false;

    // notify user that data is inconsistent
    if(have_read_con == false)
      fprintf(stderr, "%s warning: connectivity read from %s is inconsistent! Recomputing it!\n",
              __func__, filename.c_str());
  }

  if(have_read_con == false) {
    compute_full_mesh_connectivity(mesh, verbose);
    write_full_mesh_connectivity(mesh, filename);
  }
}

void elements_with_edge(const mt_meshdata & mesh, mt_int v0, mt_int v1, std::set<mt_int> & elemset)
{
  mt_int estart = mesh.n2e_dsp[v0], estop = estart + mesh.n2e_cnt[v0];
  for(mt_int i=estart; i<estop; i++)
  {
    mt_int eidx = mesh.n2e_con[i];
    mt_int nstart = mesh.e2n_dsp[eidx], nstop = nstart + mesh.e2n_cnt[eidx];
    for(mt_int j=nstart; j<nstop; j++)
    {
      mt_int nidx = mesh.e2n_con[j];
      if(nidx == v1) {
        elemset.insert(eidx);
        break;
      }
    }
  }
}

void elements_with_face(const mt_meshdata & mesh, mt_int v0, mt_int v1, mt_int v2,
                        std::set<mt_int> & elemset)
{
  // iterate over all elements connected to v0 and look for (v1, v2)
  mt_int estart = mesh.n2e_dsp[v0], estop = estart + mesh.n2e_cnt[v0];
  for(mt_int i=estart; i<estop; i++)
  {
    mt_int eidx = mesh.n2e_con[i];
    mt_int nstart = mesh.e2n_dsp[eidx], nsize = mesh.e2n_cnt[eidx];
    const mt_int* con = mesh.e2n_con.data() + nstart;
    short pos1 = -1, pos2 = -1;
    for(mt_int j=0; j<nsize; j++)
    {
      mt_int nidx = con[j];
      if(nidx == v1) pos1 = j;
      else if(nidx == v2) pos2 = j;
    }
    // (v1, v2) found -> insert eidx
    if( (pos1 > -1) && (pos2 > -1) )
      elemset.insert(eidx);
  }
}

void tet_refine_centerPoint(mt_meshdata & mesh, const std::set<mt_int> & ref_elem)
{
  size_t old_nelem = mesh.e2n_cnt.size();
  size_t old_npts  = mesh.xyz.size() / 3;

  short lonsize = mesh.lon.size() == old_nelem * 6 ? 6 : 3;

  size_t nelem = old_nelem + ref_elem.size() * 4;
  mesh_resize_elemdata(mesh, nelem, nelem*4);
  mesh.xyz.resize(old_npts*3 + ref_elem.size()*3);

  PROGRESS<size_t> prog(ref_elem.size(), "Center-Point refining progress: ");

  size_t write_pidx = old_npts;
  size_t write_eidx = old_nelem;

  for(auto it = ref_elem.begin(); it != ref_elem.end(); ++it)
  {
    prog.next();
    mt_int eidx = *it;
    mt_int v0 = mesh.e2n_con[eidx*4+0];
    mt_int v1 = mesh.e2n_con[eidx*4+1];
    mt_int v2 = mesh.e2n_con[eidx*4+2];
    mt_int v3 = mesh.e2n_con[eidx*4+3];

    // computing center point
    mt_point<mt_real> cp(0, 0, 0);
    cp += mt_point<mt_real>(mesh.xyz.data()+v0*3);
    cp += mt_point<mt_real>(mesh.xyz.data()+v1*3);
    cp += mt_point<mt_real>(mesh.xyz.data()+v2*3);
    cp += mt_point<mt_real>(mesh.xyz.data()+v3*3);
    cp /= 4.0;

    // write new point
    mesh.xyz[write_pidx*3+0] = cp.x;
    mesh.xyz[write_pidx*3+1] = cp.y;
    mesh.xyz[write_pidx*3+2] = cp.z;

    // write new elems
    mt_real* read_lon = mesh.lon.data() + eidx*lonsize;

    mesh.e2n_cnt[write_eidx] = 4;
    mesh.etags  [write_eidx] = mesh.etags[eidx];
    mesh.etype  [write_eidx] = Tetra;
    mesh.e2n_con[write_eidx*4+0] = v0;
    mesh.e2n_con[write_eidx*4+1] = v1;
    mesh.e2n_con[write_eidx*4+2] = v2;
    mesh.e2n_con[write_eidx*4+3] = write_pidx;
    memcpy(mesh.lon.data() + write_eidx*lonsize, read_lon, lonsize*sizeof(mt_real));
    write_eidx++;

    mesh.e2n_cnt[write_eidx] = 4;
    mesh.etags  [write_eidx] = mesh.etags[eidx];
    mesh.etype  [write_eidx] = Tetra;
    mesh.e2n_con[write_eidx*4+0] = v1;
    mesh.e2n_con[write_eidx*4+1] = v3;
    mesh.e2n_con[write_eidx*4+2] = v2;
    mesh.e2n_con[write_eidx*4+3] = write_pidx;
    memcpy(mesh.lon.data() + write_eidx*lonsize, read_lon, lonsize*sizeof(mt_real));
    write_eidx++;

    mesh.e2n_cnt[write_eidx] = 4;
    mesh.etags  [write_eidx] = mesh.etags[eidx];
    mesh.etype  [write_eidx] = Tetra;
    mesh.e2n_con[write_eidx*4+0] = v0;
    mesh.e2n_con[write_eidx*4+1] = v2;
    mesh.e2n_con[write_eidx*4+2] = v3;
    mesh.e2n_con[write_eidx*4+3] = write_pidx;
    memcpy(mesh.lon.data() + write_eidx*lonsize, read_lon, lonsize*sizeof(mt_real));
    write_eidx++;

    mesh.e2n_cnt[write_eidx] = 4;
    mesh.etags  [write_eidx] = mesh.etags[eidx];
    mesh.etype  [write_eidx] = Tetra;
    mesh.e2n_con[write_eidx*4+0] = v0;
    mesh.e2n_con[write_eidx*4+1] = v3;
    mesh.e2n_con[write_eidx*4+2] = v1;
    mesh.e2n_con[write_eidx*4+3] = write_pidx;
    memcpy(mesh.lon.data() + write_eidx*lonsize, read_lon, lonsize*sizeof(mt_real));
    write_eidx++;

    write_pidx++;
  }
  prog.finish();

  mt_vector<bool> keep(mesh.e2n_cnt.size(), true);

  for(auto it = ref_elem.begin(); it != ref_elem.end(); ++it)
    keep[*it] = false;

  mt_vector<mt_int> nod, eidx;
  restrict_meshdata(keep, mesh, nod, eidx);
}


void refine_uniform(mt_meshdata & mesh,
                    const MT_MAP<tuple<mt_int>, mt_int> & edge_map)
{
  mt_meshdata outmesh;

  size_t num_elem     = mesh.e2n_cnt.size();
  size_t num_elem_out = 0;
  int numfib = mesh.lon.size() == num_elem* 6 ? 6 : 3;

  for(size_t i=0; i<num_elem; i++) {
    if(mesh.etype[i] == Tri)        num_elem_out += 4;
    else if(mesh.etype[i] == Tetra) num_elem_out += 8;
    else {
      check_elem_type(mesh.etype[i], Tetra, "refine_uniform");
    }
  }

  outmesh.e2n_cnt.resize(num_elem_out);
  outmesh.etype.resize(num_elem_out);
  outmesh.etags.resize(num_elem_out);
  outmesh.lon.resize(num_elem_out*numfib);

  size_t widx = 0;
  // copy fibers and tags
  for(size_t i=0; i<num_elem; i++) {
    if(mesh.etype[i] == Tri) {
      for(int j=0; j<4; j++) {
        outmesh.etype[widx+j]   = Tri;
        outmesh.etags[widx+j]   = mesh.etags[i];
        outmesh.e2n_cnt[widx+j] = 3;
        memcpy(outmesh.lon.data()+(widx+j)*numfib, mesh.lon.data()+i*numfib, numfib*sizeof(mt_real));
      }
      widx += 4;
    }
    else {
      for(int j=0; j<8; j++) {
        outmesh.etype[widx+j]   = Tetra;
        outmesh.etags[widx+j]   = mesh.etags[i];
        outmesh.e2n_cnt[widx+j] = 4;
        memcpy(outmesh.lon.data()+(widx+j)*numfib, mesh.lon.data()+i*numfib, numfib*sizeof(mt_real));
      }
      widx += 8;
    }
  }

  // the edges are also indexing our additional vertices
  size_t num_pts     = mesh.xyz.size() / 3;
  size_t num_pts_out = num_pts + edge_map.size();
  outmesh.xyz.resize(num_pts_out*3);

  mt_real* xyz = mesh.xyz.data(), *outxyz = outmesh.xyz.data();
  // copy coordinates
  {
    memcpy(outxyz, xyz, mesh.xyz.size()*sizeof(mt_real));
    outxyz += mesh.xyz.size();

    for(auto it : edge_map) {
      vec3r p1(xyz + it.first.v1*3);
      vec3r p2(xyz + it.first.v2*3);
      vec3r np = (p1 + p2) * 0.5;

      mt_int edge_idx = it.second;
      outxyz[edge_idx*3+0] = np.x;
      outxyz[edge_idx*3+1] = np.y;
      outxyz[edge_idx*3+2] = np.z;
    }
  }

  size_t outcon_size = bucket_sort_size(outmesh.e2n_cnt);
  outmesh.e2n_con.resize(outcon_size);

  mt_int *con = mesh.e2n_con.data(), *outcon = outmesh.e2n_con.data();
  tuple<mt_int> buff[64]; // unneeded buff to fit API

  for(size_t i=0; i<num_elem; i++) {
    if(mesh.etype[i] == Tri) {
      mt_int edges[3], v1 = con[0], v2 = con[1], v3 = con[2];
      get_edges_tri(con, edge_map, edges, buff);
      con += 3;

      // a=(1,2) b=(2,3) c=(1,3)
      mt_int a = num_pts + edges[0], b = num_pts + edges[1], c = num_pts + edges[2];

      outcon[0] =v1,  outcon[1] = a, outcon[2] = c; outcon += 3;
      outcon[0] = a,  outcon[1] = b, outcon[2] = c; outcon += 3;
      outcon[0] = a,  outcon[1] =v2, outcon[2] = b; outcon += 3;
      outcon[0] = c,  outcon[1] = b, outcon[2] =v3; outcon += 3;
    }
    else {
      mt_int edges[6], v1 = con[0], v2 = con[1], v3 = con[2], v4 = con[3];
      get_edges_tet(con, edge_map, edges, buff);

      con += 4;
      // a=(1,2) b=(1,4) c=(2,4) d=(2,3) e=(3,4) f=(1,3)
      mt_int a = num_pts + edges[0], b = num_pts + edges[3], c = num_pts + edges[4],
             d = num_pts + edges[1], e = num_pts + edges[5], f = num_pts + edges[2];

      outcon[0] =v1,  outcon[1] = a, outcon[2] = f, outcon[3] = b; outcon += 4;
      outcon[0] = a,  outcon[1] =v2, outcon[2] = d, outcon[3] = c; outcon += 4;
      outcon[0] = d,  outcon[1] =v3, outcon[2] = f, outcon[3] = e; outcon += 4;
      outcon[0] = a,  outcon[1] = d, outcon[2] = f, outcon[3] = c; outcon += 4;
      outcon[0] = c,  outcon[1] = d, outcon[2] = e, outcon[3] = f; outcon += 4;
      outcon[0] = a,  outcon[1] = c, outcon[2] = b, outcon[3] = f; outcon += 4;
      outcon[0] = f,  outcon[1] = e, outcon[2] = b, outcon[3] = c; outcon += 4;
      outcon[0] = b,  outcon[1] = c, outcon[2] = e, outcon[3] =v4; outcon += 4;
    }
  }

  mesh = outmesh;
}


void add_edge(const mt_int v1,
              const mt_int v2,
              tuple<mt_int> & tupbuff,
              size_t & edgeidx,
              mt_int* ele2edge_con,
              MT_MAP<tuple<mt_int>, mt_int> & edge_map)
{
  sortTuple(v1, v2, tupbuff.v1, tupbuff.v2);
  auto it = edge_map.find(tupbuff);

  if(it != edge_map.end())
    *ele2edge_con = it->second;
  else {
    edge_map[tupbuff] = mt_int(edgeidx);
    *ele2edge_con = edgeidx;
    edgeidx++;
  }
}

void add_face(const mt_int v1,
              const mt_int v2,
              const mt_int v3,
              triple<mt_int> & tripbuff,
              size_t & faceidx,
              mt_int* ele2face_con,
              MT_MAP<triple<mt_int>, mt_int> & face_map)
{
  sortTriple(v1, v2, v3, tripbuff.v1, tripbuff.v2, tripbuff.v3);
  auto it = face_map.find(tripbuff);

  if(it != face_map.end())
    *ele2face_con = it->second;
  else {
    face_map[tripbuff] = mt_int(faceidx);
    *ele2face_con = faceidx;
    faceidx++;
  }
}

mt_int get_edge(const mt_int v1,
                const mt_int v2,
                tuple<mt_int> & tupbuff,
                const MT_MAP<tuple<mt_int>, mt_int> & edge_map)
{
  sortTuple(v1, v2, tupbuff.v1, tupbuff.v2);
  auto it = edge_map.find(tupbuff);

  if(it != edge_map.end()) return it->second;
  else                     return -1;
}

void add_edges_line(const mt_int* mesh_con,
                    mt_int* ele2edge_con,
                    size_t & edgeidx,
                    MT_MAP<tuple<mt_int>, mt_int> & edge_map)
{
  struct tuple<mt_int> tupbuff;

  mt_int v1 = mesh_con[0];
  mt_int v2 = mesh_con[1];

  add_edge(v1, v2, tupbuff, edgeidx, ele2edge_con+0, edge_map);
}

void add_edges_tri(const mt_int* mesh_con,
                   mt_int* ele2edge_con,
                   size_t & edgeidx,
                   MT_MAP<tuple<mt_int>, mt_int> & edge_map)
{
  struct tuple<mt_int> tupbuff;

  mt_int v1 = mesh_con[0];
  mt_int v2 = mesh_con[1];
  mt_int v3 = mesh_con[2];

  // (1,2)
  add_edge(v1, v2, tupbuff, edgeidx, ele2edge_con+0, edge_map);
  // (2,3)
  add_edge(v2, v3, tupbuff, edgeidx, ele2edge_con+1, edge_map);
  // (3,1)
  add_edge(v3, v1, tupbuff, edgeidx, ele2edge_con+2, edge_map);
}

void add_edges_tet(const mt_int* mesh_con,
                   mt_int* ele2edge_con,
                   size_t & edgeidx,
                   MT_MAP<tuple<mt_int>, mt_int> & edge_map)
{
  struct tuple<mt_int> tupbuff;

  mt_int v1 = mesh_con[0];
  mt_int v2 = mesh_con[1];
  mt_int v3 = mesh_con[2];
  mt_int v4 = mesh_con[3];

  // (1,2)
  add_edge(v1, v2, tupbuff, edgeidx, ele2edge_con+0, edge_map);
  // (2,3)
  add_edge(v2, v3, tupbuff, edgeidx, ele2edge_con+1, edge_map);
  // (3,1)
  add_edge(v3, v1, tupbuff, edgeidx, ele2edge_con+2, edge_map);
  // (1,4)
  add_edge(v1, v4, tupbuff, edgeidx, ele2edge_con+3, edge_map);
  // (2,4)
  add_edge(v2, v4, tupbuff, edgeidx, ele2edge_con+4, edge_map);
  // (3,4)
  add_edge(v3, v4, tupbuff, edgeidx, ele2edge_con+5, edge_map);
}

void get_edges_tri(const mt_int* mesh_con,
                   const MT_MAP<tuple<mt_int>, mt_int> & edge_map,
                   mt_int* edge_idx, tuple<mt_int>* edge)
{
  mt_int v1 = mesh_con[0];
  mt_int v2 = mesh_con[1];
  mt_int v3 = mesh_con[2];

  // (1,2)
  edge_idx[0] = get_edge(v1, v2, edge[0], edge_map);
  // (2,3)
  edge_idx[1] = get_edge(v2, v3, edge[1], edge_map);
  // (3,1)
  edge_idx[2] = get_edge(v3, v1, edge[2], edge_map);
}

void get_edges_quad(const mt_int* mesh_con,
                    const MT_MAP<tuple<mt_int>, mt_int> & edge_map,
                    mt_int* edge_idx, tuple<mt_int>* edge)
{
  mt_int v1 = mesh_con[0];
  mt_int v2 = mesh_con[1];
  mt_int v3 = mesh_con[2];
  mt_int v4 = mesh_con[2];

  // (1,2)
  edge_idx[0] = get_edge(v1, v2, edge[0], edge_map);
  // (2,3)
  edge_idx[1] = get_edge(v2, v3, edge[1], edge_map);
  // (3,4)
  edge_idx[2] = get_edge(v3, v4, edge[2], edge_map);
  // (4,1)
  edge_idx[3] = get_edge(v4, v1, edge[3], edge_map);
}

void get_edges_tet(const mt_int* mesh_con,
                   const MT_MAP<tuple<mt_int>, mt_int> & edge_map,
                   mt_int* edge_idx, tuple<mt_int>* edge)
{
  mt_int v1 = mesh_con[0];
  mt_int v2 = mesh_con[1];
  mt_int v3 = mesh_con[2];
  mt_int v4 = mesh_con[3];

  // (1,2)
  edge_idx[0] = get_edge(v1, v2, edge[0], edge_map);
  // (2,3)
  edge_idx[1] = get_edge(v2, v3, edge[1], edge_map);
  // (3,1)
  edge_idx[2] = get_edge(v3, v1, edge[2], edge_map);
  // (1,4)
  edge_idx[3] = get_edge(v1, v4, edge[3], edge_map);
  // (2,4)
  edge_idx[4] = get_edge(v2, v4, edge[4], edge_map);
  // (3,4)
  edge_idx[5] = get_edge(v3, v4, edge[5], edge_map);
}

int get_edges(const mt_meshdata & mesh, const MT_MAP<tuple<mt_int>, mt_int> & edge_map,
              const mt_int eidx, mt_vector<mt_int> & edge_indices, mt_vector<tuple<mt_int>> & edges)
{
  const mt_int* con = mesh.e2n_con.data() + mesh.e2n_dsp[eidx];

  size_t oldsize = edges.size();
  int numedges;

  switch(mesh.etype[eidx]) {
    case Tri:
      numedges = 3;
      edge_indices.resize(oldsize + numedges);
      edges.resize(oldsize + numedges);
      get_edges_tri (con, edge_map, edge_indices.data() + oldsize, edges.data() + oldsize);
      break;

    case Quad:
      numedges = 4;
      edge_indices.resize(oldsize + numedges);
      edges.resize(oldsize + numedges);
      get_edges_quad(con, edge_map, edge_indices.data() + oldsize, edges.data() + oldsize);
      break;

    case Tetra:
      numedges = 6;
      edge_indices.resize(oldsize + numedges);
      edges.resize(oldsize + numedges);
      get_edges_tet (con, edge_map, edge_indices.data() + oldsize, edges.data() + oldsize);
      break;

    default:
      fprintf(stderr, "%s error: Unimplemented element type! Aborting!\n", __func__);
      exit(1);
  }

  return numedges;
}

void add_faces_tri(const mt_int* mesh_con,
                   mt_int* ele2face_con,
                   size_t & faceidx,
                   MT_MAP<triple<mt_int>, mt_int> & face_map)
{
  struct triple<mt_int> tripbuff;

  mt_int v1 = mesh_con[0];
  mt_int v2 = mesh_con[1];
  mt_int v3 = mesh_con[2];

  add_face(v1, v2, v3, tripbuff, faceidx, ele2face_con, face_map);
}

void add_faces_tet(const mt_int* mesh_con,
                   mt_int* ele2face_con,
                   size_t & faceidx,
                   MT_MAP<triple<mt_int>, mt_int> & face_map)
{
  struct triple<mt_int> tripbuff;

  mt_int v1 = mesh_con[0];
  mt_int v2 = mesh_con[1];
  mt_int v3 = mesh_con[2];
  mt_int v4 = mesh_con[3];

  // (1,2,3)
  add_face(v1, v2, v3, tripbuff, faceidx, ele2face_con+0, face_map);
  // (2,3,4)
  add_face(v2, v3, v4, tripbuff, faceidx, ele2face_con+1, face_map);
  // (1,3,4)
  add_face(v1, v3, v4, tripbuff, faceidx, ele2face_con+2, face_map);
  // (1,2,4)
  add_face(v1, v2, v4, tripbuff, faceidx, ele2face_con+3, face_map);
}

void compute_edges(const mt_meshdata & mesh,
                   mt_mapping<mt_int> & ele2edge,
                   MT_MAP<tuple<mt_int>, mt_int> & edge_map)
{
  assert(mesh.etype.size() == mesh.e2n_cnt.size());

  ele2edge.fwd_cnt.assign(mesh.e2n_cnt.size(), mt_int(0));
  ele2edge.fwd_dsp.assign(mesh.e2n_cnt.size());
  size_t numedges = 0;

  for(size_t i=0; i<mesh.e2n_cnt.size(); i++)
  {
    mt_int n = -1;
    switch(mesh.etype[i])
    {
      case Line: n = 1;  break;
      case Tri: n = 3;   break;
      case Tetra: n = 6; break;
      default:
        std::cerr << "compute_edges: Error, element type not yet supported. Aborting!"
                  << std::endl;
        exit(1);
    }
    numedges += n;
    ele2edge.fwd_cnt[i] = n;
  }
  bucket_sort_offset(ele2edge.fwd_cnt, ele2edge.fwd_dsp);
  ele2edge.fwd_con.resize(numedges);
  size_t edgeidx=0;

  // loop over elements and construct edge_map and ele2edge map in one go
  for(size_t i=0; i<mesh.e2n_cnt.size(); i++)
  {
    const mt_int* mesh_con = mesh.e2n_con.data() + mesh.e2n_dsp[i];
    mt_int* ele2edge_con = ele2edge.fwd_con.data() + ele2edge.fwd_dsp[i];

    switch(mesh.etype[i])
    {
      case Line:
        add_edges_line(mesh_con, ele2edge_con, edgeidx, edge_map);
        break;
      case Tri:
        add_edges_tri(mesh_con, ele2edge_con, edgeidx, edge_map);
        break;
      case Tetra:
        add_edges_tet(mesh_con, ele2edge_con, edgeidx, edge_map);
        break;
      default:
        std::cerr << "compute_edges: Error, element type not yet supported. Aborting!"
                  << std::endl;
        exit(1);
    }
  }
}

void compute_faces(const mt_meshdata & mesh,
                   mt_mapping<mt_int> & ele2face,
                   MT_MAP<triple<mt_int>, mt_int> & face_map)
{
  assert(mesh.etype.size() == mesh.e2n_cnt.size());

  ele2face.fwd_cnt.assign(mesh.e2n_cnt.size(), mt_int(0));
  ele2face.fwd_dsp.assign(mesh.e2n_cnt.size());
  size_t numfaces = 0;

  for(size_t i=0; i<mesh.e2n_cnt.size(); i++)
  {
    mt_int n = -1;
    switch(mesh.etype[i])
    {
      case Tri:   n = 1; break;
      case Tetra: n = 4; break;
      default:
        std::cerr << "compute_faces: Error, element type not yet supported. Aborting!"
                  << std::endl;
        exit(1);
    }
    numfaces += n;
    ele2face.fwd_cnt[i] = n;
  }
  bucket_sort_offset(ele2face.fwd_cnt, ele2face.fwd_dsp);
  ele2face.fwd_con.resize(numfaces);
  size_t faceidx=0;

  // loop over elements and construct face_map and ele2face map in one go
  for(size_t i=0; i<mesh.e2n_cnt.size(); i++)
  {
    const mt_int* mesh_con = mesh.e2n_con.data() + mesh.e2n_dsp[i];
    mt_int* ele2face_con = ele2face.fwd_con.data() + ele2face.fwd_dsp[i];

    switch(mesh.etype[i])
    {
      case Tri:
        add_faces_tri(mesh_con, ele2face_con, faceidx, face_map);
        break;

      case Tetra:
        add_faces_tet(mesh_con, ele2face_con, faceidx, face_map);
        break;

      default:
        std::cerr << "compute_faces: Error, element type not yet supported. Aborting!"
                  << std::endl;
        exit(1);
    }
  }
}


void edgemesh_from_edgemap(const MT_MAP<tuple<mt_int>, mt_int> & edgemap,
                           mt_meshdata & edgemesh)
{
  edgemesh.etype.assign(edgemap.size(), Line);
  edgemesh.e2n_cnt.assign(edgemap.size(), 2);
  edgemesh.e2n_con.resize(edgemap.size()*2);

  size_t edge_idx=0;
  for(auto it = edgemap.begin(); it != edgemap.end(); ++it)
  {
    edgemesh.e2n_con[edge_idx*2+0] = it->first.v1;
    edgemesh.e2n_con[edge_idx*2+1] = it->first.v2;
    edge_idx++;
  }
}

void compute_line_interfaces(const mt_meshdata & mesh,
                             const mt_meshdata & surfmesh,
                             const mt_real edge_ang,
                             const bool edges_on_geom_surf,
                             mt_meshdata & linemesh)
{
  MT_MAP<tuple<mt_int>, mt_int> sel_edge_map; // edges we are interested in

  // add sharp edges to selection
  if(edge_ang > 0) {
    MT_USET<mt_int> sharp_vtx;
    MT_MAP<tuple<mt_int>, mt_int> edge_map;

    if(edges_on_geom_surf) {
      // compute surface of whole geometry
      mt_meshdata geo_surf;
      compute_surface(mesh.etype, mesh.e2n_cnt, mesh.e2n_con, geo_surf, true);
      compute_full_mesh_connectivity(geo_surf);

      // identify sharp edges on geomery surface
      mt_vector<mt_int> vtx = geo_surf.e2n_con;
      binary_sort(vtx); unique_resize(vtx);
      identify_sharp_edge_nodes(geo_surf, mesh.xyz, vtx, sharp_vtx, edge_ang);

      mt_mapping<mt_int> ele2edge;
      compute_edges(geo_surf, ele2edge, edge_map);
    }
    else {
      mt_vector<mt_int> vtx = surfmesh.e2n_con;
      binary_sort(vtx); unique_resize(vtx);
      identify_sharp_edge_nodes(surfmesh, mesh.xyz, vtx, sharp_vtx, edge_ang);

      mt_mapping<mt_int> ele2edge;
      compute_edges(surfmesh, ele2edge, edge_map);
    }

    for(auto it = edge_map.begin(); it != edge_map.end(); ++it)
      if(sharp_vtx.count(it->first.v1) > 0 && sharp_vtx.count(it->first.v2) > 0)
        sel_edge_map.insert(*it);
  }

  // add surface interfaces
  {
    // compute edges of surface
    mt_mapping<mt_int> ele2edge;
    MT_MAP<tuple<mt_int>, mt_int> edge_map;     // all edges of the surface mesh

    compute_edges(surfmesh, ele2edge, edge_map);
    ele2edge.transpose();

    // add interfaces between surfaces to selection
    for(auto it = edge_map.begin(); it != edge_map.end(); ++it){
      mt_int edge_idx = it->second;
      if(ele2edge.bwd_cnt[edge_idx] > 2)
        sel_edge_map.insert(*it);
    }
  }
  edgemesh_from_edgemap(sel_edge_map, linemesh);

  #if 0
  linemesh.xyz = mesh.xyz;
  linemesh.etags.assign(linemesh.e2n_cnt.size(), 1);
  linemesh.etype.assign(linemesh.e2n_cnt.size(), Line);
  write_mesh_selected(linemesh, "vtk", "debug_linemesh");
  #endif
}

void linear_search_vtx(const mt_vector<mt_real> & xyz,
                       const mt_point<mt_real> & ref_pt,
                       mt_int & idx, mt_real & dist)
{
  size_t nnodes = xyz.size() / 3;

  mt_point<mt_real> e;
  mt_real      mindist = INFINITY;
  size_t minidx  = 0;

  for(size_t i=0; i < nnodes; i++)
  {
    e = ref_pt - mt_point<mt_real>(xyz.data() + i*3);
    mt_real d = e.length2();

    if(mindist > d) {
      mindist = d;
      minidx = i;
    }
  }

  idx  = minidx;
  dist = mindist;
}

void linear_search_vtx(const mt_vector<mt_real> & xyz,
                       const MT_USET<mt_int> & nod_set,
                       const mt_point<mt_real> & ref_pt,
                       mt_int & idx, mt_real & dist)
{
  mt_point<mt_real> e;
  mt_real      mindist = INFINITY;
  size_t minidx  = 0;

  for(auto it = nod_set.begin(); it != nod_set.end(); ++it)
  {
    e = ref_pt - mt_point<mt_real>(xyz.data() + (*it)*3);
    mt_real d = e.length2();

    if(mindist > d) {
      mindist = d;
      minidx = *it;
    }
  }

  idx  = minidx;
  dist = mindist;
}

void neighbour_closest_to_vtx(const mt_meshdata & mesh,
                              const mt_point<mt_real> & ref_pt,
                              const mt_int nidx,
                              mt_int & end_idx,
                              mt_real & end_dist)
{
  // we create a set of (distance, vtx) tuples of all neighbours of the nidx vertex
  // thus, the closest vertex can be retreived trivially

  mt_int nstart = mesh.n2n_dsp[nidx], nstop = nstart + mesh.n2n_cnt[nidx];
  std::set<mixed_tuple<mt_real, mt_int> > neighbours; // the neighbours sorted by distance
  for(mt_int i=nstart; i<nstop; i++)
  {
    mt_int tidx = mesh.n2n_con[i];
    if( tidx != nidx) {
      mt_point<mt_real> e = mt_point<mt_real>(mesh.xyz.data() + tidx*3) - ref_pt;
      neighbours.insert({e.length2(), tidx});
    }
  }

  // the closest vertex is sorted first
  if(neighbours.size() > 0)
  {
    typename std::set<mixed_tuple<mt_real, mt_int> >::iterator closest = neighbours.begin();
    end_dist = closest->v1;
    end_idx  = closest->v2;
  }
  else {
    std::cerr << "Warning: Node " << nidx << " has no neighbours!" << std::endl;
    linear_search_vtx(mesh.xyz, ref_pt, end_idx, end_dist);
  }
}

void find_closest_vtx(const mt_meshdata & mesh,
                      const mt_point<mt_real> & ref_pt,
                      const mt_int start_idx,
                      mt_int & end_idx,
                      mt_real & end_dist)
{
  mt_int last_idx  = start_idx, curr_idx = start_idx;
  mt_real last_dist = INFINITY, curr_dist = INFINITY;

  // we iterate over the neighbourhoods of vertices until we have the one with the
  // closest distance.
  do {
    last_dist = curr_dist;
    last_idx  = curr_idx;
    neighbour_closest_to_vtx(mesh, ref_pt, last_idx, curr_idx, curr_dist);
  } while (curr_dist < last_dist);

  end_idx  = last_idx;
  end_dist = last_dist;
}

void compute_correspondance(const mt_meshdata & mesh1,
                            const mt_meshdata & mesh2,
                            mt_vector<mt_int> & corr,
                            mt_vector<mt_real> & corr_dist)
{
  size_t nnodes = mesh1.xyz.size() / 3;
  corr     .assign(nnodes, -1);
  corr_dist.assign(nnodes, 0.0f);

  PROGRESS<size_t> progress(nnodes, "Correspondance progress: ");

  kdtree vert_tree(100);
  vert_tree.build_vertex_tree(mesh2.xyz);

  int outside = 0;

  #ifdef OPENMP
  #pragma omp parallel for schedule(guided) reduction(+:outside)
  #endif
  for(size_t i = 0; i < corr.size(); i++) {
    vec3r ref(mesh1.xyz.data() + i*3);
    int nidx = -1;
    vec3r   closest_point;
    mt_real len2;

    bool inside = vert_tree.closest_vertex(ref, nidx, closest_point, len2);
    if(!inside) outside++;

    corr[i] = nidx;
    corr_dist[i] = len2;

    #ifdef OPENMP
    #pragma omp critical
    #endif
    {progress.next();}
  }

  progress.finish();

  if(outside) {
    fprintf(stderr, "%s warning: %d of %d queried vertices were not inside bbox of the used kdtree!\n",
            __func__, outside, (int)corr.size());
  }

}

void traverse_nodal_connectivity(const mt_vector<mt_int> & n2n_cnt,
                           const mt_vector<mt_int> & n2n_dsp,
                           const mt_vector<mt_int> & n2n_con,
                           mt_vector<bool> & reached)
{
  size_t nnodes = n2n_cnt.size();

  MT_USET<mt_int> cur_reached, next_reached;
  for(size_t i=0; i < nnodes; i++)
    if(reached[i]) cur_reached.insert(i);

  size_t nreached = cur_reached.size(), last_nreached = 0;

  while(last_nreached < nreached)
  {
    last_nreached = nreached;

    for(mt_int & nidx : cur_reached) {
      mt_int start = n2n_dsp[nidx], stop = start + n2n_cnt[nidx];
      for(mt_int j = start; j<stop; j++)
      {
        mt_int cidx = n2n_con[j];
        if(!reached[cidx]) {
          reached[cidx] = true;
          next_reached.insert(cidx);
          nreached++;
        }
      }
    }

    cur_reached.clear();
    cur_reached.insert(next_reached.begin(), next_reached.end());
    next_reached.clear();
  }
}

void traverse_nodal_connectivity(const mt_vector<mt_int> & n2n_cnt,
                           const mt_vector<mt_int> & n2n_dsp,
                           const mt_vector<mt_int> & n2n_con,
                           const MT_USET<mt_int> & blocked,
                           mt_vector<bool> & reached)
{
  size_t nnodes = n2n_cnt.size();

  MT_USET<mt_int> cur_reached, next_reached;
  for(size_t i=0; i < nnodes; i++)
    if(reached[i]) cur_reached.insert(i);

  size_t nreached = cur_reached.size(), last_nreached = 0;

  while(last_nreached < nreached)
  {
    last_nreached = nreached;

    for(mt_int & nidx : cur_reached)
    {
      if(blocked.count(nidx) == 0)
      {
        mt_int start = n2n_dsp[nidx], stop = start + n2n_cnt[nidx];
        for(mt_int j = start; j<stop; j++)
        {
          mt_int cidx = n2n_con[j];
          if(!reached[cidx]) {
            reached[cidx] = true;
            next_reached.insert(cidx);
            nreached++;
          }
        }
      }
    }

    cur_reached.clear();
    cur_reached.insert(next_reached.begin(), next_reached.end());
    next_reached.clear();
  }
}

void traverse_surfelem_connectivity(const mt_meshdata & surfmesh,
                                    const mt_mapping<mt_int> & ele2edge,
                                    const MT_USET<mt_int> & block_edges,
                                    mt_vector<bool> & reached)
{
  check_same_size(surfmesh.e2n_cnt, reached, __func__);
  size_t nelem = surfmesh.e2n_cnt.size();
  MT_USET<mt_int> cur_reached, next_reached;

  for(size_t i=0; i < nelem; i++)
    if(reached[i]) cur_reached.insert(i);

  size_t nreached = cur_reached.size(), last_nreached = 0;

  while(last_nreached < nreached)
  {
    last_nreached = nreached;

    // sweep over connectivity and try to reach new nodes
    for(mt_int eidx : cur_reached)
    {
      // we look at edges of reached element
      mt_int start = ele2edge.fwd_dsp[eidx], stop = start + ele2edge.fwd_cnt[eidx];
      for(mt_int j = start; j < stop; j++) {
        mt_int edge_idx = ele2edge.fwd_con[j];

        // check if current edge is sharp
        if(block_edges.count(edge_idx) == 0) {
          // if edge is not sharp, we iterate over connected elements and mark them as reached
          mt_int estart = ele2edge.bwd_dsp[edge_idx],
                 estop  = estart + ele2edge.bwd_cnt[edge_idx];

          for(mt_int k = estart; k < estop; k++) {
            mt_int reached_eidx = ele2edge.bwd_con[k];
            if(mt_int(eidx) != reached_eidx && reached[reached_eidx] == false) {
              reached[reached_eidx] = true;
              next_reached.insert(reached_eidx);
              nreached++;
            }
          }
        }
      }
    }

    cur_reached.clear();
    cur_reached.insert(next_reached.begin(), next_reached.end());
    next_reached.clear();
  }
}

void traverse_surfelem_connectivity(const mt_meshdata & surfmesh,
                                    const mt_mapping<mt_int> & ele2edge,
                                    const mt_vector<mt_real> & elem_normals,
                                    const vec3r ref_nrml,
                                    mt_real thr,
                                    mt_vector<bool> & reached)
{
  check_same_size(surfmesh.e2n_cnt, reached, __func__);
  size_t nelem = surfmesh.e2n_cnt.size();

  // thr is assumed in degrees, we convert it back to a length
  thr = cos(thr * (MT_PI/ 180.0));

  MT_USET<mt_int> cur_reached, next_reached;

  for(size_t i=0; i < nelem; i++)
    if(reached[i]) cur_reached.insert(i);

  size_t nreached = cur_reached.size(), last_nreached = 0;

  while(last_nreached < nreached)
  {
    last_nreached = nreached;

    // sweep over connectivity and try to reach new nodes
    for(mt_int eidx : cur_reached)
    {
      // we look at edges of reached element
      mt_int start = ele2edge.fwd_dsp[eidx], stop = start + ele2edge.fwd_cnt[eidx];
      for(mt_int j = start; j < stop; j++) {
        mt_int edge_idx = ele2edge.fwd_con[j];

        // we iterate over connected elements and check angle
        mt_int estart = ele2edge.bwd_dsp[edge_idx],
        estop  = estart + ele2edge.bwd_cnt[edge_idx];

        for(mt_int k = estart; k < estop; k++) {
          mt_int reached_eidx = ele2edge.bwd_con[k];
          if(mt_int(eidx) != reached_eidx && reached[reached_eidx] == false)
          {
            vec3r n(elem_normals.data() + reached_eidx*3);
            if(n.scaProd(ref_nrml) > thr) {
              reached[reached_eidx] = true;
              next_reached.insert(reached_eidx);
              nreached++;
            }
          }
        }
      }
    }

    cur_reached.clear();
    cur_reached.insert(next_reached.begin(), next_reached.end());
    next_reached.clear();
  }
}

void traverse_nodal_connectivity(const mt_vector<mt_int> & n2n_cnt,
                           const mt_vector<mt_int> & n2n_dsp,
                           const mt_vector<mt_int> & n2n_con,
                           const mt_vector<mt_real> & xyz,
                           mt_point<mt_real> & ref,
                           const mt_real rad,
                           mt_vector<bool> & reached)
{
  size_t nnodes = n2n_cnt.size();
  const mt_real squared_rad = rad*rad;

  MT_USET<mt_int> cur_reached, next_reached;
  for(size_t i=0; i < nnodes; i++)
    if(reached[i]) cur_reached.insert(i);

  size_t nreached = cur_reached.size(), last_nreached = 0;

  while(last_nreached < nreached)
  {
    last_nreached = nreached;

    for(mt_int & nidx : cur_reached)
    {
      mt_int start = n2n_dsp[nidx], stop = start + n2n_cnt[nidx];
      for(mt_int j = start; j<stop; j++)
      {
        mt_int cidx = n2n_con[j];
        mt_point<mt_real> e = ref - mt_point<mt_real>(xyz.data() + 3*cidx);
        if( (!reached[cidx]) && (e.length2() < squared_rad) ) {
          reached[cidx] = true;
          next_reached.insert(cidx);
          nreached++;
        }
      }
    }

    cur_reached.clear();
    cur_reached.insert(next_reached.begin(), next_reached.end());
    next_reached.clear();
  }
}

void nodal_connectivity_decomposition(const mt_meshdata & mesh,
                                      mt_vector<mt_meshdata> & parts,
                                      mt_vector<mt_vector<mt_int> > & parts_eidx)
{
  check_condition(mesh.n2n_dsp.size() > 0 && mesh.e2n_dsp.size() > 0,
                  "mesh n2n connectivity set up", __func__);

  size_t nelem = mesh.e2n_cnt.size();
  size_t nnode = mesh.n2e_cnt.size();
  size_t initidx, iter = 0;
  mt_vector<bool> reached;
  mt_vector<mt_int> nodal_part(nnode, -1);
  bool finished = false;

  // the vertices in the mesh
  mt_vector<mt_int> vtx(mesh.e2n_con);
  binary_sort(vtx); unique_resize(vtx);
  initidx = vtx[0];

  while(!finished) {
    reached.assign(nnode, false);
    mt_meshdata       & meshpart  = parts     .push_back(mt_meshdata());
    mt_vector<mt_int> & extr_eidx = parts_eidx.push_back(mt_vector<mt_int>());

    extr_eidx.reserve(nelem / 2);
    reached[initidx] = true;
    traverse_nodal_connectivity(mesh.n2n_cnt, mesh.n2n_dsp, mesh.n2n_con, reached);

    for(size_t eidx = 0; eidx < mesh.e2n_cnt.size(); eidx++) {
      mt_int start = mesh.e2n_dsp[eidx], stop = start + mesh.e2n_cnt[eidx];
      mt_int numreached = 0;

      for(mt_int i = start; i < stop; i++) {
        mt_int nidx = mesh.e2n_con[i];
        if(reached[nidx]) numreached++;
      }

      if(numreached == mesh.e2n_cnt[eidx])
        extr_eidx.push_back(eidx);
    }

    extract_mesh(extr_eidx, mesh, meshpart);

    for(auto v : vtx)
      if(reached[v]) nodal_part[v] = iter;

    // we finish if all nodes in the mesh have been partitioned
    finished = true;
    for(auto v : vtx) {
      if(nodal_part[v] == -1) {
        initidx = v;
        finished = false;
        break;
      }
    }

    iter++;
  }
}

void surf_connectivity_decomposition(const mt_meshdata & surfmesh,
                                     const mt_mapping<mt_int> & ele2edge,
                                     const MT_USET<mt_int> & block_edges,
                                     mt_vector<mt_meshdata> & parts,
                                     mt_vector<mt_vector<mt_int> > & parts_eidx)
{
  size_t nelem = surfmesh.e2n_cnt.size();
  size_t initidx = 0, iter = 0;
  mt_vector<bool> reached;
  mt_vector<mt_int> part(nelem, -1);
  bool finished = false;

  while(!finished) {
    reached.assign(nelem, false);
    mt_meshdata       & meshpart  = parts     .push_back(mt_meshdata());
    mt_vector<mt_int> & extr_eidx = parts_eidx.push_back(mt_vector<mt_int>());

    extr_eidx.reserve(nelem / 2);
    reached[initidx] = true;
    traverse_surfelem_connectivity(surfmesh, ele2edge, block_edges, reached);

    for(size_t i=0; i<nelem; i++) {
      if(reached[i]) {
        extr_eidx.push_back(i);
        part[i] = iter;
      }
    }

    extract_mesh(extr_eidx, surfmesh, meshpart);

    // we finish if all nodes in the mesh have been partitioned
    finished = true;

    for(size_t i=0; i<nelem; i++) {
      if(part[i] == -1) {
        initidx = i;
        finished = false;
        break;
      }
    }

    iter++;
  }
}

bool tet_fix_face_orientation(mt_meshdata & mesh, mt_int eidx, short face_idx)
{
  mt_int* nod = mesh.e2n_con.data() + mesh.e2n_dsp[eidx];
  mt_int n1 = nod[0], n2 = nod[1], n3 = nod[2], n4 = nod[3];

  vec3r ele_ctr = barycenter(4, nod, mesh.xyz.data());
  vec3r c, n;

  switch(face_idx) {
    default:
      fprintf(stderr, "%s error: Tet face index not in range [1,4]!\n", __func__);
      exit(EXIT_FAILURE);

    case 1:
      // (2,3,1)
      n = triangle_normal(n2, n3, n1, mesh.xyz);
      c = triangle_centerpoint(n2, n3, n1, mesh.xyz);
      break;

    case 2:
      // (1,4,2)
      n = triangle_normal     (n1, n4, n2, mesh.xyz);
      c = triangle_centerpoint(n1, n4, n2, mesh.xyz);
      break;

    case 3:
      // (2,4,3)
      n = triangle_normal     (n2, n4, n3, mesh.xyz);
      c = triangle_centerpoint(n2, n4, n3, mesh.xyz);
      break;

    case 4:
      // (1,3,4)
      n = triangle_normal     (n1, n3, n4, mesh.xyz);
      c = triangle_centerpoint(n1, n3, n4, mesh.xyz);
      break;
  }

  vec3r v = unit_vector(ele_ctr - c);
  mt_real sca = n.scaProd(v);

  if(sca > 0) {
    switch(face_idx) {
      default:
        fprintf(stderr, "%s error: Tet face index not in range [1,4]!\n", __func__);
        exit(EXIT_FAILURE);

      case 1:
        // (2,3,1)
        nod[1] = n1,
        nod[2] = n3,
        nod[0] = n2;
        break;

      case 2:
        // (1,4,2)
        nod[0] = n2,
        nod[3] = n4,
        nod[1] = n1;
        break;

      case 3:
        // (2,4,3)
        nod[1] = n3,
        nod[3] = n4,
        nod[2] = n2;
        break;

      case 4:
        // (1,3,4)
        nod[0] = n4,
        nod[2] = n3,
        nod[3] = n1;
        break;
    }
    return true;
  }

  return false;
}

void correct_insideOut(mt_meshdata & mesh)
{
  size_t nelem = mesh.e2n_cnt.size();

  if(mesh.e2n_dsp.size() != nelem) {
    mesh.e2n_dsp.resize(mesh.e2n_cnt.size());
    bucket_sort_offset(mesh.e2n_cnt, mesh.e2n_dsp);
  }

  size_t corr = 0, corr_old = 0, corr_pass = 0, max = 5;

  do {
    corr_old = corr;
    corr_pass++;
    for(size_t i=0; i<nelem; i++) {
      if(mesh.etype[i] == Tetra) {
        #if 1
        const mt_int edsp = mesh.e2n_dsp[i];
        const mt_int v0 = mesh.e2n_con[edsp+0];
        const mt_int v1 = mesh.e2n_con[edsp+1];
        const mt_int v2 = mesh.e2n_con[edsp+2];
        const mt_int v3 = mesh.e2n_con[edsp+3];

        mt_point<mt_real> p0(mesh.xyz.data() + v0*3);
        mt_point<mt_real> p1(mesh.xyz.data() + v1*3);
        mt_point<mt_real> p2(mesh.xyz.data() + v2*3);
        mt_point<mt_real> p3(mesh.xyz.data() + v3*3);

        const mt_real vol = signed_tet_volume(p0,p1,p2,p3);

        if(vol < 0) {
          mesh.e2n_con[edsp+0] = v2;
          mesh.e2n_con[edsp+1] = v1;
          mesh.e2n_con[edsp+2] = v0;
          mesh.e2n_con[edsp+3] = v3;
          corr++;
        }
        #else
        // in this branch we fix the faces individually. I have not come across a case
        // where it was better than the one-time check with signed_tet_volume(). -Aurel
        bool fix = false;
        fix |= tet_fix_face_orientation(mesh, i, 1);
        fix |= tet_fix_face_orientation(mesh, i, 2);
        fix |= tet_fix_face_orientation(mesh, i, 3);
        fix |= tet_fix_face_orientation(mesh, i, 4);

        if(fix) corr++;
        #endif
      }
    }
    std::cout << "Pass " << corr_pass << ": Corrected " << corr - corr_old <<
                 " inside-out elements." << std::endl;
  } while( (corr > corr_old) && (corr_pass < max) );
}

void correct_insideOut_tri(mt_meshdata & mesh, const mt_vector<mt_real> & surf_nrml)
{
  size_t nelem = mesh.e2n_cnt.size();
  size_t corr = 0, corr_old = 0, corr_pass = 0, max = 5;

  do {
    corr_old = corr;
    corr_pass++;
    for(size_t i=0; i<nelem; i++) {
      assert(mesh.etype[i] == Tri);

      const mt_int v0 = mesh.e2n_con[i*3+0];
      const mt_int v1 = mesh.e2n_con[i*3+1];
      const mt_int v2 = mesh.e2n_con[i*3+2];

      mt_point<mt_real> p0(mesh.xyz.data() + v0*3);
      mt_point<mt_real> p1(mesh.xyz.data() + v1*3);
      mt_point<mt_real> p2(mesh.xyz.data() + v2*3);
      mt_point<mt_real> n0(surf_nrml.data() + v0*3);
      mt_point<mt_real> n1(surf_nrml.data() + v1*3);
      mt_point<mt_real> n2(surf_nrml.data() + v2*3);

      mt_point<mt_real> e01 = p1 - p0;
      mt_point<mt_real> e02 = p2 - p0;

      mt_point<mt_real> n  = e01.crossProd(e02) * mt_real(-1);
      mt_point<mt_real> na = (n0 + n1 + n2) / mt_real(3.0);
      mt_real sca = n.scaProd(na);

      if(sca < 0) {
        mesh.e2n_con[i*3+0] = v2;
        mesh.e2n_con[i*3+1] = v1;
        mesh.e2n_con[i*3+2] = v0;
        corr++;
      }
    }
    std::cout << "Pass " << corr_pass << ": Corrected " << corr - corr_old <<
                 " inside-out elements." << std::endl;
  } while( (corr > corr_old) && (corr_pass < max) );
}

void correct_duplicate_elements(mt_meshdata & mesh)
{
  MT_USET<quadruple<mt_int> > tet_elems;
  MT_USET<triple<mt_int> > tri_elems;
  mt_vector<bool> keep(mesh.e2n_cnt.size(), false);
  mt_vector<mt_int> elems, qbuff(4);
  size_t dupl = 0;

  for(size_t eidx = 0; eidx < mesh.e2n_cnt.size(); eidx++)
  {
    mt_int estart = mesh.e2n_dsp[eidx];
    switch(mesh.etype[eidx]) {
      case Tetra:
      {
        qbuff.assign(mesh.e2n_con.data() + estart, mesh.e2n_con.data() + estart+4);
        binary_sort(qbuff);
        quadruple<mt_int> q = {qbuff[0], qbuff[1], qbuff[2], qbuff[3]};
        auto it = tet_elems.find(q);
        if(it == tet_elems.end()) {
          tet_elems.insert(q);
          keep[eidx] = true;
        }
        else dupl++;
        break;
      }
      case Tri:
      {
        triple<mt_int> q;
        sortTriple(mesh.e2n_con[estart+0], mesh.e2n_con[estart+1], mesh.e2n_con[estart+2],
                   q.v1, q.v2, q.v3);
        auto it = tri_elems.find(q);
        if(it == tri_elems.end()) {
          tri_elems.insert(q);
          keep[eidx] = true;
        }
        else dupl++;
        break;
      }
      default: break;
    }
  }

  if(dupl) {
    printf("Removing %lu element duplicates..\n", dupl);
    restrict_meshdata(keep, mesh, elems);
  }
}

void correct_duplicate_vertices(mt_meshdata & mesh, bool verbose)
{
  int scale_factor;
  MT_MAP<triple<mt_int>,mt_int> coords;
  int did_correct = 0;
  triple<mt_int> t;

  if(mesh.e2n_dsp.size() == 0)
    bucket_sort_offset(mesh.e2n_cnt, mesh.e2n_dsp);

  generate_coord_map(mesh, coords, scale_factor);

  for(size_t eidx = 0, ridx = 0; eidx < mesh.e2n_cnt.size(); eidx++) {
    for(mt_int i = 0; i < mesh.e2n_cnt[eidx]; i++, ridx++) {
      mt_int v = mesh.e2n_con[ridx];
      t.v1 = mt_int(mesh.xyz[v*3+0]*scale_factor);
      t.v2 = mt_int(mesh.xyz[v*3+1]*scale_factor);
      t.v3 = mt_int(mesh.xyz[v*3+2]*scale_factor);

      auto it = coords.find(t);
      if(it != coords.end()) {
        mt_int v2 = it->second;
        if(v != v2) {
          if(verbose) {
            printf("Duplicate vertex: %ld (%g, %g, %g) == %ld (%g, %g, %g)\n",
                   v,  mesh.xyz[v*3+0],  mesh.xyz[v*3+1],  mesh.xyz[v*3+2],
                   v2, mesh.xyz[v2*3+0], mesh.xyz[v2*3+1], mesh.xyz[v2*3+2]);
          }
          mesh.e2n_con[ridx] = v2;
          did_correct++;
        }
      }
      else {
        fprintf(stderr, "%s error: coord map inconsistent with mesh. Aborting!\n", __func__);
        exit(1);
      }
    }
  }

  // if we remaped nodes, we do a re-indexing
  if(did_correct > 0) {
    printf("Correcting %d duplicate nodes..\n", did_correct);
    reindex_nodes(mesh, true);
  }

  // if we remaped nodes, we also need to check if we have now ill-defined elements
  if(did_correct > 0) {
    MT_USET<mt_int> elemnodes;
    mt_vector<bool> keep(mesh.e2n_cnt.size(), true);
    int badelem = 0;

    for(size_t eidx = 0, ridx = 0; eidx < mesh.e2n_cnt.size(); eidx++) {
      elemnodes.clear();
      for(mt_int i = 0; i < mesh.e2n_cnt[eidx]; i++, ridx++)
        elemnodes.insert(mesh.e2n_con[ridx]);

      if(elemnodes.size() != size_t(mesh.e2n_cnt[eidx])) {
        if(verbose) {
          printf("Ill defined element %ld: ", (long int) eidx);
          for(mt_int i = mesh.e2n_dsp[eidx], j = 0; j<mesh.e2n_cnt[eidx]; i++, j++)
            printf("%ld ", (long int) mesh.e2n_con[i]);
          printf("\n");
        }
        badelem++;
        keep[eidx] = false;
      }
    }

    if(badelem) {
      mt_vector<mt_int> elems;
      printf("Removing %d ill-defined elements\n", badelem);
      restrict_meshdata(keep, mesh, elems);
    }
  }
}

bool has_bad_vol(const elem_t type,
                 const mt_int* con,
                 const mt_real* xyz,
                 const mt_real bad_thr)
{
  bool ret = false;

  switch(type) {
    case Tri:
    {
      const mt_int v0 = con[0], v1 = con[1], v2 = con[2];

      vec3r p0(xyz + v0*3);
      vec3r p1(xyz + v1*3);
      vec3r p2(xyz + v2*3);

      vec3r base_ctr = (p0 + p1) / 2.0;
      vec3r e01 = unit_vector(p1 - p0);
      vec3r e02 = unit_vector(p2 - p0);

      vec3r n  = unit_vector(e02 - (e01 * e01.scaProd(e02)));
      vec3r c2 = unit_vector(p2 - base_ctr);
      mt_real sca = n.scaProd(c2);

      if(sca < bad_thr)
        ret = true;

      break;
    }

    case Tetra:
    {
      const mt_int v0 = con[0], v1 = con[1], v2 = con[2], v3 = con[3];

      vec3r p0(xyz + v0*3);
      vec3r p1(xyz + v1*3);
      vec3r p2(xyz + v2*3);
      vec3r p3(xyz + v3*3);

      vec3r base_ctr = (p0 + p1 + p2) / 3.0;
      vec3r e01 = unit_vector(p1 - p0);
      vec3r e02 = unit_vector(p2 - p0);
      vec3r c3  = unit_vector(p3 - base_ctr);

      // we dont really check the volume, but rather if the "top point vector"
      // e03 has a significant component on the normal vector
      vec3r n = e01.crossProd(e02);
      mt_real sca = n.scaProd(c3);

      if(sca < bad_thr)
        ret = true;

      break;
    }

    default: break;
  }

  return ret;
}

void tet_get_angles(const mt_int* con,
                    const mt_vector<mt_real> & xyz,
                    std::set<mt_real> & angles)
{
  const mt_int v0 = con[0];
  const mt_int v1 = con[1];
  const mt_int v2 = con[2];
  const mt_int v3 = con[3];

  mt_point<mt_real> p0(xyz.data()+v0*3);
  mt_point<mt_real> p1(xyz.data()+v1*3);
  mt_point<mt_real> p2(xyz.data()+v2*3);
  mt_point<mt_real> p3(xyz.data()+v3*3);

  mt_point<mt_real> e01 = p1 - p0; e01.normalize();
  mt_point<mt_real> e02 = p2 - p0; e02.normalize();
  // mt_point<mt_real> e12 = p2 - p1; e12.normalize();

  mt_point<mt_real> e03 = p3 - p0; e03.normalize();
  mt_point<mt_real> e13 = p3 - p1; e13.normalize();
  mt_point<mt_real> e23 = p3 - p2; e23.normalize();

  mt_point<mt_real> n = e01.crossProd(e02); n.normalize();

  // mt_real a = acos(e01.scaProd(e02)) * (180.0 / MT_PI),
  //   b = 180.0 - acos(e01.scaProd(e12)) * (180.0 / MT_PI),
  //   c = 180.0 - a - b;

  // angles.insert( a );
  // angles.insert( b );
  // angles.insert( c );

  mt_real a, b, c;
  a = acos(n.scaProd(e03)) * (180.0 / MT_PI);
  a = a < 90.0 ? a : 180.0 - a;
  b = acos(n.scaProd(e13)) * (180.0 / MT_PI);
  b = b < 90.0 ? b : 180.0 - b;
  c = acos(n.scaProd(e23)) * (180.0 / MT_PI);
  c = c < 90.0 ? c : 180.0 - c;

  angles.insert( a );
  angles.insert( b );
  angles.insert( c );
}

void write_graph(const mt_vector<mt_int> & cnt,
                 const mt_vector<mt_int> & con,
                 const char* filename)
{
  size_t ncon = con.size();
  mt_vector<int>    lcnt(cnt.size());
  mt_vector<int>    lcol(ncon);
  mt_vector<double> lele(ncon);

  int row_size = htobe(int(lcnt.size()));
  int col_max = *std::max_element(lcol.begin(), lcol.end());
  int col_size = htobe(int(col_max + 2));

  for(size_t i=0; i<cnt.size(); i++)
    lcnt[i] = htobe(int(cnt[i]));

  for(size_t i=0; i<ncon; i++) {
    lcol[i] = htobe(int(con[i] + 1));
    lele[i] = htobe(double(1));
  }

  int petscMatIndex = htobe(int(1211216));
  int nele = htobe(int(ncon));

  FILE* fd = fopen(filename, MT_FOPEN_WRITE);
  if(!fd) treat_file_open_error(filename, errno);

  fwrite(&petscMatIndex, sizeof(int)   , 1         , fd);
  fwrite(&row_size     , sizeof(int)   , 1         , fd); // m
  fwrite(&col_size     , sizeof(int)   , 1         , fd); // n
  fwrite(&nele         , sizeof(int)   , 1         , fd); // nz
  fwrite(lcnt.data()   , sizeof(int)   , lcnt.size(), fd); // nnz .. nonzeros per row
  fwrite(lcol.data()   , sizeof(int)   , lcol.size(), fd); // j   .. column indices
  fwrite(lele.data()   , sizeof(double), lele.size(), fd); // s   .. element data

  fclose(fd);
}


void mt_shmem_partitioner::compute_partitions(const mt_meshdata & mesh,
      const int npart)
{
  assert(mesh.n2n_cnt.size() > 0);
  size_t nnodes = mesh.n2n_cnt.size();
  size_t psize = (nnodes) / npart;

  part.resize(npart);
  itf.assign(nnodes, false);

  // loop over partitions
  for(int part_idx=0; part_idx < npart; part_idx++)
  {
    // define partition start and end node indices
    mt_int part_start = part_idx*psize;
    mt_int part_end   = part_idx < npart - 1 ? part_start + psize : nnodes;
    // loop over nodes of current partition
    for(mt_int nidx = part_start; nidx < part_end; nidx++) {
      // add nodes to partition set
      part[part_idx].insert(nidx);
      mt_int row_start = mesh.n2n_dsp[nidx], row_end = row_start + mesh.n2n_cnt[nidx];
      for(mt_int j = row_start; j<row_end; j++)
      {
        mt_int cidx = mesh.n2n_con[j];
        // mark connected nodes outside of the partition interval as interface nodes
        if(cidx < part_start || cidx > part_end)
          itf[cidx] = true;
      }
    }
  }
}

void enclosing_element(const mt_meshdata & mesh,
                       const MT_USET<mt_int> & nod,
                       const mt_point<mt_real> & pt,
                       const mt_real edge_len,
                       mt_int & elem,
                       mt_real & a,
                       mt_real & b,
                       mt_real & c)
{
  mt_real detJ;
  elem = -1, a = -1., b = -1., c = -1.;

  for(const mt_int & nidx : nod) {
    mt_int eidx_start = mesh.n2e_dsp[nidx],
           eidx_stop  = eidx_start + mesh.n2e_cnt[nidx];

    for(mt_int i=eidx_start; i<eidx_stop; i++)
    {
      const mt_int eidx = mesh.n2e_con[i];
      const mt_int* con = mesh.e2n_con.data() + mesh.e2n_dsp[eidx];
      const mt_real* xyz = mesh.xyz.data();
      const mt_real eps = edge_len * 1e-6;

      mt_point<mt_real> p0(xyz + con[0]*3);
      mt_point<mt_real> t = pt - p0;

      bool found = false;
      switch(mesh.etype[eidx]) {
        case Tri:
        {
          // for triangles we use the barycentric coordinates of the projected
          // point.
          mt_point<mt_real> e1 = mt_point<mt_real>(xyz + con[1]*3) - p0;
          mt_point<mt_real> e2 = mt_point<mt_real>(xyz + con[2]*3) - p0;
          mt_point<mt_real> n = e1.crossProd(e2);
          mt_point<mt_real> v = e1.crossProd(t);
          mt_point<mt_real> w = t.crossProd(e2);
          mt_real nsq = n.scaProd(n);
          c = v.scaProd(n) / nsq;
          b = w.scaProd(n) / nsq;
          a = 1.0 - c - b;

          // if we find the enclosing element we set eidx and return
          if(a > -eps && b > -eps && c > -eps)
            found = true;
          break;
        }
        case Tetra:
        {
          // for tets we compute the mapping x -> xi where x is a point in eulerian
          // coordinates and xi in local element coordinates
          //
          // we have x = p0 + [e1 e2 e3] xi = p0 + J xi,
          // thus xi = J^-1 (x - p0) = J^-1 t
          mt_point<mt_real> e1 = mt_point<mt_real>(xyz + con[1]*3) - p0;
          mt_point<mt_real> e2 = mt_point<mt_real>(xyz + con[2]*3) - p0;
          mt_point<mt_real> e3 = mt_point<mt_real>(xyz + con[3]*3) - p0;
          dmat<mt_real> J(3,3);
          J[0][0] = e1.x; J[0][1] = e2.x; J[0][2] = e3.x;
          J[1][0] = e1.y; J[1][1] = e2.y; J[1][2] = e3.y;
          J[2][0] = e1.z; J[2][1] = e2.z; J[2][2] = e3.z;
          // the inverse jacobian is the mapping (x,y,z) -> (xi, eta, gamma)
          // in the 3D-triangle case, we need a pseudoinverse
          invert_3x3(J, detJ);
          // store the left pseudoinverse in J
          a = J[0][0] * t.x + J[0][1] * t.y + J[0][2] * t.z;
          b = J[1][0] * t.x + J[1][1] * t.y + J[1][2] * t.z;
          c = J[2][0] * t.x + J[2][1] * t.y + J[2][2] * t.z;
          mt_real gamma = 1 - a - b - c;
          // if we find the enclosing element we set eidx and return
          if(a > -eps && b > -eps && c > -eps && gamma > -eps)
            found = true;
          break;
        }
        case Prism:
        {
          //Check if point lies insided of the prism or not by checking whether it is behind every plane
          // see https://stackoverflow.com/questions/8877872/determining-if-a-point-is-inside-a-polyhedron

          //get the remaining 5 points of the prism
          mt_point<mt_real> p1 = mt_point<mt_real>(xyz + con[1]*3);
          mt_point<mt_real> p2 = mt_point<mt_real>(xyz + con[2]*3);
          mt_point<mt_real> p3 = mt_point<mt_real>(xyz + con[3]*3);
          mt_point<mt_real> p4 = mt_point<mt_real>(xyz + con[4]*3);
          mt_point<mt_real> p5 = mt_point<mt_real>(xyz + con[5]*3);

          //faces of the prisms have to be oriented counterclockwise.
          //faces are given by {0,2,1}, {3,5,4}, {1,2,4,5}, {2,0,3,4}, {0,1,5,4}
          mt_vector<vec3r> normals(5);
          normals[0] = ((p2-p0).crossProd(p1-p0)); normals[0].normalize();
          normals[1] = ((p5-p3).crossProd(p4-p3)); normals[1].normalize();
          normals[2] = ((p2-p1).crossProd(p4-p1)); normals[2].normalize();
          normals[3] = ((p0-p2).crossProd(p3-p2)); normals[3].normalize();
          normals[4] = ((p1-p0).crossProd(p5-p0)); normals[4].normalize();
          //Barycenters of the faces
          mt_vector<vec3r> bary(5);
          bary[0] = (p0+p2+p1) / 3.0;
          bary[1] = (p3+p5+p4) / 3.0;
          bary[2] = (p1+p2+p4+p5) / 4.0;
          bary[3] = (p2+p1+p3+p4) / 4.0;
          bary[4] = (p0+p1+p5+p3) / 4.0;
          //counter counting occurences of violation of constraint
          size_t neg_cnt = 0;
          for(size_t k=0; k < 5; k++) {
            mt_point<mt_real> p2f = bary[k] - pt;
            mt_real d = p2f.scaProd(normals[k]);
            d /= p2f.length();
            if( d < -1e-15)
              neg_cnt++;
          }
          found = neg_cnt == 0;
          break;
        }
        case Hexa:
        {
          //get the remaining 7 points of the prism
          mt_point<mt_real> p1 = mt_point<mt_real>(xyz + con[1]*3);
          mt_point<mt_real> p2 = mt_point<mt_real>(xyz + con[2]*3);
          mt_point<mt_real> p3 = mt_point<mt_real>(xyz + con[3]*3);
          mt_point<mt_real> p4 = mt_point<mt_real>(xyz + con[4]*3);
          mt_point<mt_real> p5 = mt_point<mt_real>(xyz + con[5]*3);
          mt_point<mt_real> p6 = mt_point<mt_real>(xyz + con[6]*3);
          mt_point<mt_real> p7 = mt_point<mt_real>(xyz + con[7]*3);

          //faces of the hex have to be oriented counterclockwise.
          //faces are given by {0,3,2,1}, {1,2,6,7}, {4,7,6,5}, {0,4,5,3}, {2,3,5,6}, {0,1,7,4}
          mt_vector<vec3r> normals(6);
          normals[0] = ((p3-p0).crossProd(p2-p0)); normals[0].normalize();
          normals[1] = ((p2-p1).crossProd(p6-p1)); normals[1].normalize();
          normals[2] = ((p7-p4).crossProd(p6-p4)); normals[2].normalize();
          normals[3] = ((p4-p0).crossProd(p5-p0)); normals[3].normalize();
          normals[4] = ((p3-p2).crossProd(p5-p2)); normals[4].normalize();
          normals[5] = ((p1-p0).crossProd(p7-p0)); normals[5].normalize();

          //Barycenters of the faces
          mt_vector<vec3r> bary(6);
          bary[0] = (p0+p3+p2+p1) / 4.0;
          bary[1] = (p1+p2+p6+p7) / 4.0;
          bary[2] = (p4+p7+p6+p5) / 4.0;
          bary[3] = (p0+p4+p5+p3) / 4.0;
          bary[4] = (p2+p3+p5+p6) / 4.0;
          bary[5] = (p0+p1+p7+p4) / 4.0;
          //counter counting occurences of violation of constraint
          size_t neg_cnt = 0;
          for(size_t k=0; k < 6; k++) {
            mt_point<mt_real> p2f = bary[k] - pt;
            mt_real d = p2f.scaProd(normals[k]);
            d /= p2f.length();
            if( d < -1e-15)
              neg_cnt++;
          }
          found = neg_cnt == 0;
          break;
        }
        default:
          fprintf(stderr, "enclose-check error: Unsupported elem type.\n");
          exit(EXIT_FAILURE);
      }
      if(found) {
        elem = eidx;
        return;
      }
    }
  }
}

void expand_nodeset(const mt_meshdata & mesh,
                    MT_USET<mt_int> & nod)
{
  MT_USET<mt_int> oldnod = nod;
  nod.clear();

  for(auto nidx : oldnod) {
    mt_int start = mesh.n2n_dsp[nidx],
           stop  = start + mesh.n2n_cnt[nidx];
    for(mt_int i=start; i<stop; i++) nod.insert(mesh.n2n_con[i]);
  }
}

void expand_in_radius(const mt_meshdata & mesh,
                      const mt_point<mt_real> ref,
                      const mt_real squared_rad,
                      std::set<mt_int> & nod)
{
  size_t old_nodsize;
  do {
    old_nodsize = nod.size();
    for(auto nidx : nod) {
      mt_int start = mesh.n2n_dsp[nidx],
             stop  = start + mesh.n2n_cnt[nidx];
      for(mt_int i=start; i<stop; i++) {
        const mt_int cidx = mesh.n2n_con[i];
        mt_point<mt_real> e = mt_point<mt_real>(mesh.xyz.data() + cidx*3) - ref;
        if(e.length2() < squared_rad) nod.insert(cidx);
      }
    }
  } while(old_nodsize < nod.size());
}

void element_edges_stats(const elem_t type, const mt_int* con, const mt_real* xyz,
                           mt_real & min_edge, mt_real & max_edge, mt_real & avg_edge)
{
  mt_real num_edges = 0;
  min_edge = 1e200, max_edge = 0, avg_edge = 0;

  switch(type) {
    case Prism:
    {
      const mt_int v1 = con[0], v2 = con[1], v3 = con[2], v4 = con[3], v5 = con[4], v6 = con[5];
      mt_real edge;

      const mt_real v1x = xyz[v1*3+0], v1y = xyz[v1*3+1], v1z = xyz[v1*3+2];
      const mt_real v2x = xyz[v2*3+0], v2y = xyz[v2*3+1], v2z = xyz[v2*3+2];
      const mt_real v3x = xyz[v3*3+0], v3y = xyz[v3*3+1], v3z = xyz[v3*3+2];
      const mt_real v4x = xyz[v4*3+0], v4y = xyz[v4*3+1], v4z = xyz[v4*3+2];
      const mt_real v5x = xyz[v5*3+0], v5y = xyz[v5*3+1], v5z = xyz[v5*3+2];
      const mt_real v6x = xyz[v6*3+0], v6y = xyz[v6*3+1], v6z = xyz[v6*3+2];

      // (3,1)
      edge = std::sqrt(POW2(v3x-v1x) + POW2(v3y-v1y) + POW2(v3z-v1z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      // (3,2)
      edge = std::sqrt(POW2(v3x-v2x) + POW2(v3y-v2y) + POW2(v3z-v2z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      // (2,1)
      edge = std::sqrt(POW2(v2x-v1x) + POW2(v2y-v1y) + POW2(v2z-v1z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      // (5,4)
      edge = std::sqrt(POW2(v5x-v4x) + POW2(v5y-v4y) + POW2(v5z-v4z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      // (6,5)
      edge = std::sqrt(POW2(v5x-v6x) + POW2(v5y-v6y) + POW2(v5z-v6z));
      edge = std::sqrt(POW2(v5x-v6x) + POW2(v5y-v6y) + POW2(v5z-v6z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      // (6,4)
      edge = std::sqrt(POW2(v4x-v6x) + POW2(v4y-v6y) + POW2(v4z-v6z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      // (1,4)
      edge = std::sqrt(POW2(v1x-v4x) + POW2(v1y-v4y) + POW2(v1z-v4z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      // (3,5)
      edge = std::sqrt(POW2(v3x-v5x) + POW2(v3y-v5y) + POW2(v3z-v5z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      // (2,6)
      edge = std::sqrt(POW2(v2x-v6x) + POW2(v2y-v6y) + POW2(v2z-v6z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      num_edges = 9.0;
      break;
    }
    case Hexa:
    {
      const mt_int v1 = con[0], v2 = con[1], v3 = con[2], v4 = con[3], v5 = con[4], v6 = con[5], v7 = con[6], v8 = con[7];
      mt_real edge;

      const mt_real v1x = xyz[v1*3+0], v1y = xyz[v1*3+1], v1z = xyz[v1*3+2];
      const mt_real v2x = xyz[v2*3+0], v2y = xyz[v2*3+1], v2z = xyz[v2*3+2];
      const mt_real v3x = xyz[v3*3+0], v3y = xyz[v3*3+1], v3z = xyz[v3*3+2];
      const mt_real v4x = xyz[v4*3+0], v4y = xyz[v4*3+1], v4z = xyz[v4*3+2];
      const mt_real v5x = xyz[v5*3+0], v5y = xyz[v5*3+1], v5z = xyz[v5*3+2];
      const mt_real v6x = xyz[v6*3+0], v6y = xyz[v6*3+1], v6z = xyz[v6*3+2];
      const mt_real v7x = xyz[v7*3+0], v7y = xyz[v7*3+1], v7z = xyz[v7*3+2];
      const mt_real v8x = xyz[v8*3+0], v8y = xyz[v8*3+1], v8z = xyz[v8*3+2];

      // (1,2)
      edge = std::sqrt(POW2(v2x-v1x) + POW2(v2y-v1y) + POW2(v2z-v1z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      // (3,2)
      edge = std::sqrt(POW2(v3x-v2x) + POW2(v3y-v2y) + POW2(v3z-v2z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      // (3,4)
      edge = std::sqrt(POW2(v3x-v4x) + POW2(v3y-v4y) + POW2(v3z-v4z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      // (1,4)
      edge = std::sqrt(POW2(v1x-v4x) + POW2(v1y-v4y) + POW2(v1z-v4z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      // (5,8)
      edge = std::sqrt(POW2(v5x-v8x) + POW2(v5y-v8y) + POW2(v5z-v8z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      // (8,7)
      edge = std::sqrt(POW2(v8x-v7x) + POW2(v8y-v7y) + POW2(v8z-v7z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      // (7,6)
      edge = std::sqrt(POW2(v7x-v6x) + POW2(v7y-v6y) + POW2(v7z-v6z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      // (5,6)
      edge = std::sqrt(POW2(v5x-v6x) + POW2(v5y-v6y) + POW2(v5z-v6z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      // (5,1)
      edge = std::sqrt(POW2(v5x-v1x) + POW2(v5y-v1y) + POW2(v5z-v1z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      // (8,2)
      edge = std::sqrt(POW2(v8x-v2x) + POW2(v8y-v2y) + POW2(v8z-v2z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      // (7,3)
      edge = std::sqrt(POW2(v7x-v3x) + POW2(v7y-v3y) + POW2(v7z-v3z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      // (4,6)
      edge = std::sqrt(POW2(v4x-v6x) + POW2(v4y-v6y) + POW2(v4z-v6z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      num_edges = 12.0;
      break;
    }
    case Tetra:
    {
      const mt_int v1 = con[0], v2 = con[1], v3 = con[2], v4 = con[3];
      mt_real edge;

      const mt_real v1x = xyz[v1*3+0], v1y = xyz[v1*3+1], v1z = xyz[v1*3+2];
      const mt_real v2x = xyz[v2*3+0], v2y = xyz[v2*3+1], v2z = xyz[v2*3+2];
      const mt_real v3x = xyz[v3*3+0], v3y = xyz[v3*3+1], v3z = xyz[v3*3+2];
      const mt_real v4x = xyz[v4*3+0], v4y = xyz[v4*3+1], v4z = xyz[v4*3+2];

      // (1,2)
      edge = sqrt(POW2(v1x - v2x) + POW2(v1y - v2y) + POW2(v1z - v2z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;
      // (1,3)
      edge = sqrt(POW2(v1x - v3x) + POW2(v1y - v3y) + POW2(v1z - v3z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;
      // (2,3)
      edge = sqrt(POW2(v2x - v3x) + POW2(v2y - v3y) + POW2(v2z - v3z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;
      // (2,4)
      edge = sqrt(POW2(v2x - v4x) + POW2(v2y - v4y) + POW2(v2z - v4z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;
      // (1,4)
      edge = sqrt(POW2(v1x - v4x) + POW2(v1y - v4y) + POW2(v1z - v4z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;
      // (3,4)
      edge += sqrt(POW2(v3x - v4x) + POW2(v3y - v4y) + POW2(v3z - v4z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      num_edges = 6.0;
      break;
    }

    case Tri:
    {
      const mt_int v1 = con[0], v2 = con[1], v3 = con[2];
      mt_real edge;

      const mt_real v1x = xyz[v1*3+0], v1y = xyz[v1*3+1], v1z = xyz[v1*3+2];
      const mt_real v2x = xyz[v2*3+0], v2y = xyz[v2*3+1], v2z = xyz[v2*3+2];
      const mt_real v3x = xyz[v3*3+0], v3y = xyz[v3*3+1], v3z = xyz[v3*3+2];

      // (1,2)
      edge = sqrt(POW2(v1x - v2x) + POW2(v1y - v2y) + POW2(v1z - v2z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;
      // (1,3)
      edge = sqrt(POW2(v1x - v3x) + POW2(v1y - v3y) + POW2(v1z - v3z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;
      // (2,3)
      edge = sqrt(POW2(v2x - v3x) + POW2(v2y - v3y) + POW2(v2z - v3z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      num_edges = 3.0;
      break;
    }

    case Line:
    {
      const mt_int v1 = con[0], v2 = con[1];
      mt_real edge;

      const mt_real v1x = xyz[v1*3+0], v1y = xyz[v1*3+1], v1z = xyz[v1*3+2];
      const mt_real v2x = xyz[v2*3+0], v2y = xyz[v2*3+1], v2z = xyz[v2*3+2];

      // (1,2)
      edge = sqrt(POW2(v1x - v2x) + POW2(v1y - v2y) + POW2(v1z - v2z));
      if(min_edge > edge) min_edge = edge;
      if(max_edge < edge) max_edge = edge;
      avg_edge += edge;

      num_edges = 1.0;
      break;
    }

    default:
      std::cerr << "element_edges_stats error: element not yet implemented. Aborting!" << std::endl;
      exit(EXIT_FAILURE);
  }

  avg_edge /= num_edges;
}

mt_real avrg_edgelength_estimate(const mt_meshdata & mesh, bool full)
{
  check_nonzero(mesh.e2n_dsp.size(), __func__);
  mt_real tavg = 0;

  if(full) {
    #ifdef OPENMP
    #pragma omp parallel for schedule(guided) reduction(+:tavg)
    #endif
    for(size_t eidx = 0; eidx < mesh.e2n_cnt.size(); eidx++) {
      mt_real avg, min, max;
      element_edges_stats(mesh.etype[eidx], mesh.e2n_con.data() + mesh.e2n_dsp[eidx],
          mesh.xyz.data(), min, max, avg);
      tavg += avg;
    }

    tavg /= mt_real(mesh.e2n_cnt.size());
  }
  else {
    const int numsamples = 10;

    for(int s=0; s < numsamples; s++) {
      mt_int eidx = mt_int(drand48() * (mesh.e2n_cnt.size() - 1));
      mt_real avg, min, max;
      element_edges_stats(mesh.etype[eidx], mesh.e2n_con.data() + mesh.e2n_dsp[eidx],
          mesh.xyz.data(), min, max, avg);
      tavg += avg;
    }
    tavg /= numsamples;
  }

  return tavg;
}

mt_real min_edgelength_estimate(const mt_meshdata & mesh, bool full, bool nonzeromin)
{
  //check_nonzero(mesh.e2n_dsp.size(), __func__);
  mt_real tmin = 1e100;
  //GET THE NON-ZERO MIN!
  if(full) {
    #ifdef OPENMP
    #pragma omp parallel for schedule(guided) reduction(min:tmin)
    #endif
    for(size_t eidx = 0; eidx < mesh.e2n_cnt.size(); eidx++) {
      mt_real avg, min, max;
      element_edges_stats(mesh.etype[eidx], mesh.e2n_con.data() + mesh.e2n_dsp[eidx],
          mesh.xyz.data(), min, max, avg);
      if(nonzeromin) {
        if(tmin > min && min > 0.) {
          tmin = min;
        }
      } else {
        if(tmin > min) {
          tmin = min;
        }
      }
    }
  }
  else {
    const int numsamples = 10;

    for(int s=0; s < numsamples; s++) {
      mt_int eidx = mt_int(drand48() * (mesh.e2n_cnt.size() - 1));
      mt_real avg, min, max;
      element_edges_stats(mesh.etype[eidx], mesh.e2n_con.data() + mesh.e2n_dsp[eidx],
          mesh.xyz.data(), min, max, avg);
      if(tmin > min) tmin = min;
    }
  }

  return tmin;
}


void generate_surfmesh_sizing_field(const mt_meshdata & surfmesh, mt_vector<mt_real> & sizes)
{
  assert(surfmesh.e2n_dsp.size() > 0);
  const size_t numnodes = surfmesh.xyz.size() / 3;
  const mt_real * xyz = surfmesh.xyz.data();

  sizes.resize(numnodes);

  #ifdef OPENMP
  #pragma omp parallel for schedule(guided)
  #endif
  for(size_t i=0; i < numnodes; i++) {
    vec3r np(xyz + i*3);
    mt_real local_edge_avg = 0.;
    mt_real numadd = 0.0;
    mt_int start = surfmesh.n2n_dsp[i], stop = start + surfmesh.n2n_cnt[i];
    for(mt_int k=start; k < stop; k++) {
      mt_int cidx = surfmesh.n2n_con[k];
      vec3r p(xyz + cidx*3);

      mt_real edgelen = (np - p).length();
      local_edge_avg += edgelen;
      numadd += 1.0;
    }
    numadd -= 1.0;
    if(numadd > 0)
      local_edge_avg /= numadd;
    sizes[i] = local_edge_avg;
  }
}

void identify_sharp_edge_nodes(const mt_meshdata & surfmesh,
                               const mt_vector<mt_real> & xyz,
                               const mt_vector<mt_int> & vtxlist,
                               MT_USET<mt_int> & edge_vtx,
                               mt_real thr)
{
  // edge_thr is assumed in degrees, we convert it back to a length
  thr = cos(thr * (MT_PI/ 180.0));

  for(size_t idx = 0; idx < vtxlist.size(); idx++)
  {
    mt_int nidx = vtxlist[idx];
    mt_point<mt_real> nref(0, 0, 0);
    bool sharp = false;

    mt_int start = surfmesh.n2e_dsp[nidx], end = start + surfmesh.n2e_cnt[nidx];
    for(mt_int i=start; i<end; i++)
    {
      mt_int cele = surfmesh.n2e_con[i], dsp = surfmesh.e2n_dsp[cele];
      mt_int v1 = surfmesh.e2n_con[dsp+0], v2 = surfmesh.e2n_con[dsp+1], v3 = surfmesh.e2n_con[dsp+2];
      mt_point<mt_real> n = triangle_normal(v1, v2, v3, xyz);

      if(nref.length2() > 0) {
        // a refrerence normal has been already set. we can compute a current normal
        // and do the angle check
        mt_real pro = nref.scaProd(n);
        if(fabs(pro) < thr) {
          sharp = true;
          break;
        }
      }
      else {
        // the refrerence edge has not been set. we set it now
        nref = n;
      }
    }

    if(sharp) {
      edge_vtx.insert(nidx);

      // uncomment this if you want to add the neighbours of a sharp edge to the
      // sharp edges set
      #if 0
      start = surfmesh.n2n_dsp[nidx];
      end = start + surfmesh.n2n_cnt[nidx];

      for(int i=start; i<end; i++)
        edge_vtx.insert(surfmesh.n2n_con[i]);
      #endif
    }
  }
}

void identify_sharp_edges(const mt_meshdata & surfmesh,
                          const mt_vector<mt_real> & xyz,
                          const mt_mapping<mt_int> & ele2edge,
                          const MT_MAP<tuple<mt_int>, mt_int> & edges,
                          mt_real thr,
                          MT_USET<mt_int> & sharp_edges)
{
  // edge_thr is assumed in degrees, we convert it back to a length
  thr = cos(thr * (MT_PI/ 180.0));

  for(auto eit = edges.begin(); eit != edges.end(); ++eit)
  {
    mt_int edge_idx = eit->second;
    vec3r nref(0, 0, 0);
    bool sharp = false;

    mt_int start = ele2edge.bwd_dsp[edge_idx], end = start + ele2edge.bwd_cnt[edge_idx];
    for(mt_int i=start; i<end; i++)
    {
      mt_int eidx = ele2edge.bwd_con[i], elem_dsp = surfmesh.e2n_dsp[eidx];

      mt_int v1 = surfmesh.e2n_con[elem_dsp+0], v2 = surfmesh.e2n_con[elem_dsp+1],
             v3 = surfmesh.e2n_con[elem_dsp+2];
      vec3r  n  = triangle_normal(v1, v2, v3, xyz);

      if(nref.length2() > 0) {
        // a refrerence normal has been already set. we can compute a current normal
        // and do the angle check
        mt_real pro = nref.scaProd(n);
        if(fabs(pro) < thr) {
          sharp = true;
          break;
        }
      }
      else {
        // the refrerence edge has not been set. we set it now
        nref = n;
      }
    }

    if(sharp)
      sharp_edges.insert(edge_idx);
  }
}

void identify_surface_border_nodes(const mt_meshdata & surfmesh, MT_USET<mt_int> & border)
{
  mt_mapping<mt_int> ele2edge;
  MT_MAP<tuple<mt_int>, mt_int> edges;
  compute_edges(surfmesh, ele2edge, edges);
  ele2edge.transpose();

  for(auto it = edges.begin(); it != edges.end(); ++it) {
    mt_int edge_idx = it->second;
    if(ele2edge.bwd_cnt[edge_idx] == 1) {
      // edge connected to only one element -> its a border edge of an open surface
      border.insert(it->first.v1);
      border.insert(it->first.v2);
    }
  }
}


void add_interface_to_splitlist(const mt_meshdata & mesh,
                                const mt_meshgraph & mgA,
                                const mt_meshgraph & mgB,
                                mt_int & Nidx,
                                mt_vector<split_item> & splitlist)
{
  mt_vector<mt_int> nodA(mgA.e2n_con), nodB(mgB.e2n_con);
  binary_sort(nodA); unique_resize(nodA);
  binary_sort(nodB); unique_resize(nodB);

  MT_USET<mt_int> eidxA;
  eidxA.insert(mgA.eidx.begin(), mgA.eidx.end());

  int rsize = nodA.size() > nodB.size() ? nodA.size() : nodB.size();
  mt_vector<mt_int> ifc(rsize);

  // compute set of interface nodes via set intersection
  {
    mt_int* e = std::set_intersection(nodA.begin(), nodA.end(),
        nodB.begin(), nodB.end(), ifc.begin());
    ifc.resize(e - ifc.begin());
  }

  splitlist.reserve(splitlist.size() + ifc.size());
  split_item itm;

  // iterate over interface
  for(mt_int nidx : ifc) {
    mt_int start = mesh.n2e_dsp[nidx], end = start + mesh.n2e_cnt[nidx];
    // iterate over elements connected to interface node nidx
    for(int k=start; k < end; k++) {
      mt_int eidx = mesh.n2e_con[k];
      // if element tag is in tagsA insert it into the final splitlist
      if(eidxA.count(eidx)) {
        itm.elemIdx = eidx;
        itm.oldIdx  = nidx;
        itm.newIdx  = Nidx;
        splitlist.push_back(itm);
      }
    }
    Nidx++;
  }
}

void add_surface_to_splitlist(const mt_meshdata & mesh,
                              const mt_meshdata & surf,
                              const MT_USET<mt_int> & vol_eidx,
                              mt_int & Nidx,
                              mt_vector<split_item> & splitlist)
{
  mt_vector<mt_int> ifc(surf.e2n_con);
  binary_sort(ifc); unique_resize(ifc);

  split_item itm;

  for(mt_int nidx : ifc) {
    mt_int start = mesh.n2e_dsp[nidx], end = start + mesh.n2e_cnt[nidx];
    // iterate over elements connected to interface node nidx
    for(int k=start; k < end; k++) {
      mt_int eidx = mesh.n2e_con[k];
      // if element tag is in vol_eidx insert it into the final splitlist
      if(vol_eidx.count(eidx)) {
        itm.elemIdx = eidx;
        itm.oldIdx  = nidx;
        itm.newIdx  = Nidx;
        splitlist.push_back(itm);
      }
    }
    Nidx++;
  }
}
