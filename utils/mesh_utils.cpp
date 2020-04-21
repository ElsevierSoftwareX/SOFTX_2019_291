/**
* @file mesh_utils.cpp
* @brief General mesh utility funcions
* @author Aurel Neic
* @version
* @date 2017-08-16
*/

#include <omp.h>

#include "mt_utils.h"
#include "dense_mat.hpp"


void reindex_nodes(mt_meshdata & mesh,
                   mt_vector<mt_int> & nod_out,
                   bool verbose)
{
  mt_vector<mt_int> & con = mesh.e2n_con;
  mt_vector<mt_real> & xyz = mesh.xyz;

  PROGRESS<short>* prg = NULL;
  if(verbose)
    prg = new PROGRESS<short>(5, "Reindexing progress: ");

  // create nod_out
  nod_out.assign(con.begin(), con.end());

  if(prg) prg->next();

  binary_sort(nod_out);
  unique_resize(nod_out);

  if(prg) prg->next();

  // extract new list of points
  mt_vector<mt_real> old_xyz(xyz);
  for(size_t i=0; i<nod_out.size(); i++)
  {
    mt_int n = nod_out[i];
    xyz[i*3+0] = old_xyz[n*3+0];
    xyz[i*3+1] = old_xyz[n*3+1];
    xyz[i*3+2] = old_xyz[n*3+2];
  }
  xyz.resize(nod_out.size()*3);

  if(prg) prg->next();

  // remap con to indices in nod_out -------------------------------------
  map_glob2loc(nod_out, con);

  if(prg) prg->next();

  // reallocate
  xyz.reallocate();

  if(verbose) {
    prg->finish();
    delete prg;
  }
}

void reindex_nodes(mt_meshdata & mesh, bool verbose)
{
  mt_vector<mt_int> & con = mesh.e2n_con;
  mt_vector<mt_real> & xyz = mesh.xyz;

  PROGRESS<short>* prg = NULL;
  if(verbose)
    prg = new PROGRESS<short>(6, "Reindexing progress: ");

  // create nod_out
  MT_MAP<mt_int,mt_int> g2l;
  {
    for(const mt_int & c : con)
      g2l[c] = 0;

    if(prg) prg->next();
    g2l.sort();

    if(prg) prg->next();
    auto it = g2l.begin();
    for(size_t i=0; i < g2l.size(); ++i, ++it)
      it->second = i;
  }

  if(prg) prg->next();

  // extract new list of points
  mt_vector<mt_real> old_xyz(xyz);
  for(auto & m : g2l)
  {
    mt_int glb = m.first;
    mt_int loc = m.second;
    xyz[loc*3+0] = old_xyz[glb*3+0];
    xyz[loc*3+1] = old_xyz[glb*3+1];
    xyz[loc*3+2] = old_xyz[glb*3+2];
  }
  xyz.resize(g2l.size()*3);

  if(prg) prg->next();

  // remap con to indices in nod_out -------------------------------------
  map_glob2loc(g2l, con);

  if(prg) prg->next();

  // reallocate
  xyz.reallocate();

  if(verbose) {
    prg->finish();
    delete prg;
  }
}


void restrict_elements(const mt_vector<bool> & keep,
                       mt_meshdata & mesh,
                       mt_vector<mt_int> & eidx_out,
                       bool verbose)
{
  mt_vector<mt_int> & cnt = mesh.e2n_cnt;
  mt_vector<mt_int> & con = mesh.e2n_con;
  mt_vector<mt_int> & tag = mesh.etags;

  mt_vector<elem_t> & type = mesh.etype;

  short nlon = (mesh.lon.size() / 6) == mesh.e2n_cnt.size() ? 6 : 3;

  size_t widx=0;
  size_t e_ridx=0, e_widx=0;
  size_t numelems = cnt.size();
  eidx_out.resize(numelems);

  PROGRESS<size_t>* prg = NULL;

  if(verbose)
    prg = new PROGRESS<size_t>(numelems, "Mesh elements restriction: ");

  for(size_t i=0; i<numelems; i++)
  {
    if( keep[i] )
    {
      // element index
      eidx_out[widx] = i;
      // element size
      cnt[widx] = cnt[i];
      // tag
      tag[widx] = tag[i];
      // type
      type[widx] = type[i];
      // fibers
      mt_real* lp = mesh.lon.data();
      if(lp)
        for(int j=0; j<nlon; j++) lp[widx*nlon + j] = lp[i*nlon + j];

      widx++;

      // element connectivity
      for(int j=0; j<cnt[i]; j++)
        con[e_widx++] = con[e_ridx++];

    }
    else
      e_ridx += cnt[i];

    if(prg) prg->next();
  }
  eidx_out.resize(widx);
  cnt.resize(widx);
  tag.resize(widx);
  type.resize(widx);
  mesh.lon.resize(widx*nlon);
  con.resize(e_widx);

  // reallocate
  eidx_out.reallocate();
  cnt.reallocate();
  con.reallocate();
  tag.reallocate();
  type.reallocate();

  if(verbose) {
    prg->finish();
    delete prg;
  }
}

void extract_tagged_meshgraph(const MT_USET<mt_int> & tags,
                              const mt_meshdata & mesh,
                              mt_meshgraph & graph)
{
  size_t numelem = mesh.e2n_cnt.size();
  size_t num_extr_elem = 0, num_extr_entr = 0;

  for(size_t i=0; i<numelem; i++) {
    if(tags.count(mesh.etags[i])) {
      num_extr_elem++;
      num_extr_entr += mesh.e2n_cnt[i];
    }
  }

  graph.e2n_cnt.resize(num_extr_elem);
  graph.etype  .resize(num_extr_elem);
  graph.e2n_con.resize(num_extr_entr);
  graph.eidx.resize(num_extr_elem);

  for(size_t ridx_cnt=0, ridx_con=0, widx_cnt=0, widx_con=0; ridx_cnt<numelem; ridx_cnt++) {
    if(tags.count(mesh.etags[ridx_cnt])) {
      graph.e2n_cnt[widx_cnt] = mesh.e2n_cnt[ridx_cnt];
      graph.etype[widx_cnt]   = mesh.etype[ridx_cnt];
      graph.eidx[widx_cnt]    = ridx_cnt;

      for(mt_int j=0; j<graph.e2n_cnt[widx_cnt]; j++)
        graph.e2n_con[widx_con++] = mesh.e2n_con[ridx_con++];

      widx_cnt++;
    }
    else
      ridx_con += mesh.e2n_cnt[ridx_cnt];
  }
}

void extract_meshgraph(const mt_vector<mt_int> & eidx,
                       const mt_meshdata & mesh,
                       mt_meshgraph & graph)
{
  assert(mesh.e2n_dsp.size() > 0);

  size_t num_extr_elem = eidx.size(), num_extr_entr = 0;

  for(const mt_int & e : eidx)
    num_extr_entr += mesh.e2n_cnt[e];

  graph.e2n_cnt.resize(num_extr_elem);
  graph.e2n_con.resize(num_extr_entr);
  graph.eidx.resize(num_extr_elem);
  graph.etype.resize(num_extr_elem);

  size_t widx_cnt = 0, widx_con = 0;
  for(const mt_int & ridx : eidx)
  {
    graph.e2n_cnt[widx_cnt] = mesh.e2n_cnt[ridx];
    graph.etype  [widx_cnt] = mesh.etype[ridx];
    graph.eidx   [widx_cnt] = ridx;

    mt_int ridx_con = mesh.e2n_dsp[ridx];
    for(mt_int j=0; j<graph.e2n_cnt[widx_cnt]; j++)
      graph.e2n_con[widx_con++] = mesh.e2n_con[ridx_con++];

    widx_cnt++;
  }
}

void extract_mesh(const mt_vector<mt_int> & eidx,
                  const mt_meshdata & mesh,
                  mt_meshdata & outmesh)
{
  assert(mesh.e2n_dsp.size() > 0);

  size_t num_extr_elem = eidx.size(), num_extr_entr = 0;
  const short numfibvals = mesh.e2n_cnt.size() * 6 == mesh.lon.size() ? 6 : 3;

  for(const mt_int & e : eidx)
    num_extr_entr += mesh.e2n_cnt[e];

  outmesh.e2n_cnt.resize(num_extr_elem);
  outmesh.e2n_con.resize(num_extr_entr);
  outmesh.etype.resize(num_extr_elem);
  outmesh.etags.resize(num_extr_elem, 0);
  outmesh.lon.resize(num_extr_elem * numfibvals, 0.0);

  const bool have_tags = mesh.etags.size() > 0;
  const bool have_fibs = mesh.lon.size() > 0;

  size_t widx_cnt = 0, widx_con = 0;
  for(const mt_int & ridx : eidx)
  {
    outmesh.e2n_cnt[widx_cnt] = mesh.e2n_cnt[ridx];
    outmesh.etype  [widx_cnt] = mesh.etype  [ridx];
    if(have_tags)
      outmesh.etags  [widx_cnt] = mesh.etags  [ridx];

    mt_int ridx_con = mesh.e2n_dsp[ridx];
    for(mt_int j=0; j<outmesh.e2n_cnt[widx_cnt]; j++)
      outmesh.e2n_con[widx_con++] = mesh.e2n_con[ridx_con++];

    if(have_fibs)
      for(short j=0; j<numfibvals; j++)
        outmesh.lon[widx_cnt*numfibvals + j] = mesh.lon[ridx*numfibvals + j];

    widx_cnt++;
  }

  outmesh.xyz = mesh.xyz;
}

void insert_etags(struct mt_meshdata & mesh, const mt_vector<mt_int> & etags, const mt_vector<mt_int> & eidx)
{
  for(size_t i=0; i<eidx.size(); i++) mesh.etags[eidx[i]] = etags[i];
}

void insert_fibers(struct mt_meshdata & mesh, const mt_vector<mt_real> & lon, const mt_vector<mt_int> & eidx)
{
  bool twolon = (mesh.lon.size() / 6) == mesh.e2n_cnt.size() ? true : false;

  if(twolon) {
    for(size_t i=0; i<eidx.size(); i++) {
      const mt_real* read = lon.data() + i*6;
      mt_real* write = mesh.lon.data() + eidx[i]*6;
      write[0] = read[0];
      write[1] = read[1];
      write[2] = read[2];
      write[3] = read[3];
      write[4] = read[4];
      write[5] = read[5];
    }
  }
  else {
    for(size_t i=0; i<eidx.size(); i++) {
      const mt_real* read = lon.data() + i*3;
      mt_real* write = mesh.lon.data() + eidx[i]*3;
      write[0] = read[0];
      write[1] = read[1];
      write[2] = read[2];
    }
  }
}

void insert_points(struct mt_meshdata & mesh, const mt_vector<mt_real> & xyz, const mt_vector<mt_int> & nod)
{
  for(size_t i=0; i<nod.size(); i++) {
    const mt_real* read = xyz.data() + i*3;
    mt_real* write = mesh.xyz.data() + nod[i]*3;
    write[0] = read[0];
    write[1] = read[1];
    write[2] = read[2];
  }
}

void insert_surf_tri(mt_int n1, mt_int n2, mt_int n3, size_t eidx,
                     triple<mt_int> & surf,
                     tri_sele & sele,
                     MT_MAP<triple<mt_int>, tri_sele> & surfmap)
{
  sele.v1 = n1, sele.v2 = n2, sele.v3 = n3; sele.eidx = eidx;
  sortTriple(n1, n2, n3, surf.v1, surf.v2, surf.v3);

  auto it = surfmap.find(surf);
  if(it != surfmap.end()) surfmap.erase(it);
  else surfmap[surf] = sele;
}

void insert_surf_quad(mt_int n1, mt_int n2, mt_int n3, mt_int n4, size_t eidx,
                      mt_vector<mt_int> & buff,
                      quadruple<mt_int> & surf,
                      quad_sele & sele,
                      MT_MAP<quadruple<mt_int>, quad_sele> & surfmap)
{
  buff[0] = n1, buff[1] = n2, buff[2] = n3, buff[3] = n4;
  binary_sort(buff);

  sele.v1 = n1, sele.v2 = n2, sele.v3 = n3, sele.v4 = n4, sele.eidx = eidx;
  surf.v1 = buff[0], surf.v2 = buff[1], surf.v3 = buff[2], surf.v4 = buff[3];

  auto it = surfmap.find(surf);
  if(it != surfmap.end()) surfmap.erase(it);
  else surfmap[surf] = sele;
}

void insert_surf_tet(const mt_int* nod,
                     const size_t eidx,
                     MT_MAP<triple<mt_int>, tri_sele> & surfmap)
{
  // surfaces are defined according to carpmanual in *clockwise* order
  // surfaces are (2,3,1) , (1,4,2) , (2,4,3) , (1,3,4)
  struct triple<mt_int> surf;
  struct tri_sele  sele;

  mt_int n1 = nod[0], n2 = nod[1], n3 = nod[2], n4 = nod[3];

  insert_surf_tri(n2, n3, n1, eidx, surf, sele, surfmap);
  insert_surf_tri(n1, n4, n2, eidx, surf, sele, surfmap);
  insert_surf_tri(n2, n4, n3, eidx, surf, sele, surfmap);
  insert_surf_tri(n1, n3, n4, eidx, surf, sele, surfmap);
}

void insert_surf_pyr(const mt_int* nod, const size_t eidx, mt_vector<mt_int> & buff,
                     MT_MAP<triple<mt_int>, tri_sele> & surfmap,
                     MT_MAP<quadruple<mt_int>, quad_sele> & qsurfmap)
{
  // surfaces are defined according to carpmanual in *clockwise* order
  // surfaces are (1,5,2) , (2,5,3) , (3,5,4) , (4,5,1) , (1,2,3,4)
  struct triple<mt_int>    surf; struct quadruple<mt_int> qsurf;
  struct tri_sele  sele; struct quad_sele qsele;

  mt_int n1 = nod[0], n2 = nod[1], n3 = nod[2], n4 = nod[3], n5 = nod[4];

  insert_surf_tri(n1, n5, n2, eidx, surf, sele, surfmap);
  insert_surf_tri(n2, n5, n3, eidx, surf, sele, surfmap);
  insert_surf_tri(n3, n5, n4, eidx, surf, sele, surfmap);
  insert_surf_tri(n4, n5, n1, eidx, surf, sele, surfmap);
  insert_surf_quad(n1, n2, n3, n4, eidx, buff, qsurf, qsele, qsurfmap);
}

void insert_surf_pri(const mt_int* nod, const size_t eidx, mt_vector<mt_int> & buff,
                     MT_MAP<triple<mt_int>, tri_sele> & surfmap,
                     MT_MAP<quadruple<mt_int>, quad_sele> & qsurfmap)
{
  struct triple<mt_int> surf;
  struct quadruple<mt_int> qsurf;
  struct tri_sele  sele;
  struct quad_sele qsele;

  mt_int n1 = nod[0], n2 = nod[1], n3 = nod[2], n4 = nod[3], n5 = nod[4], n6 = nod[5];

  // surfaces are defined according to carpmanual in *clockwise* order
  // surfaces are (1,2,3) , (4,5,6) , (1,4,6,2) , (3,2,6,5) , (1,3,5,4)
  insert_surf_tri(n1, n2, n3, eidx, surf, sele, surfmap);
  insert_surf_tri(n4, n5, n6, eidx, surf, sele, surfmap);
  insert_surf_quad(n1, n4, n6, n2, eidx, buff, qsurf, qsele, qsurfmap);
  insert_surf_quad(n3, n2, n6, n5, eidx, buff, qsurf, qsele, qsurfmap);
  insert_surf_quad(n1, n3, n5, n4, eidx, buff, qsurf, qsele, qsurfmap);
}

void insert_surf_hex(const mt_int* nod, const size_t eidx, mt_vector<mt_int> & buff,
                     MT_MAP<quadruple<mt_int>, quad_sele> & surfmap)
{
  // surfaces are defined according to carpmanual in *clockwise* order
  // surfaces are (1,2,3,4) , (3,2,8,7) , (4,3,7,6) , (1,4,6,5) , (2,1,5,8), (6,7,8,5)
  struct quadruple<mt_int> qsurf;
  struct quad_sele qsele;

  mt_int n1 = nod[0], n2 = nod[1], n3 = nod[2], n4 = nod[3],
         n5 = nod[4], n6 = nod[5], n7 = nod[6], n8 = nod[7];

  insert_surf_quad(n1, n2, n3, n4, eidx, buff, qsurf, qsele, surfmap);
  insert_surf_quad(n3, n2, n8, n7, eidx, buff, qsurf, qsele, surfmap);
  insert_surf_quad(n4, n3, n7, n6, eidx, buff, qsurf, qsele, surfmap);
  insert_surf_quad(n1, n4, n6, n5, eidx, buff, qsurf, qsele, surfmap);
  insert_surf_quad(n2, n1, n5, n8, eidx, buff, qsurf, qsele, surfmap);
  insert_surf_quad(n6, n7, n8, n5, eidx, buff, qsurf, qsele, surfmap);
}

void compute_surface(const mt_vector<elem_t> & etype,
                     const mt_vector<mt_int> & e2n_cnt,
                     const mt_vector<mt_int> & e2n_con,
                     mt_meshdata & surfmesh,
                     const bool full_mesh)
{
  MT_MAP<triple<mt_int>, tri_sele>  tri_surf;
  MT_MAP<quadruple<mt_int>, quad_sele> quad_surf;
  mt_vector<mt_int> elem_orig;

  compute_surface(etype, e2n_cnt, e2n_con, tri_surf, quad_surf);

  if(full_mesh)
    surfmap_to_surfmesh(tri_surf, quad_surf, surfmesh, elem_orig);
  else
    surfmap_to_vector(tri_surf, quad_surf, surfmesh.e2n_con, elem_orig);
}

void compute_surface(const mt_vector<elem_t> & etype,
                     const mt_vector<mt_int> & e2n_cnt,
                     const mt_vector<mt_int> & e2n_con,
                     const mt_vector<mt_int> & ref_eidx,
                     MT_MAP<triple<mt_int>, tri_sele>     & tri_surf,
                     MT_MAP<quadruple<mt_int>, quad_sele> & quad_surf)
{
   const mt_int* nod = e2n_con.data();
   mt_vector<mt_int> qbuff(4); // buffer if quad nodes need to be sorted

   for(size_t eidx=0; eidx<e2n_cnt.size(); eidx++) {
     switch(etype[eidx]) {
       case Tri: {
         triple<mt_int> tri;
         tri_sele trie;
         insert_surf_tri(nod[0], nod[1], nod[2], ref_eidx[eidx], tri, trie, tri_surf);
         break;
       }
       case Quad: {
         quadruple<mt_int> qd;
         quad_sele qde;
         insert_surf_quad(nod[0], nod[1], nod[2], nod[3], ref_eidx[eidx], qbuff, qd, qde, quad_surf);
         break;
       }
       case Tetra:
         insert_surf_tet(nod, ref_eidx[eidx], tri_surf);
         break;
       case Pyramid:
         insert_surf_pyr(nod, ref_eidx[eidx], qbuff, tri_surf, quad_surf);
         break;
       case Prism:
         insert_surf_pri(nod, ref_eidx[eidx], qbuff, tri_surf, quad_surf);
         break;
       case Hexa:
         insert_surf_hex(nod, ref_eidx[eidx], qbuff, quad_surf);
         break;

       default:
         fprintf(stderr, "Error: Unsupported element in surface computation!\n");
         exit(1);
     }
     nod += e2n_cnt[eidx];
   }
}

void compute_surface_seq(const mt_meshdata & mesh,
                         const MT_USET<mt_int> & tags,
                         MT_MAP<triple<mt_int>, tri_sele>     & tri_surf,
                         MT_MAP<quadruple<mt_int>, quad_sele> & quad_surf)
{
  mt_vector<mt_int> qbuff(4); // buffer if quad nodes need to be sorted
  size_t nele = mesh.e2n_cnt.size();

  mt_vector<mt_int> dsp(nele);
  bucket_sort_offset(mesh.e2n_cnt, dsp);

  for(size_t eidx=0; eidx<nele; eidx++)
  {
    if(tags.count(mesh.etags[eidx])) {
      const mt_int * nod = mesh.e2n_con.data() + dsp[eidx];
      switch(mesh.etype[eidx]) {
        case Tri: {
          triple<mt_int> tri;
          tri_sele trie;
          insert_surf_tri(nod[0], nod[1], nod[2], eidx, tri, trie, tri_surf);
          break;
        }
        case Quad: {
          quadruple<mt_int> qd;
          quad_sele qde;
          insert_surf_quad(nod[0], nod[1], nod[2], nod[3], eidx, qbuff, qd, qde, quad_surf);
          break;
        }
        case Tetra:
          insert_surf_tet(nod, eidx, tri_surf);
          break;
        case Pyramid:
          insert_surf_pyr(nod, eidx, qbuff, tri_surf, quad_surf);
          break;
        case Prism:
          insert_surf_pri(nod, eidx, qbuff, tri_surf, quad_surf);
          break;
        case Hexa:
          insert_surf_hex(nod, eidx, qbuff, quad_surf);
          break;

        default:
          fprintf(stderr, "Error: Unsupported element in surface computation!\n");
          exit(1);
      }
    }
  }
}

void compute_surface_parallel(const mt_meshdata & mesh,
                              const MT_USET<mt_int> & tags,
                              MT_MAP<triple<mt_int>, tri_sele>     & tri_surf,
                              MT_MAP<quadruple<mt_int>, quad_sele> & quad_surf)
{
  size_t nele = mesh.e2n_cnt.size();
  mt_vector<mt_int> dsp(nele);
  bucket_sort_offset(mesh.e2n_cnt, dsp);

  /*
   * The parallelization idea is to split the mesh along its longest axis into
   * nthreads number of blocks. Each thread computes the surface of its own block.
   * In the end, the threads merge their local surfaces into the final surface.
   *
   */

  bbox box;
  generate_bbox(mesh.xyz, box);

  #ifdef OPENMP
  #pragma omp parallel
  #endif
  {
    mt_vector<mt_int> qbuff(4); // buffer if quad nodes need to be sorted

    bboxAxis longest_axis     = get_longest_axis(box);
    mt_real  longest_axis_len = get_longest_axis_length(box);

    #ifdef OPENMP
    int tid = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
    #else
    int tid = 0;
    int nthreads = 1;
    #endif

    mt_real blocksize = (longest_axis_len / nthreads) * 1.01;

    MT_MAP<triple<mt_int>, tri_sele>     loc_tri_surf;
    MT_MAP<quadruple<mt_int>, quad_sele> loc_quad_surf;

    for(size_t eidx=0; eidx<nele; eidx++)
    {
      if(tags.count(mesh.etags[eidx])) {
        const mt_int * nod = mesh.e2n_con.data() + dsp[eidx];
        vec3r ctr = barycenter(mesh.e2n_cnt[eidx], nod, mesh.xyz.data());

        mt_real coord;
        switch(longest_axis) {
          default: break;
          case X: coord = ctr.x - box.bounds[0].x; break;
          case Y: coord = ctr.y - box.bounds[0].y; break;
          case Z: coord = ctr.z - box.bounds[0].z; break;
        }

        int associated_thread = coord / blocksize;

        if(associated_thread == tid) {
          switch(mesh.etype[eidx]) {
            case Tri: {
              triple<mt_int> tri;
              tri_sele trie;
              insert_surf_tri(nod[0], nod[1], nod[2], eidx, tri, trie, loc_tri_surf);
              break;
            }
            case Quad: {
              quadruple<mt_int> qd;
              quad_sele qde;
              insert_surf_quad(nod[0], nod[1], nod[2], nod[3], eidx, qbuff, qd, qde, loc_quad_surf);
              break;
            }
            case Tetra:
              insert_surf_tet(nod, eidx, loc_tri_surf);
              break;
            case Pyramid:
              insert_surf_pyr(nod, eidx, qbuff, loc_tri_surf, loc_quad_surf);
              break;
            case Prism:
              insert_surf_pri(nod, eidx, qbuff, loc_tri_surf, loc_quad_surf);
              break;
            case Hexa:
              insert_surf_hex(nod, eidx, qbuff, loc_quad_surf);
              break;

            default:
              fprintf(stderr, "Error: Unsupported element in surface computation!\n");
              exit(1);
          }
        }
      }
    }

    #ifdef OPENMP
    #pragma omp critical
    #endif
    {
      for(auto it = loc_tri_surf.begin(); it != loc_tri_surf.end(); ++it)
      {
        auto fit = tri_surf.find(it->first);
        if(fit != tri_surf.end())
          tri_surf.erase(fit);
        else
          tri_surf.insert(*it);
      }

      for(auto it = loc_quad_surf.begin(); it != loc_quad_surf.end(); ++it)
      {
        auto fit = quad_surf.find(it->first);
        if(fit != quad_surf.end())
          quad_surf.erase(fit);
        else
          quad_surf.insert(*it);
      }
    }
  }
}

void compute_surface(const mt_meshdata & mesh,
                     const MT_USET<mt_int> & tags,
                     MT_MAP<triple<mt_int>, tri_sele>     & tri_surf,
                     MT_MAP<quadruple<mt_int>, quad_sele> & quad_surf)
{
#ifdef OPENMP
  compute_surface_parallel(mesh, tags, tri_surf, quad_surf);
#else
  compute_surface_seq(mesh, tags, tri_surf, quad_surf);
#endif
}

void compute_surface(const mt_meshdata & mesh,
                     const MT_USET<mt_int> & tags,
                     mt_meshdata & surfmesh,
                     const bool full_mesh)
{
  MT_MAP<triple<mt_int>, tri_sele>  tri_surf;
  MT_MAP<quadruple<mt_int>, quad_sele> quad_surf;
  mt_vector<mt_int> elem_orig;

  compute_surface(mesh, tags, tri_surf, quad_surf);

  if(full_mesh)
    surfmap_to_surfmesh(tri_surf, quad_surf, surfmesh, elem_orig);
  else
    surfmap_to_vector(tri_surf, quad_surf, surfmesh.e2n_con, elem_orig);
}

void compute_surface(const mt_vector<elem_t> & etype,
                     const mt_vector<mt_int> & e2n_cnt,
                     const mt_vector<mt_int> & e2n_con,
                     MT_MAP<triple<mt_int>, tri_sele>     & tri_surf,
                     MT_MAP<quadruple<mt_int>, quad_sele> & quad_surf)
{
  mt_vector<mt_int> eidx(etype.size());
  for(size_t i=0; i<eidx.size(); i++) eidx[i] = i;

  compute_surface(etype, e2n_cnt, e2n_con, eidx, tri_surf, quad_surf);
}

void surfmap_to_vector(const MT_MAP<struct triple<mt_int>, struct tri_sele> & tri,
                       const MT_MAP<struct quadruple<mt_int>, struct quad_sele> & quad,
                       mt_vector<mt_int> & surf_con,
                       mt_vector<mt_int> & elem_origin)
{
   typename MT_MAP<struct triple<mt_int>, struct tri_sele>::const_iterator tri_it = tri.begin();
   typename MT_MAP<struct quadruple<mt_int>, struct quad_sele>::const_iterator quad_it = quad.begin();

   surf_con.resize(tri.size()*3 + quad.size()*6);
   elem_origin.resize(tri.size() + quad.size()*2);

   mt_int sidx = 0, eidx = 0;

   while(tri_it != tri.end()) {
     const struct tri_sele& val = tri_it->second;
     surf_con[sidx+0] = val.v1;
     surf_con[sidx+1] = val.v2;
     surf_con[sidx+2] = val.v3;

     elem_origin[eidx++] = tri_it->second.eidx;

     sidx+=3;
     ++tri_it;
   }
   while(quad_it != quad.end()) {
     const struct quad_sele& val = quad_it->second;
     surf_con[sidx+0] = val.v1;
     surf_con[sidx+1] = val.v2;
     surf_con[sidx+2] = val.v3;

     elem_origin[eidx++] = quad_it->second.eidx;

     surf_con[sidx+3] = val.v1;
     surf_con[sidx+4] = val.v3;
     surf_con[sidx+5] = val.v4;

     elem_origin[eidx++] = quad_it->second.eidx;

     sidx+=6;
     ++quad_it;
   }
}

void surfmap_to_surfmesh(const MT_MAP<struct triple<mt_int>, struct tri_sele> & tri,
                         const MT_MAP<struct quadruple<mt_int>, struct quad_sele> & quad,
                         mt_meshdata & surf_mesh,
                         mt_vector<mt_int> & elem_origin)
{
  surf_mesh.e2n_cnt.assign(tri.size(), 3);
  surf_mesh.etype  .assign(tri.size(), Tri);
  surf_mesh.etags  .assign(tri.size(), mt_int(0));

  surf_mesh.e2n_cnt.resize(tri.size() + quad.size(), 4);
  surf_mesh.etype  .resize(tri.size() + quad.size(), Quad);
  surf_mesh.etags  .resize(tri.size() + quad.size(), mt_int(0));

  elem_origin      .resize(surf_mesh.e2n_cnt.size());
  surf_mesh.e2n_con.resize(tri.size()*3 + quad.size()*4);

  auto tri_it = tri.begin();
  for(size_t k=0; k<tri.size(); k++) {
    surf_mesh.e2n_con[k*3+0] = tri_it->second.v1;
    surf_mesh.e2n_con[k*3+1] = tri_it->second.v2;
    surf_mesh.e2n_con[k*3+2] = tri_it->second.v3;

    elem_origin[k] = tri_it->second.eidx;

    ++tri_it;
  }

  auto quad_it = quad.begin();
  for(size_t k=0, dsp = tri.size()*3; k<quad.size(); k++) {
    surf_mesh.e2n_con[dsp + k*4+0] = quad_it->second.v1;
    surf_mesh.e2n_con[dsp + k*4+1] = quad_it->second.v2;
    surf_mesh.e2n_con[dsp + k*4+2] = quad_it->second.v3;
    surf_mesh.e2n_con[dsp + k*4+3] = quad_it->second.v4;

    elem_origin[tri.size()+k] = quad_it->second.eidx;

    ++quad_it;
  }
}

void vector_to_surfmap(const mt_vector<mt_int> & surf_con,
                       MT_MAP<struct triple<mt_int>, struct tri_sele> & trimap)
{
  size_t num_tri = surf_con.size() / 3;
  triple<mt_int> stri;
  tri_sele       tri;

  for(size_t i=0; i<num_tri; i++) {
    tri.v1 = surf_con[i*3+0];
    tri.v2 = surf_con[i*3+1];
    tri.v3 = surf_con[i*3+2];
    tri.eidx = i;

    sortTriple(tri.v1, tri.v2, tri.v3, stri.v1, stri.v2, stri.v3);
    trimap[stri] = tri;
  }
}

void surfmesh_to_surfmap(const mt_meshdata & surfmesh,
                         const mt_vector<mt_int> & eidx,
                         MT_MAP<triple<mt_int>, tri_sele> & tri,
                         MT_MAP<quadruple<mt_int>, quad_sele> & quad)
{
  size_t numele = surfmesh.e2n_cnt.size();
  triple<mt_int> stri; quadruple<mt_int> squad;
  tri_sele tri_ele; quad_sele quad_ele;
  mt_vector<mt_int> qbuff(4); // buffer if quad nodes need to be sorted

  check_nonzero(surfmesh.etype.size(), __func__);
  const mt_int* nod = surfmesh.e2n_con.data();

  for(size_t i=0; i<numele; i++)
  {
    switch(surfmesh.etype[i]) {
      case Tri: {
        insert_surf_tri(nod[0], nod[1], nod[2], eidx[i], stri, tri_ele, tri);
        break;
      }

      case Quad: {
        insert_surf_quad(nod[0], nod[1], nod[2], nod[3], eidx[i], qbuff, squad, quad_ele, quad);
        break;
      }

      default:
        fprintf(stderr, "%s error: Non-surface element!\n", __func__);
        exit(1);
    }

    nod += surfmesh.e2n_cnt[i];
  }
}


void surface_difference(MT_MAP<triple<mt_int>, tri_sele> & lhs_tri,
                        MT_MAP<quadruple<mt_int>, quad_sele> & lhs_quad,
                        const MT_MAP<triple<mt_int>, tri_sele> & rhs_tri,
                        const MT_MAP<quadruple<mt_int>, quad_sele> & rhs_quad)
{
   typename MT_MAP<triple<mt_int>, tri_sele>::const_iterator tri_it = rhs_tri.begin();
   typename MT_MAP<quadruple<mt_int>, quad_sele>::const_iterator quad_it = rhs_quad.begin();

   // remove elements of rhs from lhs
   while(tri_it != rhs_tri.end()) {
     typename MT_MAP<triple<mt_int>, tri_sele>::iterator it = lhs_tri.find(tri_it->first);
     if(it != lhs_tri.end()) lhs_tri.erase(it);
     ++tri_it;
   }
   while(quad_it != rhs_quad.end()) {
     typename MT_MAP<quadruple<mt_int>, quad_sele>::iterator it = lhs_quad.find(quad_it->first);
     if(it != lhs_quad.end()) lhs_quad.erase(it);
     ++quad_it;
   }
}

void surface_intersection(MT_MAP<triple<mt_int>, tri_sele> & lhs_tri,
                          MT_MAP<quadruple<mt_int>, quad_sele> & lhs_quad,
                          const MT_MAP<triple<mt_int>, tri_sele> & rhs_tri,
                          const MT_MAP<quadruple<mt_int>, quad_sele> & rhs_quad)
{
   typename MT_MAP<triple<mt_int>, tri_sele>::iterator tri_it = lhs_tri.begin();
   typename MT_MAP<quadruple<mt_int>, quad_sele>::iterator quad_it = lhs_quad.begin();

   // iterate over lhs and remove elements which do not also exist in rhs
   while(tri_it != lhs_tri.end()) {
     if(rhs_tri.count(tri_it->first) == 0) {
       typename MT_MAP<triple<mt_int>, tri_sele>::iterator it = tri_it;
       ++tri_it;
       lhs_tri.erase(it);
     }
     else ++tri_it;
   }
   while(quad_it != lhs_quad.end()) {
     if(rhs_quad.count(quad_it->first) == 0) {
       typename MT_MAP<quadruple<mt_int>, quad_sele>::iterator it = quad_it;
       ++quad_it;
       lhs_quad.erase(it);
     }
     else ++quad_it;
   }
}

void surface_union(MT_MAP<triple<mt_int>, tri_sele> & lhs_tri,
                   MT_MAP<quadruple<mt_int>, quad_sele> & lhs_quad,
                   const MT_MAP<triple<mt_int>, tri_sele> & rhs_tri,
                   const MT_MAP<quadruple<mt_int>, quad_sele> & rhs_quad)
{
  typename MT_MAP<triple<mt_int>, tri_sele>::const_iterator     tri_it  = rhs_tri.begin();
  typename MT_MAP<quadruple<mt_int>, quad_sele>::const_iterator quad_it = rhs_quad.begin();

  // iterate over rhs and add its elements to lhs
  while(tri_it != rhs_tri.end()) {
    lhs_tri.insert(*tri_it);
    ++tri_it;
  }

  while(quad_it != rhs_quad.end()) {
    lhs_quad.insert(*quad_it);
    ++quad_it;
  }
}

void generate_nbc_data(const mt_meshdata & mesh,
                       const mt_meshdata & surf,
                       const mt_vector<mt_int> & elem_orig,
                       struct nbc_data & nbc)
{
  size_t ssize = surf.e2n_cnt.size();
  nbc.eidx.resize(ssize); nbc.sp_vtx.resize(ssize); nbc.tag.resize(ssize);

  for(size_t i = 0; i<ssize; i++) {
    mt_int mesh_eidx = elem_orig[i];

    const mt_int* s = surf.e2n_con.data() + surf.e2n_dsp[i];
    const mt_int* e = mesh.e2n_con.data() + mesh.e2n_dsp[mesh_eidx];

    MT_USET<mt_int> si;
    si.insert(s, s + surf.e2n_cnt[i]);

    mt_int v = -1;
    for(mt_int j = 0; j < mesh.e2n_cnt[mesh_eidx]; j++)
      if(si.count(e[j]) == 0) {
        v = e[j];
        break;
      }

    // now we can write the nbc data
    nbc.sp_vtx[i] = v;
    nbc.eidx[i]   = mesh_eidx;
    nbc.tag[i]    = mesh.etags[mesh_eidx];
  }
}

void unified_surface_from_tags(const mt_meshdata & mesh,
                               const mt_vector<MT_USET<mt_int> > & tags,
                               mt_meshdata & surfmesh,
                               MT_USET<mt_int> * vtx_set)
{
  MT_MAP<triple<mt_int>, tri_sele>     full_trimap;
  MT_MAP<quadruple<mt_int>, quad_sele> full_quadmap;

  for(size_t i=0; i<tags.size(); i++)
  {
    MT_MAP<triple<mt_int>, tri_sele>     cur_trimap;
    MT_MAP<quadruple<mt_int>, quad_sele> cur_quadmap;

    if(vtx_set != NULL) {
      mt_meshgraph mg;
      extract_tagged_meshgraph(tags[i], mesh, mg);
      vtx_set->insert(mg.e2n_con.begin(), mg.e2n_con.end());
      compute_surface(mg.etype, mg.e2n_cnt, mg.e2n_con, cur_trimap, cur_quadmap);
    }
    else {
      compute_surface(mesh, tags[i], cur_trimap, cur_quadmap);
    }

    surface_union(full_trimap, full_quadmap, cur_trimap, cur_quadmap);
  }
  mt_vector<mt_int> elem_orig;
  surfmap_to_vector(full_trimap, full_quadmap, surfmesh.e2n_con, elem_orig);
  surfmesh.e2n_cnt.assign(surfmesh.e2n_con.size() / 3, 3);
  surfmesh.etype.assign(surfmesh.e2n_cnt.size(), Tri);

  #if 0
  surfmesh.xyz = mesh.xyz;
  surfmesh.etags.assign(surfmesh.e2n_cnt.size(), mt_int(0));
  write_mesh_selected(surfmesh, "vtk", "debug_surfmesh");
  #endif
}

void unified_surface_from_list(const std::string & list,
                               const char delimiter,
                               mt_meshdata & surfmesh)
{
  struct triple<mt_int> surf;
  struct tri_sele       sele;
  MT_MAP<triple<mt_int>, tri_sele>     trimap;
  MT_MAP<quadruple<mt_int>, quad_sele> quadmap;
  MT_MAP<triple<mt_int>, tri_sele>::iterator it;

  mt_vector<std::string> surfs;
  split_string(list, delimiter, surfs);

  for(size_t i=0; i<surfs.size(); i++)
  {
    fixBasename(surfs[i]);
    std::string curr_surf_file = surfs[i] + SURF_EXT;
    std::cout << "Reading surface: " << curr_surf_file << std::endl;

    mt_meshdata curr_surf;
    readElements(curr_surf, curr_surf_file);
    size_t nele = curr_surf.e2n_cnt.size();

    // generate unique surface
    for(size_t eidx=0; eidx < nele; eidx++)
    {
      mt_int n1 = curr_surf.e2n_con[eidx*3+0];
      mt_int n2 = curr_surf.e2n_con[eidx*3+1];
      mt_int n3 = curr_surf.e2n_con[eidx*3+2];

      sele.v1 = n1, sele.v2 = n2, sele.v3 = n3; sele.eidx = eidx;
      sortTriple(n1, n2, n3, surf.v1, surf.v2, surf.v3);

      it = trimap.find(surf);
      if(it == trimap.end()) trimap[surf] = sele;
    }
  }
  mt_vector<mt_int> elem_orig;
  surfmap_to_vector(trimap, quadmap, surfmesh.e2n_con, elem_orig);
  surfmesh.e2n_cnt.assign(surfmesh.e2n_con.size() / 3, 3);
  surfmesh.etype.assign(surfmesh.e2n_cnt.size(), Tri);
  surfmesh.etags.assign(surfmesh.e2n_cnt.size(), mt_int(0));
}

void mesh_resize_elemdata(mt_meshdata & mesh, size_t size_elem, size_t size_con)
{
  bool twoFib = mesh.lon.size() == (mesh.e2n_cnt.size() * 6);

  mesh.e2n_cnt.resize(size_elem);
  mesh.e2n_con.resize(size_con);

  mesh.etype.resize(size_elem);
  mesh.etags.resize(size_elem);
  mesh.lon  .resize(twoFib ? size_elem * 6 : size_elem * 3);
}

void mesh_add_elem(mt_meshdata & mesh, const elem_t type, const mt_int* con, const mt_int tag)
{
  short numfibvals = mesh.e2n_cnt.size() * 6 == mesh.lon.size() ? 6 : 3;
  mt_int numcon = -1;

  switch(type) {
    case Line:
      numcon = 2;
      break;

    case Tri:
      numcon = 3;
      break;

    case Quad:
    case Tetra:
      numcon = 4;
      break;

    case Pyramid:
      numcon = 5;
      break;

    case Prism:
      numcon = 6;
      break;

    case Hexa:
      numcon = 8;
      break;
  }

  mesh.e2n_cnt.push_back(numcon);
  mesh.etype.push_back(type);
  mesh.etags.push_back(tag);

  for(mt_int i=0; i<numcon; i++)
    mesh.e2n_con.push_back(con[i]);

  for(short i=0; i<numfibvals; i++)
    mesh.lon.push_back(0.0);
}

void surface_to_mesh(const mt_meshdata & mesh,
                     const mt_vector<mt_int> & surf_con,
                     const nbc_data & nbc,
                     mt_meshdata & surfmesh)
{
  bool lon2 = mesh.lon.size() == mesh.e2n_cnt.size()*6;
  size_t nelem = surf_con.size() / 3;

  // make sure we have proper nbc data
  assert(nelem = nbc.eidx.size());

  // all elements are triangles
  surfmesh.e2n_cnt.assign(nelem, 3);
  surfmesh.etype.assign(nelem, Tri);
  surfmesh.etags.assign(nbc.tag.begin(), nbc.tag.end());

#if 0
  surfmesh.e2n_con.assign(surf_con.begin(), surf_con.end());
#else
  surfmesh.e2n_con.resize(surf_con.size());
  // flip surface connectivity to have outward facing normals
  for(size_t i=0; i<nelem; i++) {
    surfmesh.e2n_con[i*3+0] = surf_con[i*3+2];
    surfmesh.e2n_con[i*3+1] = surf_con[i*3+1];
    surfmesh.e2n_con[i*3+2] = surf_con[i*3+0];
  }
#endif

  mt_vector<mt_int> nod(surf_con);
  binary_sort(nod); unique_resize(nod);
  nod.reallocate();

  // remap surface indices
  map_glob2loc(nod, surfmesh.e2n_con);

  // copy vertex coords
  size_t nnodes = nod.size();
  surfmesh.xyz.resize(nnodes*3);
  for(size_t i=0; i<nnodes; i++) {
    mt_int n = nod[i];
    surfmesh.xyz[i*3+0] = mesh.xyz[n*3+0];
    surfmesh.xyz[i*3+1] = mesh.xyz[n*3+1];
    surfmesh.xyz[i*3+2] = mesh.xyz[n*3+2];
  }

  // copy fibers
  surfmesh.lon.resize( lon2 ? nelem*6 : nelem*3 );
  if(lon2) {
    for(size_t i=0; i<nelem; i++) {
      mt_int eidx = nbc.eidx[i];
      surfmesh.lon[i*6+0] = mesh.lon[eidx*6+0];
      surfmesh.lon[i*6+1] = mesh.lon[eidx*6+1];
      surfmesh.lon[i*6+2] = mesh.lon[eidx*6+2];
      surfmesh.lon[i*6+3] = mesh.lon[eidx*6+3];
      surfmesh.lon[i*6+4] = mesh.lon[eidx*6+4];
      surfmesh.lon[i*6+5] = mesh.lon[eidx*6+5];
    }
  }
  else {
    for(size_t i=0; i<nelem; i++) {
      mt_int eidx = nbc.eidx[i];
      surfmesh.lon[i*3+0] = mesh.lon[eidx*3+0];
      surfmesh.lon[i*3+1] = mesh.lon[eidx*3+1];
      surfmesh.lon[i*3+2] = mesh.lon[eidx*3+2];
    }
  }
}

void nbc_data_into_surface(const mt_meshdata & mesh, const nbc_data & nbc,
                           mt_meshdata & surfmesh)
{
  bool lon2 = mesh.lon.size() == mesh.e2n_cnt.size()*6;
  size_t nelem = surfmesh.e2n_cnt.size();

  // make sure we have proper nbc data
  assert(nelem == nbc.eidx.size());
  surfmesh.etags.assign(nbc.tag.begin(), nbc.tag.end());

#if 1
  // flip surface connectivity to have outward facing normals in meshalyzer
  for(size_t i=0, ridx=0, widx=0; i<nelem; i++) {
    switch(surfmesh.e2n_cnt[i]) {
      case 3:
      {
        mt_int v1 = surfmesh.e2n_con[ridx++];
        mt_int v2 = surfmesh.e2n_con[ridx++];
        mt_int v3 = surfmesh.e2n_con[ridx++];
        surfmesh.e2n_con[widx++] = v3;
        surfmesh.e2n_con[widx++] = v2;
        surfmesh.e2n_con[widx++] = v1;
        break;
      }

      case 4: {
        mt_int v1 = surfmesh.e2n_con[ridx++];
        mt_int v2 = surfmesh.e2n_con[ridx++];
        mt_int v3 = surfmesh.e2n_con[ridx++];
        mt_int v4 = surfmesh.e2n_con[ridx++];
        surfmesh.e2n_con[widx++] = v4;
        surfmesh.e2n_con[widx++] = v3;
        surfmesh.e2n_con[widx++] = v2;
        surfmesh.e2n_con[widx++] = v1;
        break;
      }
    }
  }
#endif

  mt_vector<mt_int> nod(surfmesh.e2n_con);
  binary_sort(nod); unique_resize(nod);
  nod.reallocate();

  // remap surface indices
  map_glob2loc(nod, surfmesh.e2n_con);

  // copy vertex coords
  size_t nnodes = nod.size();
  surfmesh.xyz.resize(nnodes*3);
  for(size_t i=0; i<nnodes; i++) {
    mt_int n = nod[i];
    surfmesh.xyz[i*3+0] = mesh.xyz[n*3+0];
    surfmesh.xyz[i*3+1] = mesh.xyz[n*3+1];
    surfmesh.xyz[i*3+2] = mesh.xyz[n*3+2];
  }

  // copy fibers
  surfmesh.lon.resize( lon2 ? nelem*6 : nelem*3 );
  if(lon2) {
    for(size_t i=0; i<nelem; i++) {
      mt_int eidx = nbc.eidx[i];
      surfmesh.lon[i*6+0] = mesh.lon[eidx*6+0];
      surfmesh.lon[i*6+1] = mesh.lon[eidx*6+1];
      surfmesh.lon[i*6+2] = mesh.lon[eidx*6+2];
      surfmesh.lon[i*6+3] = mesh.lon[eidx*6+3];
      surfmesh.lon[i*6+4] = mesh.lon[eidx*6+4];
      surfmesh.lon[i*6+5] = mesh.lon[eidx*6+5];
    }
  }
  else {
    for(size_t i=0; i<nelem; i++) {
      mt_int eidx = nbc.eidx[i];
      surfmesh.lon[i*3+0] = mesh.lon[eidx*3+0];
      surfmesh.lon[i*3+1] = mesh.lon[eidx*3+1];
      surfmesh.lon[i*3+2] = mesh.lon[eidx*3+2];
    }
  }
}

void remove_nodes_from_surf(const MT_USET<mt_int> & nodes,
                            mt_vector<mt_int> & surf_cnt,
                            mt_vector<mt_int> & surf_con,
                            struct nbc_data & nbc)
{
  size_t nele = surf_cnt.size();
  size_t widx=0, widx_con=0;
  size_t ridx_con=0;

  for(size_t i=0; i<nele; i++) {
    switch(surf_cnt[i]) {
      case 3:
      {
        mt_int vtx0 = surf_con[ridx_con++];
        mt_int vtx1 = surf_con[ridx_con++];
        mt_int vtx2 = surf_con[ridx_con++];

        bool ok0 = (nodes.count(vtx0) == 0);
        bool ok1 = (nodes.count(vtx1) == 0);
        bool ok2 = (nodes.count(vtx2) == 0);

        if(ok0 && ok1 && ok2) {
          surf_con[widx_con++] = vtx0;
          surf_con[widx_con++] = vtx1;
          surf_con[widx_con++] = vtx2;

          surf_cnt[widx]   = surf_cnt[i];

          nbc.eidx[widx]   = nbc.eidx[i];
          nbc.sp_vtx[widx] = nbc.sp_vtx[i];
          nbc.tag[widx]    = nbc.tag[i];

          widx++;
        }
        break;
      }

      case 4:
      {
        mt_int vtx0 = surf_con[ridx_con++];
        mt_int vtx1 = surf_con[ridx_con++];
        mt_int vtx2 = surf_con[ridx_con++];
        mt_int vtx3 = surf_con[ridx_con++];

        bool ok0 = (nodes.count(vtx0) == 0);
        bool ok1 = (nodes.count(vtx1) == 0);
        bool ok2 = (nodes.count(vtx2) == 0);
        bool ok3 = (nodes.count(vtx3) == 0);

        if(ok0 && ok1 && ok2 && ok3) {
          surf_con[widx_con++] = vtx0;
          surf_con[widx_con++] = vtx1;
          surf_con[widx_con++] = vtx2;
          surf_con[widx_con++] = vtx3;

          surf_cnt[widx]   = surf_cnt[i];

          nbc.eidx[widx]   = nbc.eidx[i];
          nbc.sp_vtx[widx] = nbc.sp_vtx[i];
          nbc.tag[widx]    = nbc.tag[i];

          widx++;
        }
        break;
      }

      default: break;
    }
  }

  surf_con.resize(widx_con);
  surf_cnt.resize(widx);

  nbc.eidx.resize(widx);
  nbc.sp_vtx.resize(widx);
  nbc.tag.resize(widx);
}


void remove_nodes_from_surf(const MT_USET<mt_int> & nodes,
                            mt_vector<mt_int> & surf_cnt,
                            mt_vector<mt_int> & surf_con)
{
  size_t nele = surf_cnt.size();
  size_t widx=0, widx_con=0;
  size_t ridx_con=0;

  for(size_t i=0; i<nele; i++) {
    switch(surf_cnt[i]) {
      case 3:
      {
        mt_int vtx0 = surf_con[ridx_con++];
        mt_int vtx1 = surf_con[ridx_con++];
        mt_int vtx2 = surf_con[ridx_con++];

        bool ok0 = (nodes.count(vtx0) == 0);
        bool ok1 = (nodes.count(vtx1) == 0);
        bool ok2 = (nodes.count(vtx2) == 0);

        if(ok0 && ok1 && ok2) {
          surf_con[widx_con++] = vtx0;
          surf_con[widx_con++] = vtx1;
          surf_con[widx_con++] = vtx2;

          surf_cnt[widx]   = surf_cnt[i];

          widx++;
        }
        break;
      }

      case 4:
      {
        mt_int vtx0 = surf_con[ridx_con++];
        mt_int vtx1 = surf_con[ridx_con++];
        mt_int vtx2 = surf_con[ridx_con++];
        mt_int vtx3 = surf_con[ridx_con++];

        bool ok0 = (nodes.count(vtx0) == 0);
        bool ok1 = (nodes.count(vtx1) == 0);
        bool ok2 = (nodes.count(vtx2) == 0);
        bool ok3 = (nodes.count(vtx3) == 0);

        if(ok0 && ok1 && ok2 && ok3) {
          surf_con[widx_con++] = vtx0;
          surf_con[widx_con++] = vtx1;
          surf_con[widx_con++] = vtx2;
          surf_con[widx_con++] = vtx3;

          surf_cnt[widx]   = surf_cnt[i];

          widx++;
        }
        break;
      }

      default: break;
    }
  }

  surf_con.resize(widx_con);
  surf_cnt.resize(widx);
}

void remove_elems_from_surf(const mt_vector<bool> & keep,
                            mt_vector<mt_int> & surf_cnt,
                            mt_vector<mt_int> & surf_con,
                            struct nbc_data & nbc)
{
  size_t nele = surf_cnt.size();
  size_t widx=0, widx_con=0;
  size_t ridx_con=0;

  for(size_t i=0; i<nele; i++) {
    switch(surf_cnt[i]) {
      case 3:
      {
        mt_int vtx0 = surf_con[ridx_con++];
        mt_int vtx1 = surf_con[ridx_con++];
        mt_int vtx2 = surf_con[ridx_con++];

        if(keep[i]) {
          surf_con[widx_con++] = vtx0;
          surf_con[widx_con++] = vtx1;
          surf_con[widx_con++] = vtx2;

          surf_cnt[widx]   = surf_cnt[i];

          nbc.eidx[widx]   = nbc.eidx[i];
          nbc.sp_vtx[widx] = nbc.sp_vtx[i];
          nbc.tag[widx]    = nbc.tag[i];

          widx++;
        }
        break;
      }

      case 4:
      {
        mt_int vtx0 = surf_con[ridx_con++];
        mt_int vtx1 = surf_con[ridx_con++];
        mt_int vtx2 = surf_con[ridx_con++];
        mt_int vtx3 = surf_con[ridx_con++];

        if(keep[i]) {
          surf_con[widx_con++] = vtx0;
          surf_con[widx_con++] = vtx1;
          surf_con[widx_con++] = vtx2;
          surf_con[widx_con++] = vtx3;

          surf_cnt[widx]   = surf_cnt[i];

          nbc.eidx[widx]   = nbc.eidx[i];
          nbc.sp_vtx[widx] = nbc.sp_vtx[i];
          nbc.tag[widx]    = nbc.tag[i];

          widx++;
        }
        break;
      }

      default: break;
    }
  }

  surf_con.resize(widx_con);
  surf_cnt.resize(widx);

  nbc.eidx.resize(widx);
  nbc.sp_vtx.resize(widx);
  nbc.tag.resize(widx);
}




void compute_element_surface_normals(const mt_meshdata & surfmesh,
                                     const mt_vector<mt_real> & xyz,
                                     mt_vector<mt_real> & snrml)
{
  assert(surfmesh.e2n_dsp.size() > 0);

  size_t nele_surf = surfmesh.e2n_cnt.size();
  snrml.resize(nele_surf*3);

  for(size_t eidx = 0; eidx < nele_surf; eidx++)
  {
    mt_int dsp = surfmesh.e2n_dsp[eidx];

    mt_int v1 = surfmesh.e2n_con[dsp+0];
    mt_int v2 = surfmesh.e2n_con[dsp+1];
    mt_int v3 = surfmesh.e2n_con[dsp+2];

    // also for quads, we compute a triangle surface normal
    vec3r n = triangle_normal(v1, v2, v3, xyz);

    snrml[eidx*3+0] = n.x;
    snrml[eidx*3+1] = n.y;
    snrml[eidx*3+2] = n.z;
  }
}

void compute_nodal_surface_normals(const mt_meshdata & surfmesh,
                                   const mt_vector<mt_real> & xyz,
                                   mt_vector<mt_real> & snrml)
{
  // full mesh connectivity needs to be set up
  check_nonzero(surfmesh.e2n_dsp.size(), __func__);
  check_nonzero(surfmesh.n2e_cnt.size(), __func__);
  check_nonzero(surfmesh.n2n_cnt.size(), __func__);

  size_t nele = surfmesh.e2n_cnt.size(), nnod = surfmesh.n2n_cnt.size();
  mt_vector<vec3r> enrmls(nele), nnrmls(nnod);

  for(size_t eidx = 0; eidx < nele; eidx++)
  {
    mt_int dsp = surfmesh.e2n_dsp[eidx];

    mt_int v1 = surfmesh.e2n_con[dsp+0];
    mt_int v2 = surfmesh.e2n_con[dsp+1];
    mt_int v3 = surfmesh.e2n_con[dsp+2];

    // also for quads, we compute a triangle surface normal
    enrmls[eidx] = triangle_normal(v1, v2, v3, xyz);
  }

#if 1
  //Angle averaged vertex normals
  // according to "A Comparison of Algorithms for Vertex Normal Computation" this yields
  // smoother results than the default are weights
  #ifdef OPENMP
  #pragma omp parallel for schedule(guided)
  #endif
  for(size_t nidx = 0; nidx < nnod; nidx++)
  {
    if(surfmesh.n2e_cnt[nidx]) {
      vec3r avrg(0,0,0), temp;
      mt_vector<vec3r> edges(3);
      mt_int start = surfmesh.n2e_dsp[nidx], stop = start + surfmesh.n2e_cnt[nidx];
      for(mt_int i=start; i < stop; i++)
      {
        mt_int eidx = surfmesh.n2e_con[i];
        const mt_int* con = surfmesh.e2n_con.data() + surfmesh.e2n_dsp[eidx];
        vec3r refp(xyz.data() + 3 * nidx);

        // the algorithm picks the edges to the first two vertices that are not nidx of
        // each connected facet for computing the coefficients
        short widx = 0;
        for(short k=0; k < 3; k++) {
          if(con[k] != mt_int(nidx)) {
            temp.get(xyz.data() + 3 * con[k]);
            edges[widx++] = temp;
          }
        }
        mt_real alpha = angle_between_line(refp, edges[0], edges[1]);
        avrg += enrmls[eidx] * alpha;
      }
      nnrmls[nidx] = avrg;
    }
  }
#else
  elemData_to_nodeData(surfmesh, enrmls, nnrmls);
#endif

  size_t num_zero = 0, num_fixed = 0;

  snrml.assign(nnod * 3, 0.0);
  for(size_t nidx = 0; nidx < nnod; nidx++)
  {
    if(surfmesh.n2e_cnt[nidx] > 0) {
      if(nnrmls[nidx].length2() == 0.0)
      {
        num_zero++;

        vec3r avrg(0,0,0);
        mt_int start = surfmesh.n2e_dsp[nidx], stop = start + surfmesh.n2e_cnt[nidx];
        for(mt_int i=start; i < stop; i++)
        {
          mt_int eidx = surfmesh.n2e_con[i];
          avrg += enrmls[eidx];
        }

        if(avrg.length2() > 0) {
          // we take the averaged normal
          nnrmls[nidx] = avrg;
          num_fixed++;
        } else {
          // we take the normal of the first element
          nnrmls[nidx] = enrmls[surfmesh.n2e_con[start]];
        }
      }

      nnrmls[nidx].normalize();
      snrml[nidx*3+0] = nnrmls[nidx].x;
      snrml[nidx*3+1] = nnrmls[nidx].y;
      snrml[nidx*3+2] = nnrmls[nidx].z;
    }
  }

  if(num_zero) {
    fprintf(stderr, "%s warning: %ld zero normals were computed. %ld were fixed by using homogenous weights.\n"
                    "Maybe the surface element orientation is not consistent.\n", __func__,
                    num_zero, num_fixed);
  }
}

bool apply_normal_orientation(mt_meshdata & surf, const bool inward)
{
  if(surf.e2n_dsp.size() == 0)
    bucket_sort_offset(surf.e2n_cnt, surf.e2n_dsp);

  float edge_len_est = min_edgelength_estimate(surf, true, true);
  kdtree surf_tree(10);
  surf_tree.build_tree(surf);
  bool have_flipped = false;

  for(size_t i = 0; i < surf.e2n_cnt.size(); i++) {
    mt_int* con = surf.e2n_con.data() + surf.e2n_dsp[i];
    mt_int v0 = con[0], v1 = con[1], v2 = con[2];

    vec3f n = triangle_normal(v0, v1, v2, surf.xyz);
    vec3f c = triangle_centerpoint(v0, v1, v2, surf.xyz);
    vec3f testpoint = c + (n * (0.01f * edge_len_est));

    bool is_inside = inside_closed_surface(surf_tree, testpoint);
    bool needs_flip = is_inside != inward;

    // if the current orientation does not match the desired one, we flip the normals
    if(needs_flip) {
      con[0] = v2;
      con[1] = v1;
      con[2] = v0;
      have_flipped = true;
    }
  }

  return have_flipped;
}

bool apply_normal_orientation(mt_meshdata & surf, kdtree & surf_tree, mt_vector<mt_int> & eidx, const bool inward)
{
  if(surf.e2n_dsp.size() == 0)
    bucket_sort_offset(surf.e2n_cnt, surf.e2n_dsp);

  float edge_len_est = min_edgelength_estimate(surf, true, true);
  bool have_flipped = false;

  for(mt_int i : eidx) {
    mt_int* con = surf.e2n_con.data() + surf.e2n_dsp[i];
    mt_int v0 = con[0], v1 = con[1], v2 = con[2];

    vec3f n = triangle_normal(v0, v1, v2, surf.xyz);
    vec3f c = triangle_centerpoint(v0, v1, v2, surf.xyz);
    vec3f testpoint = c + (n * (0.01f * edge_len_est));

    bool is_inside = inside_closed_surface(surf_tree, testpoint);
    bool needs_flip = is_inside != inward;

    // if the current orientation does not match the desired one, we flip the normals
    if(needs_flip) {
      con[0] = v2;
      con[1] = v1;
      con[2] = v0;
      have_flipped = true;
    }
  }

  return have_flipped;
}

void remove_bad_surface_edges(const mt_meshdata & mesh,
                              mt_meshdata & surfmesh)
{
  // full mesh connectivity needs to be set up
  assert(surfmesh.n2n_cnt.size() > 0 && surfmesh.n2e_cnt.size() > 0);

  mt_vector<mt_real> snrml;
  compute_nodal_surface_normals(surfmesh, mesh.xyz, snrml);

  mt_vector<mt_int> & cnt = surfmesh.n2n_cnt;
  mt_vector<mt_int> & con = surfmesh.n2n_con;

  size_t ridx=0, widx=0;

  for(size_t i=0; i<cnt.size(); i++)
  {
    mt_int c = 0;  // current count
    mt_point<mt_real> cn(snrml.data() + i*3);  // current normal

    for(mt_int j=0; j<cnt[i]; j++)
    {
      mt_point<mt_real> tn(snrml.data() + con[ridx]*3);

      if( cn.scaProd(tn) > mt_real(0.0) )
      {
        c++;
        con[widx++] = con[ridx];
      }
      ridx++;
    }
    cnt[i] = c;
  }
  con.resize(widx);
  bucket_sort_offset(surfmesh.n2n_cnt, surfmesh.n2n_dsp);
}


void psdata_to_mesh(const mt_psdata & ps, mt_meshdata & psmesh)
{
  struct triple_double_hash {
    static inline bool cmp(triple<double> a, triple<double> b)
    {
      return int(a.v1*1000) == int(b.v1*1000) &&
             int(a.v2*1000) == int(b.v2*1000) &&
             int(a.v3*1000) == int(b.v3*1000);
    }

    static inline unsigned int hash(triple<double> a)
    {
      unsigned int h = hashmap::mkhash_init;
      h = hashmap::mkhash(h, hashmap::hash_ops<int32_t>::hash((int)a.v1*1000));
      h = hashmap::mkhash(h, hashmap::hash_ops<int32_t>::hash((int)a.v2*1000));
      h = hashmap::mkhash(h, hashmap::hash_ops<int32_t>::hash((int)a.v3*1000));
      return h;
    }
  };

  hashmap::unordered_map<triple<double>, mt_int, triple_double_hash> vtx_map;

  size_t num_elem = 0;

  // first generate vertex indexing
  mt_int gidx = 0;
  for(size_t cidx = 0; cidx < ps.cables.size(); cidx++)
  {
    size_t nvtx = ps.cables[cidx].pts.size() / 3;
    num_elem += nvtx - 1;

    for(size_t i=0; i<nvtx; i++)
    {
      triple<double> pt = {ps.cables[cidx].pts[i*3+0],
                           ps.cables[cidx].pts[i*3+1],
                           ps.cables[cidx].pts[i*3+2]};
      auto it = vtx_map.find(pt);
      if(it == vtx_map.end())
        vtx_map[pt] = gidx++;
    }
  }

  // generate elem data
  psmesh.etype.assign(num_elem, Line);
  psmesh.etags.resize(num_elem);
  psmesh.lon.assign(num_elem * 3, mt_real(0));
  psmesh.e2n_cnt.assign(num_elem, 2);
  psmesh.e2n_con.resize(num_elem*2);

  for(size_t cidx = 0, eidx = 0; cidx < ps.cables.size(); cidx++)
  {
    size_t nvtx = ps.cables[cidx].pts.size() / 3;
    triple<double> pt;
    for(size_t i=1; i<nvtx; i++, eidx++)
    {
      const mt_real* pts = ps.cables[cidx].pts.data() + (i-1)*3;

      pt = {pts[0], pts[1], pts[2]};
      mt_int v1 = vtx_map[pt];
      pt = {pts[3], pts[4], pts[5]};
      mt_int v2 = vtx_map[pt];

      // connectivity
      psmesh.e2n_con[eidx*2+0] = v1;
      psmesh.e2n_con[eidx*2+1] = v2;
      // tag
      psmesh.etags[eidx] = cidx;
      // fiber
      psmesh.lon[eidx*3] = 1;
    }
  }

  // generate vertex data
  psmesh.xyz.resize(vtx_map.size()*3);
  for(auto it = vtx_map.begin(); it != vtx_map.end(); ++it)
  {
    psmesh.xyz[it->second*3+0] = it->first.v1;
    psmesh.xyz[it->second*3+1] = it->first.v2;
    psmesh.xyz[it->second*3+2] = it->first.v3;
  }
}

void generate_coord_map(mt_meshdata & mesh,
                        MT_MAP<triple<mt_int>,mt_int> & cmap,
                        int & scale_factor)
{
  // we have to take minimal edge lenght here in order to compute the correct scale factor
  // since the average can be enough oders of magnitude higher than the minimum for us to
  // loose coordinates
  // TODO: Elias
  mt_real edge_len = min_edgelength_estimate(mesh, true);

  // we treat edge_len == 0 with an error
  if(edge_len == 0.0) {
    fprintf(stderr, "%s WARNING: Zero minimum edge length detected!\nRetry edge length estimate with nonzero minimum!\n", __func__);
    // We can take the non-zero minimum, aka the first edge that has non zero length and base the estimate on that
    edge_len = min_edgelength_estimate(mesh, true, true);
  }

  int dp = get_dec_power_estimate(edge_len);
  int scale_power = 3;

  if(dp < 0)
    scale_power += -dp;

  scale_factor = pow(10.0, scale_power);

  for(size_t i=0; i<mesh.xyz.size() / 3; i++) {
    triple<mt_int> t;
    t.v1 = mt_int(mesh.xyz[i*3+0]*scale_factor);
    t.v2 = mt_int(mesh.xyz[i*3+1]*scale_factor);
    t.v3 = mt_int(mesh.xyz[i*3+2]*scale_factor);
    cmap[t] = i;
  }
}

void mesh_union(mt_meshdata & mesh1, const mt_meshdata & mesh2, const bool error_on_empty_intf)
{
  // we want to make sure that the dsp is always up to date
  bucket_sort_offset(mesh1.e2n_cnt, mesh1.e2n_dsp);

  // We are identifying equal vertex coordinates in both meshes by converting them
  // to integers and putting them into a hashmap.
  int scale_factor;
  MT_MAP<triple<mt_int>,mt_int> mesh1_coords;
  generate_coord_map(mesh1, mesh1_coords, scale_factor);

  MT_MAP<mt_int,mt_int> mesh2_to_mesh1_interf;
  for(size_t i=0; i<mesh2.xyz.size() / 3; i++) {
    triple<mt_int> t;
    t.v1 = mt_int(mesh2.xyz[i*3+0]*scale_factor);
    t.v2 = mt_int(mesh2.xyz[i*3+1]*scale_factor);
    t.v3 = mt_int(mesh2.xyz[i*3+2]*scale_factor);

    auto it = mesh1_coords.find(t);

    if(it != mesh1_coords.end()) {
      mt_int mesh2_idx = i, mesh1_idx = it->second;

      vec3r m1_pt(mesh1.xyz.data()+mesh1_idx*3), m2_pt(mesh2.xyz.data()+mesh2_idx*3);
      mt_real len = (m1_pt-m2_pt).length();

      if(len)
        fprintf(stderr, "%s warning: unifying vertices of distance %f\n", __func__, len);

      mesh2_to_mesh1_interf[mesh2_idx] = mesh1_idx;
    }
  }

  std::cout << "Interface size: " << mesh2_to_mesh1_interf.size() << std::endl;
  if(error_on_empty_intf && mesh2_to_mesh1_interf.size() == 0) {
    fprintf(stderr, "%s error: Empty mesh interface is not allowed! Aborting!\n", __func__);
    return;
  }

  // construct whole mapping vector for mesh2
  size_t nnod_msh1 = mesh1.xyz.size() / 3, nnod_msh2 = mesh2.xyz.size() / 3;
  mt_vector<mt_int> mesh2_to_wholemesh(nnod_msh2);
  for(size_t n=0; n<nnod_msh2; n++) {
    if(mesh2_to_mesh1_interf.count(n)) // interface node
      mesh2_to_wholemesh[n] = mesh2_to_mesh1_interf[n];
    else  // other nodes
      mesh2_to_wholemesh[n] = nnod_msh1 + n;
  }

  mt_vector<mt_int> mesh2_e2n_con = mesh2.e2n_con;
  for(size_t i=0; i<mesh2_e2n_con.size(); i++)
    mesh2_e2n_con[i] = mesh2_to_wholemesh[mesh2_e2n_con[i]];

  mesh1.e2n_cnt.append(mesh2.e2n_cnt.begin(), mesh2.e2n_cnt.end());
  mesh1.e2n_con.append(mesh2_e2n_con.begin(), mesh2_e2n_con.end());

  mesh1.etype.append(mesh2.etype.begin(), mesh2.etype.end());
  mesh1.etags.append(mesh2.etags.begin(), mesh2.etags.end());
  mesh1.lon.append(mesh2.lon.begin(), mesh2.lon.end());

  mesh1.xyz.append(mesh2.xyz.begin(), mesh2.xyz.end());

  // reindex to get rid of duplicated vertices
  reindex_nodes(mesh1, true);

  bucket_sort_offset(mesh1.e2n_cnt, mesh1.e2n_dsp);
  correct_duplicate_elements(mesh1);
}


/// function evaluating if nodes are inside a closed surface or not
void nodes_in_surface(const mt_meshdata & mesh,
                      const kdtree & tree,
                      mt_vector<bool> & in_surf)
{
  mt_vector<vec3f> points;
  array_to_points(mesh.xyz, points);
  in_surf.resize(points.size());

  #ifdef OPENMP
  #pragma omp parallel for schedule(guided)
  #endif
  for(size_t i=0; i<points.size(); i++) {
    if(tree.in_root_bbox(points[i]))
      in_surf[i] = inside_closed_surface(tree, points[i]);
    else
      in_surf[i] = false;
  }
}


void sample_elem_tags(mt_meshdata & mesh, mt_mask & tag_found, const kdtree & tree, const mt_int newtag,
                     const int sampling_type)
{
  if(mesh.e2n_dsp.size() == 0)
    bucket_sort_offset(mesh.e2n_cnt, mesh.e2n_dsp);

  switch(sampling_type)
  {
    default:
    case 0:
    {
      // this version samples elements with all nodes inside the surface
      mt_vector<bool> nod_in_surf;
      nodes_in_surface(mesh, tree, nod_in_surf);
      for(size_t eidx=0; eidx < mesh.e2n_cnt.size(); eidx++)
      {
        if(tag_found.count(eidx) == 0) {
          mt_int start = mesh.e2n_dsp[eidx], stop = start + mesh.e2n_cnt[eidx];
          mt_int ncon = 0;
          for(mt_int j=start; j<stop; j++) {
            mt_int nidx = mesh.e2n_con[j];
            if(nod_in_surf[nidx]) ncon++;
          }

          if(ncon == mesh.e2n_cnt[eidx]) {
            mesh.etags[eidx] = newtag;
            tag_found.insert(eidx);
          }
        }
      }
      break;
    }
    case 1:
    {
      // this version samples elements with their center-point inside the surface

      #ifdef OPENMP
      #pragma omp parallel for schedule(guided)
      #endif
      for(size_t eidx=0; eidx < mesh.e2n_cnt.size(); eidx++)
      {
        if(tag_found.count(eidx) == 0) {
          vec3r ctr(0,0,0);
          mt_int start = mesh.e2n_dsp[eidx], stop = start + mesh.e2n_cnt[eidx];
          for(mt_int j=start; j<stop; j++)
          {
            mt_int nidx = mesh.e2n_con[j];
            ctr += vec3r(mesh.xyz.data() + nidx*3);
          }
          ctr /= mt_real(mesh.e2n_cnt[eidx]);

          if(inside_closed_surface(tree, ctr)) {
            mesh.etags[eidx] = newtag;
            tag_found.insert(eidx);
          }
        }
      }
      break;
    }
  }
}

void tet_gradient(const mt_int* con, const mt_vector<mt_real> & xyz,
                  const mt_vector<mt_real> & data, mt_real* grad)
{
  mt_real rhs[4];
  dmat<mt_real> M(4,4);

  for (int i=0;i<4;i++) {
    mt_int v = con[i];
    M[i][0] = 1.;
    M[i][1] = xyz[v*3+0];
    M[i][2] = xyz[v*3+1];
    M[i][3] = xyz[v*3+2];

    rhs[i]  = data[v];
  }

  M.lu_decomp();
  M.lu_solve(rhs);

  grad[0] = rhs[1];
  grad[1] = rhs[2];
  grad[2] = rhs[3];
}

void compute_gradient(const mt_meshdata & mesh, const mt_vector<mt_real> & data,
                      const bool nodal_input, const bool nodal_output,
                      mt_vector<vec3r> & grad, mt_vector<mt_real> & mag)
{
  mt_vector<mt_real> emag, ndata, *_mag = NULL;
  const mt_vector<mt_real> *_data = NULL;
  mt_vector<vec3r>   egrad, *_grad = NULL;
  size_t numele = mesh.e2n_cnt.size();

  // the core algorithm computes gradients of a nodal function into the
  // element centers. is we want different input / output representation
  // we need to interpolate
  if(nodal_input) {
    _data = &data;
  }
  else {
    elemData_to_nodeData(mesh, data, ndata);
    _data = &ndata;
  }

  if(nodal_output) {
    _grad = &egrad;
    _mag  = &emag;
  }
  else {
    _grad = &grad;
    _mag  = &mag;
  }

  // references to the elment based data
  const mt_vector<mt_real> & d = *_data;
  mt_vector<mt_real> & m = *_mag;
  mt_vector<vec3r>   & g = *_grad;

  m.resize(numele), g.resize(numele);

  #ifdef OPENMP
  #pragma omp parallel for schedule(guided)
  #endif
  for(size_t eidx=0; eidx<mesh.e2n_cnt.size(); eidx++) {
    mt_real grad_vec[3];
    switch(mesh.etype[eidx]) {
      case Tetra:
      {
        const mt_int* con = mesh.e2n_con.data() + mesh.e2n_dsp[eidx];
        tet_gradient(con, mesh.xyz, d, grad_vec);
        g[eidx] = vec3r(grad_vec);
        m[eidx] = g[eidx].length();
      }
      default: break;
    }
  }

  if(nodal_output) {
    elemData_to_nodeData(mesh, egrad, grad);
    elemData_to_nodeData(mesh, emag,  mag);
  }
}

void read_surf_info(const mt_meshdata & mesh, mt_meshdata & surf, mt_vector<mt_int> & eidx,
                    std::string basefile)
{
  fixBasename(basefile);
  std::string surffile = basefile + SURF_EXT;
  std::string vtxfile  = basefile + SURF_EXT + VTX_EXT;
  std::string nbcfile  = basefile + NBC_EXT;
  nbc_data nbc;

  if(file_exists(nbcfile)) {
    std::cout << "Reading nbc file: " << nbcfile << std::endl;
    read_nbc(nbc, surf.e2n_con, nbcfile);
    surf.e2n_cnt.assign(surf.e2n_con.size() / 3, mt_int(3));
    surf.etype  .assign(surf.e2n_con.size() / 3, Tri);
    eidx = nbc.eidx;
  }
  else if(file_exists(surffile)) {
    check_nonzero(mesh.n2e_cnt.size(), __func__);

    std::cout << "Reading surface: " << surffile << std::endl;
    readElements(surf, surffile);
    bucket_sort_offset(surf.e2n_cnt, surf.e2n_dsp);

    size_t nelems = surf.e2n_cnt.size();
    eidx.resize(nelems);

    std::cout << "Recovering volumetric element indices .. " << std::endl;

    #ifdef OPENMP
    #pragma omp parallel for schedule(dynamic, 10)
    #endif
    for(size_t e = 0; e < nelems; e++) {
      std::set<mt_int> elemset;
      mt_int * con = surf.e2n_con.data() + surf.e2n_dsp[e];

      elements_with_face(mesh, con[0], con[1], con[2], elemset);

      switch(elemset.size()) {
        default: fprintf(stderr, "%s error: More than 2 elements seem to share a face! Aborting!\n", __func__);
                 exit(1);

        case 0: eidx[e] = 0; break;

        case 1: eidx[e] = *elemset.begin(); break;

        case 2: {
          auto it = elemset.begin();
          mt_int e1 = *it; ++it; mt_int e2 = *it;
          vec3r ctr1   = barycenter(mesh, e1);
          vec3r ctr_tr = triangle_centerpoint(con[0], con[1], con[2], mesh.xyz);
          vec3r n      = triangle_normal(con[0], con[1], con[2], mesh.xyz);

          if(n.scaProd(unit_vector(ctr1 - ctr_tr)) < 0.0)
            eidx[e] = e1;
          else
            eidx[e] = e2;

          break;
        }
      }
    }
  }
  else if(file_exists(vtxfile)) {
    std::cout << "Reading vertex file: " << vtxfile << std::endl;
    mt_vector<mt_int> vtx;
    read_vtx(vtx, vtxfile);
    mt_mask vtxnod(mesh.xyz.size() / 3);
    vtxnod.insert(vtx.begin(), vtx.end());

    std::cout << "Generating surface data .. " << std::endl;

    mt_mapping<mt_int> ele2face;
    MT_MAP<triple<mt_int>, mt_int> face_map;
    MT_USET<triple<mt_int> > sel_faces;

    compute_faces(mesh, ele2face, face_map);

    for(auto it = face_map.begin(); it != face_map.end(); ++it)
    {
      const triple<mt_int> & trp = it->first;
      if(vtxnod.count(trp.v1) && vtxnod.count(trp.v2) && vtxnod.count(trp.v3))
        sel_faces.insert(trp);
    }

    size_t nelems = sel_faces.size();
    surf.e2n_cnt.assign(nelems, mt_int(3));
    surf.etype  .assign(nelems, Tri);
    surf.e2n_con.resize(nelems*3);
    eidx.resize(nelems);

    size_t widx = 0;
    for(const triple<mt_int> & t : sel_faces) {
      surf.e2n_con[widx+0] = t.v1;
      surf.e2n_con[widx+1] = t.v2;
      surf.e2n_con[widx+2] = t.v3;
      widx += 3;
    }
    bucket_sort_offset(surf.e2n_cnt, surf.e2n_dsp);

    std::cout << "Recovering volumetric element indices .. " << std::endl;

    #ifdef OPENMP
    #pragma omp parallel for schedule(dynamic, 10)
    #endif
    for(size_t e = 0; e<nelems; e++) {
      std::set<mt_int> elemset;
      mt_int * con = surf.e2n_con.data() + surf.e2n_dsp[e];

      elements_with_face(mesh, con[0], con[1], con[2], elemset);

      switch(elemset.size()) {
        case 0:
          fprintf(stderr, "%s error: Surf element does not belong to mesh! Aborting!\n", __func__);
          exit(1);

        default:
        case 1: {
          eidx[e] = *elemset.begin();
          mt_int v1 = con[0], v2 = con[1], v3 = con[2];

          vec3r ctr    = barycenter(mesh, eidx[e]);
          vec3r ctr_tr = triangle_centerpoint(v1, v2, v3, mesh.xyz);
          vec3r n      = triangle_normal(v1, v2, v3, mesh.xyz);

          if(n.scaProd(unit_vector(ctr - ctr_tr)) > 0.0) {
            con[0] = v3, con[1] = v2, con[2] = v1;
          }
          break;
        }
      }
    }
  }
  else {
    fprintf(stderr, "%s error: Neither .surf, .vtx nor .neubc file can be found! Aborting!\n", __func__);
    exit(1);
  }
}

void write_surf_info(const mt_meshdata & surface, const nbc_data* nbc,
                     const size_t numele, const size_t npts, std::string bname)
{
  fixBasename(bname);
  std::string surfname = bname + SURF_EXT;
  std::string vtxname  = bname + SURF_EXT + VTX_EXT;
  std::string nbcname  = bname + NBC_EXT;

  // generate unique nodelist from surface
  mt_vector<mt_int> vtx;
  vtx.assign(surface.e2n_con.begin(), surface.e2n_con.end());
  binary_sort(vtx); unique_resize(vtx);

  std::cout << "Writing surface " << surfname << std::endl;
  write_surf(surface.e2n_cnt, surface.e2n_con, surfname);

  std::cout << "Writing vertices " << vtxname << std::endl;
  write_vtx(vtx, vtxname);

  if(nbc) {
    std::cout << "Writing neubc " << nbcname << std::endl;
    write_nbc(*nbc, surface.e2n_con, npts, numele, nbcname);
  }
}
