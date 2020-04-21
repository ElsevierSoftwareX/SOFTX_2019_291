#include "io_utils.h"
#include "mesh_generation.h"

void mesh_to_tetgenio(const mt_meshdata & mesh, tetgenio & io)
{
  io.firstnumber = 0;
  io.numberofpoints = mesh.xyz.size() / 3;
  io.pointlist = new REAL[mesh.xyz.size()];

  for(size_t i=0; i<mesh.xyz.size(); i++) io.pointlist[i] = mesh.xyz[i];
  io.numberoffacets = mesh.e2n_cnt.size();
  io.facetlist = new tetgenio::facet[io.numberoffacets];

  for(size_t eidx = 0, ridx=0; eidx < mesh.e2n_cnt.size(); eidx++)
  {
    tetgenio::facet* f = &io.facetlist[eidx];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;

    tetgenio::polygon* p = &f->polygonlist[0];
    p->numberofvertices = mesh.e2n_cnt[eidx];
    p->vertexlist = new int[p->numberofvertices];

    for(mt_int i=0; i<mesh.e2n_cnt[eidx]; i++)
      p->vertexlist[i] = mesh.e2n_con[ridx+i];

    ridx += mesh.e2n_cnt[eidx];
  }

}

void tetgenio_to_mesh(const tetgenio & io, mt_meshdata & mesh)
{
  mesh.xyz.resize(io.numberofpoints * 3);
  for(size_t i=0; i<mesh.xyz.size(); i++) mesh.xyz[i] = io.pointlist[i];

  if(io.numberoftetrahedra == 0) {
    if(io.numberoftrifaces) {
      size_t numele = io.numberoftrifaces;
      mesh.e2n_cnt.assign(numele, mt_int(3));
      mesh.etype.assign(numele, Tri);
      mesh.etags.assign(numele, mt_int(0));

      mesh.e2n_con.resize(numele * 3);
      for(size_t i=0; i<mesh.e2n_con.size(); i++)
        mesh.e2n_con[i] = io.trifacelist[i];
    }
    else {
      fprintf(stderr, "%s error: tetgenio holds neither tets nor triangles! Aborting!\n",
              __func__);
      exit(1);
    }
  }
  else {
    size_t numele = io.numberoftetrahedra;
    mesh.e2n_cnt.assign(numele, mt_int(4));
    mesh.etype.assign(numele, Tetra);
    mesh.etags.assign(numele, mt_int(0));

    mesh.e2n_con.resize(numele * 4);
    for(size_t i=0; i<mesh.e2n_con.size(); i++)
      mesh.e2n_con[i] = io.tetrahedronlist[i];
  }
}

void mesh_with_tetgen(mt_meshdata & surfmesh,
                      mt_meshdata & volmesh,
                      mt_vector<float> & holes_xyz,
                      float edge_length,
                      float tetgen_qual,
                      mt_vector<mt_real> * size_field,
                      bool preserve_boundary)
{
  tetgenio in, out;

  mesh_to_tetgenio(surfmesh, in);

  // include holes into tetgenio
  if(holes_xyz.size()) {
    in.numberofholes = holes_xyz.size() / 3;
    in.holelist = new REAL[holes_xyz.size()];

    for(size_t i=0; i<holes_xyz.size(); i++)
      in.holelist[i] = holes_xyz[i];
  }

  char tetgen_par_str[1024];

  // include sizing field into tetgenio
  if(size_field) {
    check_size(*size_field, surfmesh.xyz.size() / 3, __func__);
    in.numberofpointmtrs = 1;
    in.pointmtrlist      = new double[size_field->size()];
    for(size_t i=0; i<size_field->size(); i++) in.pointmtrlist[i] = (*size_field)[i];
    sprintf(tetgen_par_str, "pq%fm", tetgen_qual);
  }
  else {
    float desired_volume = edge_length*edge_length*edge_length / 8.48528f;
    sprintf(tetgen_par_str, "pq%fa%f", tetgen_qual, desired_volume);
  }

  if(preserve_boundary) {
    int old_len = strlen(tetgen_par_str);
    tetgen_par_str[old_len] = 'Y';
    tetgen_par_str[old_len+1] = '\0';
  }

  tetrahedralize(tetgen_par_str, &in, &out, (tetgenio*)NULL, (tetgenio*)NULL);
  tetgenio_to_mesh(out, volmesh);

  // generate fibers. We use the same amount of directions as are stored in the
  // input surface
  int numfib = surfmesh.lon.size() == surfmesh.e2n_cnt.size() * 6 ? 2 : 1;
  printf("Generating empty fibers for output mesh (%d fiber directions) ..\n", numfib);
  volmesh.lon.assign(numfib * 3 * volmesh.e2n_cnt.size(), 0.0);
}

void generate_boundary_layers(const mt_meshdata & outer_surf,
                              const int num_layers,
                              const mt_real transition_scale,
                              const mt_real transition_inc,
                              mt_meshdata & inner_surf,
                              mt_meshdata & output_mesh)
{
  check_nonzero(outer_surf.xyz.size(), "generate_boundary_layers");

  mt_real t_scale = transition_scale;

  mt_vector<mt_real> edge_sizes;
  generate_surfmesh_sizing_field(outer_surf, edge_sizes);

  // compute inwards facing surface normals.
  inner_surf = outer_surf;
  apply_normal_orientation(inner_surf, true);

  mt_vector<mt_real> surf_normals;
  compute_nodal_surface_normals(inner_surf, inner_surf.xyz, surf_normals);

  mt_vector<mt_int> smth_nod = inner_surf.e2n_con;
  binary_sort(smth_nod); unique_resize(smth_nod);

  size_t num_tri  = inner_surf.e2n_cnt.size();
  size_t num_vert = inner_surf.xyz.size() / 3;

  mt_vector<mt_int> e2n_cnt(num_tri, 6), e2n_con(num_tri*6), tags(num_tri);
  mt_vector<elem_t> type(num_tri, Prism);

  // directional_smoothing(inner_surf, smth_nod, surf_normals, 1.0, true, 0.1, 50);
  output_mesh.xyz.append(inner_surf.xyz.begin(), inner_surf.xyz.end());

  mt_int  iter = 6;
  mt_real smth = 0.12;

  for(int layer = 0; layer < num_layers; layer++) {

    tags.assign(num_tri, layer);
    for(size_t i=0; i<num_tri; i++) {
      e2n_con[i*6+0] = inner_surf.e2n_con[i*3+0] + (num_vert * (layer+0));
      e2n_con[i*6+1] = inner_surf.e2n_con[i*3+1] + (num_vert * (layer+0));
      e2n_con[i*6+2] = inner_surf.e2n_con[i*3+2] + (num_vert * (layer+0));
      e2n_con[i*6+3] = inner_surf.e2n_con[i*3+0] + (num_vert * (layer+1));
      e2n_con[i*6+5] = inner_surf.e2n_con[i*3+1] + (num_vert * (layer+1));
      e2n_con[i*6+4] = inner_surf.e2n_con[i*3+2] + (num_vert * (layer+1));
    }
    output_mesh.e2n_cnt.append(e2n_cnt.begin(), e2n_cnt.end());
    output_mesh.e2n_con.append(e2n_con.begin(), e2n_con.end());
    output_mesh.etags.append(tags.begin(), tags.end());
    output_mesh.etype.append(type.begin(), type.end());

    for(size_t i=0; i<inner_surf.xyz.size() / 3; i++) {
      inner_surf.xyz[i*3+0] += surf_normals[i*3+0] * edge_sizes[i] * t_scale;
      inner_surf.xyz[i*3+1] += surf_normals[i*3+1] * edge_sizes[i] * t_scale;
      inner_surf.xyz[i*3+2] += surf_normals[i*3+2] * edge_sizes[i] * t_scale;
    }

    directional_smoothing(inner_surf, smth_nod, surf_normals, 0.5, 0.8, true, smth, iter);

    t_scale *= transition_inc;

    output_mesh.xyz.append(inner_surf.xyz.begin(), inner_surf.xyz.end());
    compute_nodal_surface_normals(inner_surf, inner_surf.xyz, surf_normals);
  }

  int numfib = 1;
  if(outer_surf.lon.size())
    numfib = outer_surf.lon.size() == outer_surf.e2n_cnt.size() * 6 ? 6 : 3;

  output_mesh.lon.assign(output_mesh.e2n_cnt.size() * numfib);
}

void generate_bbox_mesh(const mt_meshdata & inp_mesh, const mt_real scale, mt_meshdata & out_mesh)
{
  bbox imesh_bbox;
  cartesian_csys csys;

  mt_vector<vec3r> points;
  array_to_points(inp_mesh.xyz, points);
  generate_bbox(points, imesh_bbox);
  csys.bbox_scale_centered(imesh_bbox, scale);

  // Hex faces: (0,1,2,3) (0,4,7,1) (1,7,6,2) (2,6,5,3) (3,5,4,0) (4,5,6,7)
  // Associated triangle faces:
  // (0,1,2) (0,2,3) (0,4,7) (0,7,1) (1,7,6) (1,6,2) (2,6,5) (2,5,3) (3,5,4) (3,4,0) (4,5,6) (4,6,7)
  const vec3r pmin  = imesh_bbox.bounds[0];
  const vec3r delta = vec3r(imesh_bbox.bounds[1]) - pmin;

  points.resize(8);
  points[0] = pmin;
  points[1] = pmin + vec3r(delta.x,       0, 0);
  points[2] = pmin + vec3r(delta.x, delta.y, 0);
  points[3] = pmin + vec3r(0,       delta.y, 0);
  points[4] = pmin + vec3r(0,             0, delta.z);
  points[5] = pmin + vec3r(0,       delta.y, delta.z);
  points[6] = pmin + vec3r(delta.x, delta.y, delta.z);
  points[7] = pmin + vec3r(delta.x,       0, delta.z);

  points_to_array(points, out_mesh.xyz);
  out_mesh.e2n_cnt.assign(size_t(12), 3);
  out_mesh.etype.assign(size_t(12), Tri);
  out_mesh.etags.assign(size_t(12), mt_int(0));

  out_mesh.e2n_con.resize(size_t(12*3));
  mt_int* con = out_mesh.e2n_con.data();

  con[ 0*3+0] = 0, con[ 0*3+1] = 1, con[ 0*3+2] = 2; // (0,1,2)
  con[ 1*3+0] = 0, con[ 1*3+1] = 2, con[ 1*3+2] = 3; // (0,2,3)
  con[ 2*3+0] = 0, con[ 2*3+1] = 4, con[ 2*3+2] = 7; // (0,4,7)
  con[ 3*3+0] = 0, con[ 3*3+1] = 7, con[ 3*3+2] = 1; // (0,7,1)
  con[ 4*3+0] = 1, con[ 4*3+1] = 7, con[ 4*3+2] = 6; // (1,7,6)
  con[ 5*3+0] = 1, con[ 5*3+1] = 6, con[ 5*3+2] = 2; // (1,6,2)
  con[ 6*3+0] = 2, con[ 6*3+1] = 6, con[ 6*3+2] = 5; // (2,6,5)
  con[ 7*3+0] = 2, con[ 7*3+1] = 5, con[ 7*3+2] = 3; // (2,5,3)
  con[ 8*3+0] = 3, con[ 8*3+1] = 5, con[ 8*3+2] = 4; // (3,5,4)
  con[ 9*3+0] = 3, con[ 9*3+1] = 4, con[ 9*3+2] = 0; // (3,4,0)
  con[10*3+0] = 4, con[10*3+1] = 5, con[10*3+2] = 6; // (4,5,6)
  con[11*3+0] = 4, con[11*3+1] = 6, con[11*3+2] = 7; // (4,6,7)
}

