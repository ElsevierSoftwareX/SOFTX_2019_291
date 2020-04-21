/**
 * @file volumefractions.cpp
 * @brief standalone for computing volume fraction tags for ficticious domain permeability computation
 * @author Jana Fuchsberger
 * @version
 * @date 2019-03-28
 */
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <functional>
#include <algorithm>

#include "mt_modes_base.h"

const std::string surf_disp_par = "-surf_pos=";
const std::string mesh_disp_par = "-mesh_pos=";
const std::string comp_vel_par = "-compute_velocity=";
const std::string out_frac_fname_par = "-ofracs=";
const std::string out_vels_fname_par = "-ovels=";
const std::string use_rbf_par ="-rbf=";

static const mt_int NW = 48;
static const mt_int NQ = static_cast<mt_int>(NW * 1.5);

struct vf_options{
  std::string msh_base;
  mt_vector<std::string> surf;
  std::string ifmt;
  mt_vector<std::string> surf_xdynpt;
  std::string mesh_xdynpt;
  std::string compute_velocity;
  std::string fname_fracs;
  std::string fname_vels;
  std::string rbf;
};

void print_vf_help()
{
  fprintf(stderr, "Compute volume fraction tags for ficticious domain permeability computation.");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of volume mesh.\n", mesh_par.c_str());
  fprintf(stderr, "%s<path1>,<path2>\t (input) list of surface meshes.\n", surf_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh input format. may be: %s\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "%s<path1>,<path2>\t (optional) list of surface mesh position x.dynpt files.\n", surf_disp_par.c_str());
  fprintf(stderr, "%s<path>\t (optional) path to volume mesh positions x.dynpt file.\n", mesh_disp_par.c_str());
  fprintf(stderr, "%s<path>\t (optional) output fractions basename.\n", out_frac_fname_par.c_str());
  fprintf(stderr, "%s<path>\t (optional) output obstacle velocity basename.\n", out_vels_fname_par.c_str());
  fprintf(stderr, "%s<bool>\t (optional) flag for computing the velocity of the movement.\n", comp_vel_par.c_str());
  fprintf(stderr, "%s<bool>\t (optional) toggle RBF interpolation in moving meshes. Default off.\n", use_rbf_par.c_str());
  fprintf(stderr, "\n");
}

int vf_parse_options(int argc, char** argv, struct vf_options & opts)
{
  if(argc < 2) {
    print_vf_help();
    return 1;
  }

  // parse all parameters -----------------------------------------------------------------
  std::string surf_list, surf_xdynpt_list;
  for(int i=1; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh_par, opts.msh_base);
    if(!match) match = parse_param(param, surf_par, surf_list);
    if(!match) match = parse_param(param, inp_format_par, opts.ifmt);
    if(!match) match = parse_param(param, surf_disp_par, surf_xdynpt_list);
    if(!match) match = parse_param(param, mesh_disp_par, opts.mesh_xdynpt);
    if(!match) match = parse_param(param, comp_vel_par, opts.compute_velocity);
    if(!match) match = parse_param(param, out_frac_fname_par, opts.fname_fracs);
    if(!match) match = parse_param(param, out_vels_fname_par, opts.fname_vels);
    if(!match) match = parse_param(param, use_rbf_par, opts.rbf);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 3;
    }
  }

  // check if all relevant parameters have been set ---------------------------------------------------
  bool msh_ok = opts.msh_base.size() > 0, surf_ok = surf_list.size() > 0;


  if( !(msh_ok && surf_ok) )
  {
    std::cerr << "Error: Insufficient parameters provided." << std::endl;
    print_vf_help();
    return 2;
  }
  fixBasename(opts.msh_base);

  // prepare vector of surf base names and optionally corresponding displacement files
  split_string(surf_list, ',', opts.surf);
  for(size_t o=0; o<opts.surf.size(); o++) fixBasename(opts.surf[o]);

  if(surf_xdynpt_list.size() > 0){
    split_string(surf_xdynpt_list, ',', opts.surf_xdynpt);
    for(size_t o=0; o<opts.surf_xdynpt.size(); o++) fixBasename(opts.surf_xdynpt[o]);

    if(opts.surf_xdynpt.size() != opts.surf.size())
    {
      std::cerr << "Error: The number of surface position x.dynpt files has to match the number of surface meshes." << std::endl;
      print_vf_help();
      return 2;
    }
  }

  return 0;
}

// calculates the volume fraction of a tetraeder inside a closed surface which is represented as a kd tree
mt_real tet_enclosed_volume(const mt_vector< mt_point<mt_real> > & in_pts,
                            const mt_vector< mt_point<mt_real> > & out_pts,
                            const kdtree & tree)
{
  float t;
  bool close_cut = false;
  bool no_inters = false;
  bool did_inters;
  float rel_close_dist = 0.02f;
  mt_real enclosed_vol = 0.0;

  switch(in_pts.size())
  {
    case 1: // case 1 point inside the surface -> calculate volume of sub-tetraeder
    {
      mt_point<mt_real> frac_pts[4];
      frac_pts[0] = in_pts[0];

      for(size_t i=0; i<out_pts.size(); i++)
      {
        vec3r e = out_pts[i] - in_pts[0];
        tri_elem hit_ele;
        vec3r    hit_pos;
        did_inters = tree.closest_intersect(ray(in_pts[0], e), 0.0f+rel_close_dist, 1.0f-rel_close_dist, hit_ele, hit_pos);
        if(!did_inters) no_inters = true;
        vec3r point(in_pts[0]);
        t = (hit_pos - point).length() / e.length();
        if(t < rel_close_dist) close_cut = true;

        frac_pts[i+1]=hit_pos;
      }
      if(close_cut||no_inters) enclosed_vol = 0.0;
      else enclosed_vol = tet_volume(frac_pts);
      break;
    }
    case 2: // case 2 points inside -> calculate prism volume
    {
      mt_point<mt_real> frac_pts[6];
      /*
       * We cut the tet so that in_pts are together on the inside and out_pts are together in the outside
       *
       *        out[0]
       *        . .  .
       *      .    .   .
       *     .      .   out[1]
       * - -.- - - - .- . - - -
       *   .          . .
       * in[0]. . . . . in[1]
       */
      vec3r p00, e00 = out_pts[0] - in_pts[0];
      vec3r p10, e10 = out_pts[1] - in_pts[0];
      vec3r p01, e01 = out_pts[0] - in_pts[1];
      vec3r p11, e11 = out_pts[1] - in_pts[1];

      tri_elem hit_ele;
      did_inters = tree.closest_intersect(ray(in_pts[0], e00), 0.0f+rel_close_dist, 1.0f-rel_close_dist, hit_ele, p00);
      if(!did_inters) no_inters = true;
      vec3r point(in_pts[0]);
      t = (p00 - point).length() / e00.length();
      if(t < rel_close_dist) close_cut = true;

      did_inters = tree.closest_intersect(ray(in_pts[0], e10), 0.0f+rel_close_dist, 1.0f-rel_close_dist, hit_ele, p10);
      if(!did_inters) no_inters = true;
      t = (p10 - point).length() / e10.length();
      if(t < rel_close_dist) close_cut = true;

      did_inters = tree.closest_intersect(ray(in_pts[1], e01), 0.0f+rel_close_dist, 1.0f-rel_close_dist, hit_ele, p01);
      if(!did_inters) no_inters = true;
      point = in_pts[1];
      t = (p01 - point).length() / e01.length();
      if(t < rel_close_dist) close_cut = true;

      did_inters = tree.closest_intersect(ray(in_pts[1], e11), 0.0f+rel_close_dist, 1.0f-rel_close_dist, hit_ele, p11);
      if(!did_inters) no_inters = true;
      t = (p11 - point).length() / e11.length();
      if(t < rel_close_dist) close_cut = true;

      if(close_cut||no_inters) enclosed_vol = 0.0;
      else
      {
        frac_pts[0] = in_pts[0], frac_pts[1] = p00, frac_pts[2] = p10;
        frac_pts[3] = in_pts[1], frac_pts[4] = p01, frac_pts[5] = p11;
        enclosed_vol = prism_volume(frac_pts);
      }
      break;
    }
    case 3: // case 3 points inside the surface ->  calculate tetraeder volume - sub-tetraeder volume
    {
      vec3r frac_pts[4];
      vec3r ele_pts[4];
      frac_pts[0] = out_pts[0];
      ele_pts[0] = out_pts[0];

      for(size_t i=0; i<in_pts.size(); i++)
      {
        vec3r e = in_pts[i] - out_pts[0];
        tri_elem hit_ele;
        vec3r    hit_pos;

        did_inters = tree.closest_intersect(ray(out_pts[0], e), 0.0f+rel_close_dist, 1.0f-rel_close_dist, hit_ele, hit_pos);
        if(!did_inters) no_inters = true;
        vec3r point(out_pts[0]);
        t = (hit_pos - point).length() / e.length();
        if(t < rel_close_dist||no_inters) close_cut = true;

        frac_pts[i+1] = hit_pos;
        ele_pts[i+1] =  in_pts[i];
      }

      if(close_cut) enclosed_vol = tet_volume(ele_pts)*0.9999f; //bugfix to prevent getting a normalized volfrac > 1
      else enclosed_vol = tet_volume(ele_pts)-tet_volume(frac_pts);
      break;
    }
  }
  return enclosed_vol;
}

inline void sort_in_out(const bool is_in,
                        const mt_point<mt_real> point,
                        mt_vector< mt_point<mt_real> > & in_pts,
                        mt_vector< mt_point<mt_real> > & out_pts)
{
  if(is_in) in_pts.push_back(point);
  else out_pts.push_back(point);
}

// calculates the volume fraction of a pyramid inside a closed surface
// by splitting the pyramid into 4 tetraeder and using tet_enclosed_volume
mt_real pyr_enclosed_volume(mt_point<mt_real> * points,
                            mt_vector<bool> in_surf,
                            const kdtree & tree)
{
  // we construct a tetraeder for each edge. the third point of the tet base is the
  // center point of the pyramid's base (fc)
  mt_point<mt_real> fc = (points[0] + points[1] + points[2] + points[3]) * mt_real(0.25);
  bool fc_in_surf = inside_closed_surface(tree,fc);

  mt_real enclosed_vol = 0.0;
  mt_vector< mt_point<mt_real> > in_pts, out_pts;

  //if at least one point of the sub-tet is inside the closed surface the
  //tet volume fraction is calculated and added to the pyramid volume fraction

  // sub-tetraeder consisting of edge (0,1) center point of the pyramids base (fc) and the pyramids tip (4)
  if(in_surf[0]||in_surf[1]||in_surf[4]||fc_in_surf)
  {
    sort_in_out(in_surf[0], points[0], in_pts, out_pts);
    sort_in_out(in_surf[1], points[1], in_pts, out_pts);
    sort_in_out(fc_in_surf, fc, in_pts, out_pts);
    sort_in_out(in_surf[4], points[4], in_pts, out_pts);

    enclosed_vol += tet_enclosed_volume(in_pts, out_pts, tree);
  }

  in_pts.resize(0), out_pts.resize(0);
  // sub-tet: (1,2) (fc) (4)
  if(in_surf[1]||in_surf[2]||in_surf[4]||fc_in_surf)
  {
    sort_in_out(in_surf[1], points[1], in_pts, out_pts);
    sort_in_out(in_surf[2], points[2], in_pts, out_pts);
    sort_in_out(fc_in_surf, fc, in_pts, out_pts);
    sort_in_out(in_surf[4], points[4], in_pts, out_pts);

    enclosed_vol += tet_enclosed_volume(in_pts, out_pts, tree);
  }

  in_pts.resize(0), out_pts.resize(0);
  // sub-tet: (2,3) (fc) (4)
  if(in_surf[2]||in_surf[3]||in_surf[4]||fc_in_surf)
  {
    sort_in_out(in_surf[2], points[2], in_pts, out_pts);
    sort_in_out(in_surf[3], points[3], in_pts, out_pts);
    sort_in_out(fc_in_surf, fc, in_pts, out_pts);
    sort_in_out(in_surf[4], points[4], in_pts, out_pts);

    enclosed_vol += tet_enclosed_volume(in_pts, out_pts, tree);
  }

  in_pts.resize(0), out_pts.resize(0);
  // sub-tet: (3,0) (fc) (4)
  if(in_surf[3]||in_surf[0]||in_surf[4]||fc_in_surf)
  {
    sort_in_out(in_surf[0], points[0], in_pts, out_pts);
    sort_in_out(in_surf[3], points[3], in_pts, out_pts);
    sort_in_out(fc_in_surf, fc, in_pts, out_pts);
    sort_in_out(in_surf[4], points[4], in_pts, out_pts);

    enclosed_vol += tet_enclosed_volume(in_pts, out_pts, tree);
  }

  return enclosed_vol;
}

// calculates the volume fraction of a prism inside a closed surface
// by splitting the pyramid into 3 pyramids and 2 tetraeder
// and consequently using tet_enclosed_volume and pyr_enclosed_volume
mt_real prism_enclosed_volume(mt_point<mt_real> * prism_pts,
                              mt_vector<bool> inside,
                              const kdtree & tree)
{
  mt_point<mt_real> ctr(0,0,0);
  for(mt_int i=0; i < 6; i++) ctr += prism_pts[i];
  ctr *= (1. / 6.);
  bool ctr_in_surf = inside_closed_surface(tree,ctr);

  mt_real enclosed_vol = 0.0;
  mt_vector<bool> pyr_inside(5);
  mt_vector<bool> tet_inside(4);

  mt_point<mt_real> pyr_pts[5];
  pyr_pts[4] = ctr; pyr_inside[4] = ctr_in_surf;

  mt_point<mt_real> tet_pts[4];
  tet_pts[3] = ctr; tet_inside[3] = ctr_in_surf;

  // we construct a pyramid for each quad face. the pyramid tip is the centerpoint.
  // face (0,3,5,1)
  pyr_pts[0] = prism_pts[0]; pyr_inside[0] = inside[0];
  pyr_pts[1] = prism_pts[3]; pyr_inside[1] = inside[3];
  pyr_pts[2] = prism_pts[5]; pyr_inside[2] = inside[5];
  pyr_pts[3] = prism_pts[1]; pyr_inside[3] = inside[1];
  enclosed_vol += pyr_enclosed_volume(pyr_pts, pyr_inside, tree);

  // face (1,5,4,2)
  pyr_pts[0] = prism_pts[1]; pyr_inside[0] = inside[1];
  pyr_pts[1] = prism_pts[5]; pyr_inside[1] = inside[5];
  pyr_pts[2] = prism_pts[4]; pyr_inside[2] = inside[4];
  pyr_pts[3] = prism_pts[2]; pyr_inside[3] = inside[2];
  enclosed_vol += pyr_enclosed_volume(pyr_pts, pyr_inside, tree);

  // face (0,2,4,3)
  pyr_pts[0] = prism_pts[0]; pyr_inside[0] = inside[0];
  pyr_pts[1] = prism_pts[2]; pyr_inside[1] = inside[2];
  pyr_pts[2] = prism_pts[4]; pyr_inside[2] = inside[4];
  pyr_pts[3] = prism_pts[3]; pyr_inside[3] = inside[3];
  enclosed_vol += pyr_enclosed_volume(pyr_pts, pyr_inside, tree);

  // we construct a tet for each tri face. the tet tip is the centerpoint.
  mt_vector< mt_point<mt_real> > in_pts, out_pts;
  // face (0,1,2)
  tet_pts[0] = prism_pts[0]; tet_inside[0] = inside[0];
  tet_pts[1] = prism_pts[1]; tet_inside[1] = inside[1];
  tet_pts[2] = prism_pts[2]; tet_inside[2] = inside[2];

  //if at least one point of the sub-tet is inside the closed surface
  //the tet volume fraction is calculated and added to the pyramid volume fraction
  if(tet_inside[0]||tet_inside[1]||tet_inside[2]||tet_inside[3])
  {
    for(mt_int i=0; i<4; i++) sort_in_out(tet_inside[i], tet_pts[i], in_pts, out_pts);
    enclosed_vol += tet_enclosed_volume(in_pts, out_pts, tree);
  }

  in_pts.resize(0), out_pts.resize(0);
  // face (3,4,5)
  tet_pts[0] = prism_pts[3]; tet_inside[0] = inside[3];
  tet_pts[1] = prism_pts[4]; tet_inside[1] = inside[4];
  tet_pts[2] = prism_pts[5]; tet_inside[2] = inside[5];

  if(tet_inside[0]||tet_inside[1]||tet_inside[2]||tet_inside[3])
  {
    for(mt_int i=0; i<4; i++) sort_in_out(tet_inside[i], tet_pts[i], in_pts, out_pts);
    enclosed_vol += tet_enclosed_volume(in_pts, out_pts, tree);
  }
  return enclosed_vol;
}

// calculates the volume fraction of a hexaeder inside a closed surface
// by splitting the hexaeder into 6 pyramids and consequently using pyr_enclosed_volume
mt_real hex_enclosed_volume(mt_point<mt_real> * hexa_pts,
                            mt_vector<bool> inside,
                            const kdtree & tree)
{
  mt_point<mt_real> ctr(0,0,0);
  for(mt_int i=0; i < 8; i++) ctr += hexa_pts[i];
  ctr *= (1. / 8.);

  mt_real enclosed_vol = 0.0;
  mt_point<mt_real> pyr_pts[5];
  mt_vector<bool> pyr_inside(5);

  pyr_pts[4] = ctr;
  pyr_inside[4] = inside_closed_surface(tree,ctr);

  // we construct a pyramid for each face. the pyramid tip is the centerpoint.
  // face (0,1,2,3)
  pyr_pts[0] = hexa_pts[0]; pyr_inside[0] = inside[0];
  pyr_pts[1] = hexa_pts[1]; pyr_inside[1] = inside[1];
  pyr_pts[2] = hexa_pts[2]; pyr_inside[2] = inside[2];
  pyr_pts[3] = hexa_pts[3]; pyr_inside[3] = inside[3];
  enclosed_vol += pyr_enclosed_volume(pyr_pts, pyr_inside, tree);

  // face (0,4,7,1)
  pyr_pts[0] = hexa_pts[0]; pyr_inside[0] = inside[0];
  pyr_pts[1] = hexa_pts[4]; pyr_inside[1] = inside[4];
  pyr_pts[2] = hexa_pts[7]; pyr_inside[2] = inside[7];
  pyr_pts[3] = hexa_pts[1]; pyr_inside[3] = inside[1];
  enclosed_vol += pyr_enclosed_volume(pyr_pts, pyr_inside,tree);

  // face (1,7,6,2)
  pyr_pts[0] = hexa_pts[1]; pyr_inside[0] = inside[1];
  pyr_pts[1] = hexa_pts[7]; pyr_inside[1] = inside[7];
  pyr_pts[2] = hexa_pts[6]; pyr_inside[2] = inside[6];
  pyr_pts[3] = hexa_pts[2]; pyr_inside[3] = inside[2];
  enclosed_vol += pyr_enclosed_volume(pyr_pts, pyr_inside, tree);

  // face (2,6,5,3)
  pyr_pts[0] = hexa_pts[2]; pyr_inside[0] = inside[2];
  pyr_pts[1] = hexa_pts[6]; pyr_inside[1] = inside[6];
  pyr_pts[2] = hexa_pts[5]; pyr_inside[2] = inside[5];
  pyr_pts[3] = hexa_pts[3]; pyr_inside[3] = inside[3];
  enclosed_vol += pyr_enclosed_volume(pyr_pts, pyr_inside, tree);

  // face (3,5,4,0)
  pyr_pts[0] = hexa_pts[3]; pyr_inside[0] = inside[3];
  pyr_pts[1] = hexa_pts[5]; pyr_inside[1] = inside[5];
  pyr_pts[2] = hexa_pts[4]; pyr_inside[2] = inside[4];
  pyr_pts[3] = hexa_pts[0]; pyr_inside[3] = inside[0];
  enclosed_vol += pyr_enclosed_volume(pyr_pts, pyr_inside, tree);

  // face (4,5,6,7)
  pyr_pts[0] = hexa_pts[4]; pyr_inside[0] = inside[4];
  pyr_pts[1] = hexa_pts[5]; pyr_inside[1] = inside[5];
  pyr_pts[2] = hexa_pts[6]; pyr_inside[2] = inside[6];
  pyr_pts[3] = hexa_pts[7]; pyr_inside[3] = inside[7];
  enclosed_vol += pyr_enclosed_volume(pyr_pts, pyr_inside, tree);

  return enclosed_vol;
}

// calculates normalized volume fractions of tet, hexa, pyr and prism elements inside a closed surface
void volume_fractions(const mt_meshdata & mesh,
                      const mt_meshdata & surf,
                      mt_vector<mt_real> & fractions)
{
  fractions.resize(mesh.e2n_cnt.size());
  kdtree tree(50);
  tree.build_tree(surf);
  MT_USET<mt_int> split_elem;

  // determine if nodes are inside or outside closed surface
  mt_vector<bool> inside;
  nodes_in_surface(mesh,tree,inside);

  // loop over all elements to sort into inside, outside and split
  #ifdef OPENMP
  #pragma omp parallel for schedule(guided)
  #endif
  for(unsigned int elem=0; elem<mesh.e2n_cnt.size(); elem++) {
    mt_int start=mesh.e2n_dsp[elem], stop = start + mesh.e2n_cnt[elem];
    mt_int in_count=0;
    for(mt_int i=start; i<stop; i++){
      mt_int node = mesh.e2n_con[i];
      if (inside[node]) in_count++;
    }
    if(in_count==mesh.e2n_cnt[elem]) fractions[elem] = 1.0; // set fraction value for element inside to 1
    else if(in_count==0) fractions[elem] = 0.0; // set fraction value for element outside to o
    else {
      #ifdef OPENMP
      #pragma omp critical
      #endif
      split_elem.insert(elem); // remember all split elements
    }
  }

  mt_vector<mt_int> temp;
  temp.assign(split_elem.begin(), split_elem.end());

  // loop over split elements to calculate volume fractions
  #ifdef OPENMP
  #pragma omp parallel for schedule(guided)
  #endif
  for(size_t eit = 0; eit < temp.size(); eit++)
  //for(auto eit = split_elem.begin(); eit != split_elem.end(); ++eit)
  {
    mt_int eidx = temp[eit];
    // depending on the element type the element is split up in sub-tets to calculate the volume fraction
    switch(mesh.etype[eidx]) {
      case Tri:
      case Quad:
        fprintf(stderr, "Error: Unsupported element!\n");
        exit(EXIT_FAILURE);
      case Tetra:
      {
        mt_vector< mt_point<mt_real> > in_pts, out_pts;
        mt_point<mt_real> tet_pts[4];
        mt_int start=mesh.e2n_dsp[eidx], stop = start + mesh.e2n_cnt[eidx];
        // sort nodes of current tet element into inside and outside nodes
        // and store the corresponding coordinates for the volume fraction calculation
        for(mt_int i=start; i<stop; i++)
        {
          mt_int node = mesh.e2n_con[i];
          vec3r point(mesh.xyz.data() + node * 3);
          if(inside[node]) in_pts.push_back(point);
          else out_pts.push_back(point);
          tet_pts[i-start] = point;
        }

        mt_real enclosed_vol = tet_enclosed_volume(in_pts, out_pts, tree);
        mt_real tet_vol = tet_volume(tet_pts);
        // calculate volume fraction
        fractions[eidx] = clamp(enclosed_vol/tet_vol,0.0,1.0, true);
        break;
      }
      case Pyramid:
      {
        mt_vector<bool> in_surf;
        mt_point<mt_real> pyr_pts[5];
        mt_int start=mesh.e2n_dsp[eidx], stop = start + mesh.e2n_cnt[eidx];
        // get coords of nodes of current pyr element
        for(mt_int i=start; i<stop; i++)
        {
          mt_int node = mesh.e2n_con[i];
          vec3r point(mesh.xyz.data() + node * 3);
          if(inside[node]) in_surf.push_back(true);
          else in_surf.push_back(false);
          pyr_pts[i-start] = point;
        }

        mt_real enclosed_vol = pyr_enclosed_volume(pyr_pts, in_surf, tree);
        mt_real pyr_vol = pyr_volume(pyr_pts);
        // calculate volume fraction
        fractions[eidx] = clamp(enclosed_vol/pyr_vol,0.0,1.0, true);
        break;
      }
      case Prism:
      {
        mt_vector<bool> in_surf;
        mt_point<mt_real> prism_pts[6];
        mt_int start=mesh.e2n_dsp[eidx], stop = start + mesh.e2n_cnt[eidx];
        // get coords of nodes of current prism element
        for(mt_int i=start; i<stop; i++)
        {
          mt_int node = mesh.e2n_con[i];
          vec3r point(mesh.xyz.data() + node * 3);
          if(inside[node]) in_surf.push_back(true);
          else in_surf.push_back(false);
          prism_pts[i-start] = point;
        }

        mt_real enclosed_vol = prism_enclosed_volume(prism_pts, in_surf, tree);
        mt_real prism_vol = prism_volume(prism_pts);
        // calculate volume fraction
        fractions[eidx] = clamp(enclosed_vol/prism_vol,0.0,1.0, true);
        break;
      }
      case Hexa:
      {
        mt_vector<bool> in_surf;
        mt_point<mt_real> hexa_pts[8];
        mt_int start=mesh.e2n_dsp[eidx], stop = start + mesh.e2n_cnt[eidx];
        // get coords of nodes of current hexa element
        for(mt_int i=start; i<stop; i++)
        {
          mt_int node = mesh.e2n_con[i];
          vec3r point(mesh.xyz.data() + node * 3);
          if(inside[node]) in_surf.push_back(true);
          else in_surf.push_back(false);
          hexa_pts[i-start] = point;
        }

        mt_real enclosed_vol = hex_enclosed_volume(hexa_pts, in_surf, tree);
        mt_real hex_vol = hex_volume(hexa_pts);
        // calculate volume fraction
        fractions[eidx] = clamp(enclosed_vol/hex_vol,0.0,1.0, true);
        break;
      }
      default:
        fprintf(stderr, "Error: Unsupported element!\n");
        exit(EXIT_FAILURE);
    }
  }
}

//combines volumefractions and velocities of all provided surfaces
void combine_outputs(const mt_meshdata & mesh,
                     mt_vector<mt_real> & fracs,
                     const mt_vector<mt_vector<mt_real>> & all_fracs,
                     mt_vector<mt_real> & vels,
                     const mt_vector<mt_vector<mt_real>> & all_vels)
{
  mt_int n_elems=mesh.e2n_cnt.size();
  fracs.assign(all_fracs[0].size(), 0.0);
  vels.assign(all_vels[0].size(), 0.0);

  for(mt_int elem=0; elem<n_elems; elem++)
  {
    for(size_t o=0; o<all_fracs.size(); o++)
    {
      if(all_fracs[o][elem]>fracs[elem]) fracs[elem] = all_fracs[o][elem];
    }

    if(fracs[elem]>0.0)
    {
      mt_int start=mesh.e2n_dsp[elem], stop = start + mesh.e2n_cnt[elem];
      for(mt_int i=start; i<stop; i++)
      {
        mt_int node = mesh.e2n_con[i];
        mt_real * vels_ptr = vels.data() + 3*node;
        for(size_t o=0; o<all_vels.size(); o++)
        {
          const mt_real * all_vels_ptr = all_vels[o].data() + 3*node;

          vels_ptr[0] += all_vels_ptr[0];
          vels_ptr[1] += all_vels_ptr[1];
          vels_ptr[2] += all_vels_ptr[2];
        }
        vels_ptr[0] /= all_vels.size();
        vels_ptr[1] /= all_vels.size();
        vels_ptr[2] /= all_vels.size();
      }
    }
  }
}

//combines volumefractions of all provided surfaces
void combine_outputs(mt_vector<mt_real> & fracs,
                     const mt_vector<mt_vector<mt_real>> & all_fracs)
{
  mt_int n_elems=all_fracs[0].size();
  fracs.assign(n_elems, 0.0);

  for(mt_int elem=0; elem<n_elems; elem++)
  {
    for(size_t o=0; o<all_fracs.size(); o++)
    {
      if(all_fracs[o][elem]>fracs[elem]) fracs[elem] = all_fracs[o][elem];
    }
  }
}

//calculates volumefractions and velocities for the case of moving surfaces and a moving mesh
//and writes them to .igb files
void calc_fracs_moving_mesh_moving_surface(mt_meshdata & mesh,
                                           mt_vector<mt_meshdata> & surf,
                                           const vf_options opts,
                                           const mt_real scale_factor)
{
  igb_header mesh_igb, frac_igb, vel_igb;
  mt_vector<igb_header> surf_igb(opts.surf_xdynpt.size());
  for(size_t o=0; o<opts.surf_xdynpt.size(); o++){
    init_igb_header(opts.surf_xdynpt[o], surf_igb[o]);
    read_igb_header(surf_igb[o]);
  }

  init_igb_header(opts.mesh_xdynpt, mesh_igb);
  read_igb_header(mesh_igb);

  const bool comp_vel = opts.compute_velocity.size() > 0 ? static_cast<bool>(atoi(opts.compute_velocity.c_str())) : true;
  const bool use_rbf = opts.rbf.size() > 0;

  frac_igb = mesh_igb;
  frac_igb.filename = opts.fname_fracs.size() > 0 ? opts.fname_fracs + IGB_EXT : opts.msh_base + ".fracs" + IGB_EXT;
  frac_igb.v_x = mesh.e2n_cnt.size();
  set_igb_header_datatype("float", frac_igb);
  frac_igb.fileptr  = NULL;
  write_igb_header(frac_igb);
  std::vector<std::vector<mt_real> > mesh_pos, surf_pos;
  std::vector<std::vector<float> > igb_fracs(1), igb_vels(1);
  igb_fracs[0].resize(mesh.e2n_cnt.size());

  const int dpn = 3;
  PointCloudInterpolator* cloud_interp = nullptr;

  mt_meshdata mesh_base(mesh);
  mt_vector<mt_real> mesh_dsp(mesh.xyz.size());
  mt_vector<mt_real> interp_mesh_disp;
  mt_vector<mt_real> surf_vel, interp_surf_vel(mesh.xyz.size());
  mt_vector<mt_vector<mt_real>> all_vels(opts.surf.size());

  if(comp_vel) {
    igb_vels[0].resize(mesh.xyz.size());
    vel_igb = mesh_igb;
    vel_igb.filename = opts.fname_vels.size() > 0 ? opts.fname_vels + IGB_EXT : opts.msh_base + ".vels" + IGB_EXT;
    vel_igb.v_x = mesh.n2e_cnt.size();
    set_igb_header_datatype("vec3f", vel_igb);
    vel_igb.fileptr  = NULL;
    write_igb_header(vel_igb);
    cloud_interp = new PointCloudInterpolator;
  }


  mt_vector<mt_meshdata> surf_base(opts.surf.size());
  mt_vector<mt_vector<mt_real>> old_surf_xyz(opts.surf.size());
  mt_vector<mt_vector<mt_real>> old_old_surf_xyz(opts.surf.size());

  float dt = surf_igb[0].v_inc_t;
  size_t tsteps = surf_igb[0].v_t;
  size_t N_surf_disp = size_t(surf_igb[0].v_x);

  for(size_t o=0; o<opts.surf.size(); o++){
    check_size(tsteps, surf_igb[o].v_t, "Surface mesh displacement calculation: #timesteps has to be equal for all obstacles.");
    check_size(dt, surf_igb[o].v_inc_t, "Surface mesh displacement calculation: dt has to be equal for all obstacles.");
    check_size(surf_igb[o].v_x, surf[o].xyz.size() / 3, "Surface mesh displacement calculation: spacial dimensions of displacement and surface base mesh have to match for all obstacles.");
    surf_base[o] = surf[o];
    old_surf_xyz[o].assign(surf[o].xyz.begin(), surf[o].xyz.end());
    old_old_surf_xyz[o].assign(surf[o].xyz.begin(), surf[o].xyz.end());
  }

  size_t N_mesh_disp = size_t(mesh_igb.v_x);
  check_size(N_mesh_disp, mesh.xyz.size() / 3, "Volume mesh displacement calculation");


  mt_vector<mt_point<mt_real> > idat_vec, odat_vec;
  mt_vector<mt_int> corr;
  mt_vector<mt_real> corr_dist, idat(mesh.e2n_cnt.size());
  mt_manifold imnfld, omnfld;
  imnfld.on_mnfld.resize(mesh.xyz.size() / 3),
  omnfld.on_mnfld.resize(surf[0].xyz.size() / 3);

  PROGRESS<mt_int> progress(tsteps,"Computing time dependent volume fractions:");
  for(size_t t=0; t<tsteps; t++)
  {
    progress.next();
    read_igb_block(mesh_pos, 1, mesh_igb);
    // Apply scale factor to all entries
    std::transform(mesh_pos[0].begin(), mesh_pos[0].end(), mesh_pos[0].begin(),
                     std::bind1st(std::multiplies<mt_real>(), scale_factor));
    mesh.xyz.assign(mesh_pos[0].begin(), mesh_pos[0].end());
    // displace the mesh
    for(size_t j=0; j<N_mesh_disp; j++)
    {
      mt_real const * _xyz_new = mesh.xyz.data() + 3*j;
      mt_real const * _xyz_base = mesh_base.xyz.data() + 3*j;
      mt_real * _disp = mesh_dsp.data() + 3*j;
      _disp[0] = _xyz_new[0] -_xyz_base[0];
      _disp[1] = _xyz_new[1] -_xyz_base[1];
      _disp[2] = _xyz_new[2] -_xyz_base[2];
    }

    mt_vector<mt_real> fracs;
    mt_vector<mt_vector<mt_real>> all_fracs(opts.surf.size());

    for(size_t o=0; o<opts.surf.size(); o++) //loop over all obstacles
    {
      if(comp_vel) {
        old_old_surf_xyz[o].assign(old_surf_xyz[o].begin(), old_surf_xyz[o].end());
        old_surf_xyz[o].assign(surf[o].xyz.begin(), surf[o].xyz.end());
      }
      read_igb_block(surf_pos, 1, surf_igb[o]);
      // Apply scale factor to all entries
      std::transform(surf_pos[0].begin(), surf_pos[0].end(), surf_pos[0].begin(),
                       std::bind1st(std::multiplies<mt_real>(), scale_factor));
      surf[o].xyz.assign(surf_pos[0].begin(), surf_pos[0].end());
      if(t == 0) {
        old_old_surf_xyz[o].assign(surf[o].xyz.begin(), surf[o].xyz.end());
        old_surf_xyz[o].assign(surf[o].xyz.begin(), surf[o].xyz.end());
      }
      compute_correspondance(surf_base[o], mesh_base, corr, corr_dist);
      idat.assign(mesh_dsp.begin(), mesh_dsp.end());
      array_to_points(idat, idat_vec);
      nodal_interpolation(mesh_base, surf_base[o], imnfld, omnfld, corr, corr_dist, idat_vec, odat_vec);
      interp_mesh_disp.resize(surf[o].xyz.size());
      points_to_array(odat_vec, interp_mesh_disp);

      //displace the surface considering the mesh displacement
      surf_vel.resize(surf[0].xyz.size(), 0.0);
      interp_surf_vel.resize(mesh.xyz.size(), 0.0);
      N_surf_disp = (size_t)surf_igb[o].v_x;
      for(size_t j=0; j<N_surf_disp; j++)
      {
        mt_real * _surf_xyz_ptr = surf[o].xyz.data() + 3*j;
        mt_real const * _interp_mesh_disp_ptr = interp_mesh_disp.data() + 3*j;
        _surf_xyz_ptr[0] += _interp_mesh_disp_ptr[0];
        _surf_xyz_ptr[1] += _interp_mesh_disp_ptr[1];
        _surf_xyz_ptr[2] += _interp_mesh_disp_ptr[2];
      }
      if(t == 0) {
        old_old_surf_xyz[o].assign(surf[o].xyz.begin(), surf[o].xyz.end());
        old_surf_xyz[o].assign(surf[o].xyz.begin(), surf[o].xyz.end());
      }
      if(comp_vel && t > 0) {
        for(size_t j=0; j<N_surf_disp; j++) {
          mt_real * _surf_vel_ptr = surf_vel.data() + 3*j;
          mt_real const * _surf_xyz_ptr = surf[o].xyz.data() + 3*j;
          mt_real const * _old_xyz_ptr = old_surf_xyz[o].data() + 3*j;
          mt_real const * _old_old_xyz_ptr = old_old_surf_xyz[o].data() + 3*j;
          if(t > 1) {
            _surf_vel_ptr[0] = (3 * _surf_xyz_ptr[0] - 4 * _old_xyz_ptr[0] + _old_old_xyz_ptr[0]) / (2 * dt);
            _surf_vel_ptr[1] = (3 * _surf_xyz_ptr[1] - 4 * _old_xyz_ptr[1] + _old_old_xyz_ptr[1]) / (2 * dt);
            _surf_vel_ptr[2] = (3 * _surf_xyz_ptr[2] - 4 * _old_xyz_ptr[2] + _old_old_xyz_ptr[2]) / (2 * dt);
          } else {
            _surf_vel_ptr[0] = (_surf_xyz_ptr[0] - _old_xyz_ptr[0]) / dt;
            _surf_vel_ptr[1] = (_surf_xyz_ptr[1] - _old_xyz_ptr[1]) / dt;
            _surf_vel_ptr[2] = (_surf_xyz_ptr[2] - _old_xyz_ptr[2]) / dt;
          }
        }
      }
      volume_fractions(mesh, surf[o], fracs);
      all_fracs[o]=fracs;
      if(comp_vel && t > 0) {
        cloud_interp->set_pts(surf[o].xyz);
        cloud_interp->setup_roi(NW);
        if(use_rbf) {
        cloud_interp->setup_rbf_coeffs(surf_vel, dpn, NQ);
        cloud_interp->interpolate_data(mesh.xyz, interp_surf_vel);
        } else
          cloud_interp->interpolate_shepard(mesh.xyz, surf_vel, dpn, false, interp_surf_vel, false);
        // reverse scaling
        for(size_t k=0; k < interp_surf_vel.size(); k++)
          interp_surf_vel[k] /= scale_factor;
      }
      all_vels[o] = interp_surf_vel;
    }
    if(comp_vel)
    {
      combine_outputs(mesh, fracs, all_fracs, interp_surf_vel, all_vels);
      igb_fracs[0].assign(fracs.begin(), fracs.end());
      write_igb_block(igb_fracs, frac_igb);
      igb_vels[0].assign(interp_surf_vel.begin(), interp_surf_vel.end());
      write_igb_block(igb_vels, vel_igb);
    }
    else
    {
      combine_outputs(fracs, all_fracs);
      igb_fracs[0].assign(fracs.begin(), fracs.end());
      write_igb_block(igb_fracs, frac_igb);
    }
  }
  progress.finish();
  fclose(frac_igb.fileptr);
  for(size_t o=0; o<opts.surf_xdynpt.size(); o++) fclose(surf_igb[o].fileptr);
  fclose(mesh_igb.fileptr);
  if(comp_vel) {
    fclose(vel_igb.fileptr);
    delete cloud_interp;
  }
}

//calculates volumefractions and velocities for the case of moving surfaces and a fixed mesh
//and writes them to .igb files
void calc_fracs_fixed_mesh_moving_surface(const mt_meshdata & mesh,
                                          mt_vector<mt_meshdata> & surf,
                                          const vf_options opts,
                                          const mt_real scale_factor)
{
  mt_vector<mt_vector<mt_real>> surf_base_xyz(opts.surf.size());
  mt_vector<igb_header> surf_dynpt(opts.surf_xdynpt.size());

  for(size_t o=0; o<opts.surf_xdynpt.size(); o++){
    surf_base_xyz[o].assign(surf[o].xyz.begin(), surf[o].xyz.end());
    init_igb_header(opts.surf_xdynpt[o], surf_dynpt[o]);
    read_igb_header(surf_dynpt[o]);
  }

  float dt = surf_dynpt[0].v_inc_t;
  size_t tsteps = surf_dynpt[0].v_t;
  size_t N_surf_disp = size_t(surf_dynpt[0].v_x);

  for(size_t o=0; o<opts.surf_xdynpt.size(); o++){
    check_size(tsteps, surf_dynpt[o].v_t, "Surface mesh displacement calculation: #timesteps has to be equal for all obstacles.");
    check_size(dt, surf_dynpt[o].v_inc_t, "Surface mesh displacement calculation: dt has to be equal for all obstacles.");
    check_size(surf_dynpt[o].v_x, surf[o].xyz.size() / 3, "Surface mesh displacement calculation: spacial dimensions of displacement and surface base mesh have to match for all obstacles.");
  }

  igb_header frac_igb, vel_igb;
  frac_igb = surf_dynpt[0];
  frac_igb.filename = opts.fname_fracs.size() > 0 ? opts.fname_fracs + IGB_EXT : opts.msh_base + ".fracs" + IGB_EXT;
  frac_igb.v_x = mesh.e2n_cnt.size();
  set_igb_header_datatype("float", frac_igb);
  frac_igb.fileptr  = NULL;
  write_igb_header(frac_igb);

  // surf_pos holds the position of the surface at time
  std::vector<std::vector<mt_real> > surf_pos;
  std::vector<std::vector<float> > igb_fracs(1), igb_vels(1);
  igb_fracs[0].resize(mesh.e2n_cnt.size());
  mt_vector<mt_real> interp_surf_vel;
  mt_vector<mt_real> surf_vel;
  mt_vector<mt_vector<mt_real>> all_vels(opts.surf.size());

  mt_vector<mt_vector<mt_real>> surf_pos_n1(opts.surf.size());
  mt_vector<mt_vector<mt_real>> surf_pos_n2(opts.surf.size());
  const int dpn = 3;
  PointCloudInterpolator * cloud_interp = nullptr;

  const bool comp_vel = opts.compute_velocity.size() > 0 ? static_cast<bool>(atoi(opts.compute_velocity.c_str())) : true;
  const bool use_rbf = opts.rbf.size() > 0;
  if(comp_vel) {
    for(size_t k=0; k < opts.surf.size(); k++) {
      surf_pos_n1[k].assign(surf[k].xyz.begin(), surf[k].xyz.end());
      surf_pos_n2[k].assign(surf[k].xyz.begin(), surf[k].xyz.end());
    }
    interp_surf_vel.resize(mesh.xyz.size());
    vel_igb = surf_dynpt[0];
    vel_igb.filename = opts.fname_vels.size() > 0 ? opts.fname_vels + IGB_EXT : opts.msh_base + ".vels" + IGB_EXT;
    vel_igb.v_x = mesh.n2e_cnt.size();
    set_igb_header_datatype("vec3f", vel_igb);
    vel_igb.fileptr  = NULL;
    write_igb_header(vel_igb);
    igb_vels[0].resize(mesh.xyz.size());
    cloud_interp = new PointCloudInterpolator;
  }


  PROGRESS<mt_int> progress(tsteps,"Computing time dependent volume fractions:");
  for(size_t t=0; t<tsteps; t++)
  {
    progress.next();

    mt_vector<mt_real> fracs;
    mt_vector<mt_vector<mt_real>> all_fracs(opts.surf.size());
    for(size_t o=0; o<opts.surf.size(); o++) //loop over all obstacles
    {
      //shift time steps
      if(comp_vel) {
        surf_pos_n2[o].assign(surf_pos_n1[o].begin(), surf_pos_n1[o].end());
        surf_pos_n1[o].assign(surf[o].xyz.begin(), surf[o].xyz.end());
      }
      read_igb_block(surf_pos, 1, surf_dynpt[o]);
      // Apply scale factor to all entries
      std::transform(surf_pos[0].begin(), surf_pos[0].end(), surf_pos[0].begin(),
                     std::bind1st(std::multiplies<mt_real>(), scale_factor));
      surf[o].xyz.assign(surf_pos[0].begin(), surf_pos[0].end());
      //initial position is reference for velocity calculation
      if(t == 0) {
        surf_pos_n2[o].assign(surf[o].xyz.begin(), surf[o].xyz.end());
        surf_pos_n1[o].assign(surf[o].xyz.begin(), surf[o].xyz.end());
      }
      // compute volume fractions with new coords
      volume_fractions(mesh, surf[o], fracs);
      all_fracs[o]=fracs;
      surf_vel.resize(surf[0].xyz.size(), 0.0);
      interp_surf_vel.assign(interp_surf_vel.size(), 0.0);
      // compute velocity only for t > 0 since we assume zero velocity in the beginning
      if(t > 0 && comp_vel) {
        N_surf_disp = (size_t)(surf_dynpt[o].v_x);
        for(size_t j=0; j<N_surf_disp; j++)
        {
          mt_real const * _surf_xyz_ptr = surf[o].xyz.data() + dpn*j;
          mt_real const * _old_xyz_ptr = surf_pos_n1[o].data() + dpn*j;
          mt_real const * _old_old_xyz_ptr = surf_pos_n2[o].data() + dpn*j;
          mt_real * _surf_vel_ptr = surf_vel.data() + 3*j;
          // Higher order velocity approximation for tsteps > 1
          if(t > 1) {
            _surf_vel_ptr[0] = (3 * _surf_xyz_ptr[0] - 4 * _old_xyz_ptr[0] + _old_old_xyz_ptr[0]) / (2 * dt);
            _surf_vel_ptr[1] = (3 * _surf_xyz_ptr[1] - 4 * _old_xyz_ptr[1] + _old_old_xyz_ptr[1]) / (2 * dt);
            _surf_vel_ptr[2] = (3 * _surf_xyz_ptr[2] - 4 * _old_xyz_ptr[2] + _old_old_xyz_ptr[2]) / (2 * dt);
          } else {
            _surf_vel_ptr[0] = (_surf_xyz_ptr[0] - _old_xyz_ptr[0]) / (dt);
            _surf_vel_ptr[1] = (_surf_xyz_ptr[1] - _old_xyz_ptr[1]) / (dt);
            _surf_vel_ptr[2] = (_surf_xyz_ptr[2] - _old_xyz_ptr[2]) / (dt);
          }
        }
        cloud_interp->set_pts(surf[o].xyz);
        cloud_interp->setup_roi(NW);
        if(use_rbf) {
          cloud_interp->setup_rbf_coeffs(surf_vel, dpn, NQ);
          cloud_interp->interpolate_data(mesh.xyz, interp_surf_vel);
        } else
          cloud_interp->interpolate_shepard(mesh.xyz, surf_vel, dpn, false, interp_surf_vel, false);
        //Reverse scaling
        for(size_t k=0; k < interp_surf_vel.size(); k++)
          interp_surf_vel[k] /= scale_factor;
      }
      all_vels[o] = interp_surf_vel;
    }
    if(comp_vel)
    {
      combine_outputs(mesh, fracs, all_fracs, interp_surf_vel, all_vels);
      igb_fracs[0].assign(fracs.begin(), fracs.end());
      write_igb_block(igb_fracs, frac_igb);
      igb_vels[0].assign(interp_surf_vel.begin(), interp_surf_vel.end());
      write_igb_block(igb_vels, vel_igb);
    }
    else
    {
      combine_outputs(fracs, all_fracs);
      igb_fracs[0].assign(fracs.begin(), fracs.end());
      write_igb_block(igb_fracs, frac_igb);
    }
  }
  progress.finish();
  fclose(frac_igb.fileptr);
  for(size_t o=0; o<opts.surf.size(); o++) fclose(surf_dynpt[o].fileptr);
  if(comp_vel)
  {
    fclose(vel_igb.fileptr);
    delete cloud_interp;
  }
}

//calculates volumefractions for the case of fixed surfaces and a moving mesh
//and writes them to an .igb file
void calc_fracs_moving_mesh_fixed_surface(mt_meshdata & mesh,
                                          mt_vector<mt_meshdata> & surf,
                                          const vf_options opts,
                                          const mt_real scale_factor)
{
  igb_header mesh_igb, frac_igb, vel_igb;
  init_igb_header(opts.mesh_xdynpt, mesh_igb);
  read_igb_header(mesh_igb);
  frac_igb = mesh_igb;
  frac_igb.filename = opts.fname_fracs.size() > 0 ? opts.fname_fracs + IGB_EXT : opts.msh_base + ".fracs" + IGB_EXT;
  frac_igb.v_x = mesh.e2n_cnt.size();
  set_igb_header_datatype("float", frac_igb);
  frac_igb.fileptr  = NULL;
  write_igb_header(frac_igb);

  size_t tsteps = mesh_igb.v_t;
  float dt = mesh_igb.v_inc_t;
  size_t N_surf_disp = surf[0].xyz.size() / 3;
  size_t N_mesh_disp = size_t(mesh_igb.v_x);
  check_size(N_mesh_disp, mesh.xyz.size() / 3, "Volume mesh displacement calculation");

  std::vector<std::vector<mt_real> > mesh_pos;
  std::vector<std::vector<float> > igb_fracs(1), igb_vels(1);
  igb_fracs[0].resize(mesh.e2n_cnt.size());
  mt_vector<mt_meshdata> surf_base(opts.surf.size());
  mt_vector<mt_real> mesh_xyz_old(mesh.xyz.size());
  mt_vector<mt_real> mesh_xyz_old_old(mesh.xyz.size());
  mesh_xyz_old.assign(mesh.xyz.begin(), mesh.xyz.end());
  mt_vector<mt_real> mesh_vel(mesh.xyz.size());
  mt_vector<mt_real> mesh_dsp(mesh.xyz.size());


  for(size_t o=0; o<opts.surf.size(); o++) surf_base[o] = surf[o];
  mt_meshdata mesh_base(mesh);

  mt_vector<mt_real> interp_mesh_disp;
  mt_vector<mt_point<mt_real>> idat_vec, odat_vec;
  mt_vector<mt_int> corr;
  mt_vector<mt_real> corr_dist, idat(mesh.e2n_cnt.size());
  mt_manifold imnfld, omnfld;
  imnfld.on_mnfld.resize(mesh.xyz.size() / 3),
  omnfld.on_mnfld.resize(surf[0].xyz.size() / 3);

  const bool comp_vel = opts.compute_velocity.size() > 0 ? static_cast<bool>(atoi(opts.compute_velocity.c_str())) : true;
  if(comp_vel) {
    vel_igb = mesh_igb;
    vel_igb.filename = opts.fname_vels.size() > 0 ? opts.fname_vels + IGB_EXT : opts.msh_base + ".vels" + IGB_EXT;
    vel_igb.v_x = mesh.n2e_cnt.size();
    set_igb_header_datatype("vec3f", vel_igb);
    vel_igb.fileptr  = NULL;
    write_igb_header(vel_igb);
    igb_vels[0].resize(mesh.xyz.size());
  }
  PROGRESS<mt_int> progress(tsteps, "Computing time dependent volume fractions:");
  for(size_t t=0; t<tsteps; t++)
  {
    progress.next();
    // Shift solution by one timestep
    if(comp_vel) {
      mesh_xyz_old_old.assign(mesh_xyz_old.begin(), mesh_xyz_old.end());
      mesh_xyz_old.assign(mesh.xyz.begin(), mesh.xyz.end());
    }

    read_igb_block(mesh_pos, 1, mesh_igb);
    mt_vector<mt_real> fracs, vel;
    mt_vector<mt_vector<mt_real>> all_fracs(opts.surf.size());
    // Apply scale factor to all entries
    std::transform(mesh_pos[0].begin(), mesh_pos[0].end(), mesh_pos[0].begin(),
                     std::bind1st(std::multiplies<mt_real>(), scale_factor));
    mesh.xyz.assign(mesh_pos[0].begin(), mesh_pos[0].end());
    if(comp_vel && t == 0) {
      mesh_xyz_old_old.assign(mesh.xyz.begin(), mesh.xyz.end());
      mesh_xyz_old.assign(mesh.xyz.begin(), mesh.xyz.end());
    }
    for(size_t j=0; j<N_mesh_disp; j++)
    {
      mt_real const * _xyz_new = mesh.xyz.data() + 3*j;
      mt_real const * _xyz_base = mesh_base.xyz.data() + 3*j;
      mt_real * disp = mesh_dsp.data() + 3*j;
      disp[0] = _xyz_new[0] - _xyz_base[0];
      disp[1] = _xyz_new[1] - _xyz_base[1];
      disp[2] = _xyz_new[2] - _xyz_base[2];
    }

    for(size_t o=0; o<opts.surf.size(); o++)
    {
      compute_correspondance(surf_base[o], mesh_base, corr, corr_dist);
      idat.assign(mesh_dsp.begin(), mesh_dsp.end());
      array_to_points(idat, idat_vec);
      nodal_interpolation(mesh_base, surf_base[o], imnfld, omnfld, corr, corr_dist, idat_vec, odat_vec);
      interp_mesh_disp.resize(surf[0].xyz.size());
      points_to_array(odat_vec, interp_mesh_disp);
      N_surf_disp = surf[o].xyz.size() / 3;
      for(size_t j=0; j<N_surf_disp; j++)
      {
        mt_real * _surf_xyz_ptr = surf[o].xyz.data() + 3*j;
        mt_real const * _surf_base_xyz_ptr = surf_base[o].xyz.data() + 3*j;
        mt_real const * _interp_mesh_disp_ptr = interp_mesh_disp.data() + 3*j;

        _surf_xyz_ptr[0] = _surf_base_xyz_ptr[0] + _interp_mesh_disp_ptr[0];
        _surf_xyz_ptr[1] = _surf_base_xyz_ptr[1] + _interp_mesh_disp_ptr[1];
        _surf_xyz_ptr[2] = _surf_base_xyz_ptr[2] + _interp_mesh_disp_ptr[2];
      }

      volume_fractions(mesh, surf[o], fracs);
      all_fracs[o]=fracs;
    }
    if(comp_vel)
    {
      combine_outputs(fracs, all_fracs);
      mesh_vel.assign(mesh.xyz.size(), 0.0);
      if(t > 0) {
        for(size_t elem=0; elem<mesh.e2n_cnt.size(); elem++)
          {
            if(fracs[elem]>0.0)
            {
              mt_int start=mesh.e2n_dsp[elem], stop = start + mesh.e2n_cnt[elem];
              for(mt_int i=start; i<stop; i++)
              {
                mt_int node = mesh.e2n_con[i];
                mt_real * _msh_vel_ptr = mesh_vel.data() + 3*node;
                mt_real const * _xyz_new = mesh.xyz.data() + 3*node;
                mt_real const * _xyz_old = mesh_xyz_old.data() + 3*node;
                mt_real const * _xyz_old_old = mesh_xyz_old_old.data() + 3*node;
                if(t > 1) {
                  _msh_vel_ptr[0] = (3 * _xyz_new[0] - 4 * _xyz_old[0] + _xyz_old_old[0]) / (2 * dt * scale_factor);
                  _msh_vel_ptr[1] = (3 * _xyz_new[1] - 4 * _xyz_old[1] + _xyz_old_old[1]) / (2 * dt * scale_factor);
                  _msh_vel_ptr[2] = (3 * _xyz_new[2] - 4 * _xyz_old[2] + _xyz_old_old[2]) / (2 * dt * scale_factor);
                } else {
                  _msh_vel_ptr[0] = (_xyz_new[0] - _xyz_old[0]) / (dt * scale_factor);
                  _msh_vel_ptr[1] = (_xyz_new[1] - _xyz_old[1]) / (dt * scale_factor);
                  _msh_vel_ptr[2] = (_xyz_new[2] - _xyz_old[2]) / (dt * scale_factor);
                }
              }
            }
         }
      }
      igb_fracs[0].assign(fracs.begin(), fracs.end());
      write_igb_block(igb_fracs, frac_igb);
      if(comp_vel) {
        igb_vels[0].assign(mesh_vel.begin(), mesh_vel.end());
        write_igb_block(igb_vels, vel_igb);
      }
    }
    else
    {
      combine_outputs(fracs, all_fracs);
      igb_fracs[0].assign(fracs.begin(), fracs.end());
      write_igb_block(igb_fracs, frac_igb);
    }
  }
  progress.finish();
  fclose(frac_igb.fileptr);
  fclose(mesh_igb.fileptr);
  if(comp_vel)
    fclose(vel_igb.fileptr);
}

//calculates volumefractions for the case of fixed surfaces and a fixed mesh
//and writes them to a .dat file
void calc_fracs_fixed_mesh_fixed_surface(const mt_meshdata & mesh,
                                         const mt_vector<mt_meshdata> & surf,
                                         const vf_options opts)
{
  std::cout << "Computing volume fractions." << std::endl;
  mt_vector<mt_real> fracs;
  mt_vector<mt_vector<mt_real>> all_fracs(opts.surf.size());
  for(size_t o=0; o<opts.surf.size(); o++)
  {
    volume_fractions(mesh, surf[o], fracs);
    all_fracs[o]=fracs;
  }
  combine_outputs(fracs, all_fracs);

  std::string outfile = opts.fname_fracs.size() > 0 ? opts.fname_fracs + DAT_EXT : opts.msh_base + ".fracs" + DAT_EXT;
  write_vector_ascii(fracs, outfile);
}

//main
int main(int argc, char** argv)
{
  struct timeval t1, t2;
  struct vf_options opts;
  mt_meshdata mesh;

  int ret = vf_parse_options(argc, argv, opts);
  if (ret != 0) return 1;

  gettimeofday(&t1, NULL);

  mt_filename file(opts.msh_base, opts.ifmt);
  std::cout << "Reading mesh: " << file.base << std::endl;
  read_mesh_selected(mesh, file.format, file.base, CRP_READ_ELEM | CRP_READ_PTS);
  std::cout << "Mesh consists of " << mesh.e2n_cnt.size() << " elements, "
    << mesh.xyz.size()/3 << " nodes." << std::endl;

  std::cout << "Setting up n2e / n2n graphs .. " << std::endl;
  compute_full_mesh_connectivity(mesh, file.base);

  mt_real edge_len = avrg_edgelength_estimate(mesh);
  int dp = get_dec_power_estimate(edge_len);
  mt_real scale_factor = 1.0;
  bool scaling = false;

  if(dp < 0) {
    dp = -dp;
    scaling = true;
    scale_factor = pow(10.0, double(dp));
    std::cout << "Spatial unit of mesh is too big. ";
    std::cout << "Scaling mesh by " << scale_factor << std::endl;

    for(size_t i=0; i<mesh.xyz.size(); i++) mesh.xyz[i] *= scale_factor;
  }

  mt_vector<mt_meshdata> surf(opts.surf.size());
  for(size_t o=0; o<opts.surf.size(); o++)
    {
      mt_filename surf_file(opts.surf[o], opts.ifmt);
      std::cout << "Reading mesh: " << surf_file.base << std::endl;
      read_mesh_selected(surf[o], surf_file.format, surf_file.base, CRP_READ_ELEM | CRP_READ_PTS);
      std::cout << "Surface mesh consists of " << surf[o].e2n_cnt.size() << " elements, "
      << surf[o].xyz.size()/3 << " nodes." << std::endl;

      if(surf[o].etype[0] != Tri) {
        fprintf(stderr, "Surface mesh must be a triangle mesh! Aborting!\n");
        exit(1);
      }

      std::cout << "Setting up n2e / n2n graphs .. " << std::endl;
      compute_full_mesh_connectivity(surf[o]);

      if(scaling) {
        std::cout << "Scaling surface by " << scale_factor << std::endl;
        for(size_t i=0; i<surf[o].xyz.size(); i++) surf[o].xyz[i] *= scale_factor;
      }
    }


  if(opts.surf_xdynpt.size() > 0 && opts.mesh_xdynpt.size() > 0)
  {
    calc_fracs_moving_mesh_moving_surface(mesh, surf, opts, scale_factor);
  }
  else if(opts.surf_xdynpt.size() > 0 && opts.mesh_xdynpt.size() == 0)
  {
    calc_fracs_fixed_mesh_moving_surface(mesh, surf, opts, scale_factor);
  }
  else if(opts.surf_xdynpt.size() == 0 && opts.mesh_xdynpt.size() > 0)
  {
    calc_fracs_moving_mesh_fixed_surface(mesh, surf, opts, scale_factor);
  }
  else
  {
    calc_fracs_fixed_mesh_fixed_surface(mesh, surf, opts);
  }

  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  return 0;
}
