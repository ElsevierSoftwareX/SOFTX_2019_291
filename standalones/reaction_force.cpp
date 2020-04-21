#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <string>

#include "mt_modes_base.h"


struct reaction_force_options {

  mt_filename mesh;
  std::string surf;
  std::string stress;
  std::string out;

};

static const std::string stress_par = "-stress=";
static const std::string out_par    = "-out=";


void print_reaction_force_help()
{
  fprintf(stderr, "reaction_force: Compute reaction force on a surface.");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to the basename of the mesh.\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to the surface file.\n", surf_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to the stress force file.\n", stress_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to the force output file.\n", out_par.c_str());
  fprintf(stderr, "\n");
}

int reaction_force_parse_options(int argc, char** argv, reaction_force_options & opts)
{
  if(argc < 2) {
    print_reaction_force_help();
    return 1;
  }

  std::string msh_base;

  // parse all parameters -----------------------------------------------------------------
  for(int i=1; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh_par, msh_base);
    if(!match) match = parse_param(param, surf_par, opts.surf);
    if(!match) match = parse_param(param, stress_par, opts.stress);
    if(!match) match = parse_param(param, out_par,  opts.out);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 3;
    }
  }
  opts.mesh.assign(msh_base, "");
  fixBasename(opts.surf);
  fixBasename(opts.out);

  // check if all relevant parameters have been set ---------------------------------------------------

  if( !(opts.mesh.isSet() && opts.surf.size() && opts.stress.size() && opts.out.size()) )
  {
    std::cerr << "Error: Insufficient parameters provided." << std::endl;
    print_reaction_force_help();
    return 2;
  }

  return 0;
}

int main(int argc, char** argv)
{
  struct timeval t1, t2;
  struct reaction_force_options opts;
  struct mt_meshdata mesh;
  mt_vector<mt_int> surf_con;
  nbc_data          nbc;

  int ret = reaction_force_parse_options(argc, argv, opts);
  if (ret != 0) return 1;

  // read mesh
  std::cout << "Reading mesh: " << opts.mesh.base << std::endl;
  gettimeofday(&t1, NULL);
  read_mesh_selected(mesh, opts.mesh.format, opts.mesh.base);
  for(size_t i=0; i<mesh.xyz.size(); i++) mesh.xyz[i] *= 1e-3;
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  std::cout << "Reading surface: " << opts.surf << std::endl;
  gettimeofday(&t1, NULL);
  read_nbc(nbc, surf_con, opts.surf + NBC_EXT);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  std::vector<std::vector<double> > stress_tensor;
  igb_header igbhead;

  std::cout << "Reading stress tensor: " << opts.stress << std::endl;
  gettimeofday(&t1, NULL);
  init_igb_header(opts.stress, igbhead);
  read_igb_data(stress_tensor, igbhead, opts.stress);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  // number of read tensors must match mesh size
  assert(stress_tensor.size() > 0 && stress_tensor[0].size() == mesh.e2n_cnt.size()*9);

  std::cout << "Computing reaction force .." << std::endl;
  gettimeofday(&t1, NULL);

  std::vector<std::vector<double> > rf(stress_tensor.size());

  for(size_t tm=0; tm < stress_tensor.size(); tm++)
  {
    rf[tm].resize(nbc.eidx.size());
    for(size_t i=0; i<nbc.eidx.size(); i++)
    {
      mt_int eidx = nbc.eidx[i];
      double* s = &stress_tensor[tm][eidx*9];
      mt_point<double> p0(mesh.xyz.data() + surf_con[i*3+0]*3);
      mt_point<double> p1(mesh.xyz.data() + surf_con[i*3+1]*3);
      mt_point<double> p2(mesh.xyz.data() + surf_con[i*3+2]*3);

      mt_point<double> p01 = p1 - p0;
      mt_point<double> p02 = p2 - p0;
      mt_point<double> n   = p01.crossProd(p02);
      double len = n.length();
      double vol = len * 0.5;
      n /= len;

      mt_point<double> prod;
      prod.x = s[0]*n.x + s[1]*n.y + s[2]*n.z;
      prod.y = s[3]*n.x + s[4]*n.y + s[5]*n.z;
      prod.z = s[6]*n.x + s[7]*n.y + s[8]*n.z;

      rf[tm][i] = prod.length() * 1e-3 * vol;
    }
  }
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  std::cout << "Writing output: " << opts.out << std::endl;
  gettimeofday(&t1, NULL);

  std::string datfile = opts.out + DAT_EXT;
  std::string igbfile = opts.out + IGB_EXT;

  FILE* fd = fopen(datfile.c_str(), MT_FOPEN_WRITE);
  if(fd) {
    for(size_t tm=0; tm < rf.size(); tm++)
    {
      double time = tm * igbhead.v_inc_t;
      fprintf(fd, "%.2f ", time);

      double sum = 0.0;
      std::vector<double> seq;
      seq = rf[tm]; std::sort(seq.begin(), seq.end());
      for(auto s : seq) sum += s;
      fprintf(fd, "%.10lf\n", double(sum));
    }

    fclose(fd);
  }
  else
    treat_file_open_error(datfile, errno);


  bool write_igb = true;
  if(write_igb) {
    for(size_t tm=0; tm < rf.size(); tm++)
      for(size_t i=0; i < rf[tm].size(); i++) rf[tm][i] *= 1e6;

    igb_header wigb;
    init_igb_header(igbfile, wigb);
    set_igb_header_dim(rf[0].size(), 1, 1, rf.size(), wigb);
    set_igb_header_units("um", "ms", "micro-newton", wigb);
    set_igb_header_datatype("float", wigb);
    set_igb_header_systeme("little_endian", wigb);
    write_igb_header(wigb);
    write_igb_data(rf, wigb);
  }

  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
}





