/**
* @file deactivate_igb_based_on_tags.cpp
* @brief Zeroize IGB data in a submesh defined by a tag set
* @author Elias Karabelas
* @version
* @date 2019-12-16
*/


#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>

#include "mt_modes_base.h"
#include "igb_utils.hpp"
#include <functional> 
#include <numeric>


struct deactivate_igb_options
{
  std::string msh_base;
  std::string input_data;
  std::string output_data;
  std::string ifmt;
  std::set<int> tags;
};

void print_deactivate_igb_help()
 {
  fprintf(stderr, "deactivate_igb_tags: Zeroize IGB data in a submesh based on a tag set\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to input IGB file\n", idat_par.c_str());
  fprintf(stderr, "%stag1%ctag2\t (input) \"%c\"-seperated list of tags\n", tags_par.c_str(), tag_separator, tag_separator);
  fprintf(stderr, "%s<path>\t (output) path to output IGB file\n", odat_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh input format.\n", inp_format_par.c_str());
  fprintf(stderr, "The supported input formats are:\n%s\n", input_formats.c_str());
  fprintf(stderr, "\n");
}

int deactivate_igb_parse_options(int argc, char** argv, struct deactivate_igb_options & opts)
{
  if(argc < 3) {
    print_deactivate_igb_help();
    return 1;
  }

  // parse all enclose parameters -----------------------------------------------------------------
  for(int i=1; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh_par, opts.msh_base);
    if(!match) match = parse_param(param, idat_par, opts.input_data);
    if(!match) match = parse_param(param, odat_par, opts.output_data);
    if(!match) match = parse_param(param, inp_format_par, opts.ifmt);

    if( (!match) && param.compare(0, tags_par.size(), tags_par) == 0)
    {
      std::string tags_str;
      tags_str.assign(param.begin()+tags_par.size(), param.end());
      mt_vector<std::string> taglist;
      split_string(tags_str, tag_separator, taglist);
      for(size_t t=0; t<taglist.size(); t++)
        opts.tags.insert( atoi(taglist[t].c_str()) );
      match = true;
    }

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 3;
    }
  }
  fixBasename(opts.msh_base);


  // check if all relevant parameters have been set ---------------------------------------------------
  const bool mshok = opts.msh_base.size() > 0,
       inputok = opts.input_data.size() > 0,
       outputok = opts.output_data.size() > 0,
       tagsok = opts.tags.size() > 0;
  if( !(mshok && inputok && outputok && tagsok) )
  {
    std::cerr << "Error: Insufficient parameters provided." << std::endl;
    print_deactivate_igb_help();
    return 2;
  }

  return 0;
}

mt_real calc_lumped_l2_scalar_product(const struct mt_meshdata & mesh, const mt_vector<mt_real> & u, 
                                      const mt_vector<mt_real> & v, const int dpn)
{
  mt_real result = 0.0;
  #ifdef OPENMP
  #pragma omp parallel for schedule(dynamic, 100) reduction(+:result)
  #endif
  for(size_t eidx = 0; eidx < mesh.e2n_cnt.size(); eidx++)
  {
    const short nnodes = get_nnodes(mesh.etype[eidx]);
    const mt_int* con = mesh.e2n_con.data() + mesh.e2n_dsp[eidx];
    mt_real vol = volume(mesh.etype[eidx], con, mesh.xyz.data());
    mt_real local_result = 0.0;
    for(short k=0; k < nnodes; k++) {
      const mt_real * uptr = u.data() + con[k] * dpn;
      const mt_real * vptr = v.data() + con[k] * dpn;
      for(short j=0; j < dpn; j++) {
        local_result += uptr[j] * vptr[j];
      }
    }
    local_result *= static_cast<mt_real>(vol / nnodes);
    result += local_result;
  }
  return result;
}

void fill_rigid_body_vectors(const struct mt_meshdata & mesh, mt_vector<mt_vector<mt_real> > & rbm)
{
  rbm.resize(6);
  const size_t nnodes = mesh.xyz.size() / 3;
  for(size_t i=0; i < rbm.size(); i++) {
    rbm[i].resize(mesh.xyz.size());
  }  
  
  for(size_t k=0; k < nnodes; k++) {
    mt_real * rbm0 = rbm[0].data() + k * 3;
    mt_real * rbm1 = rbm[1].data() + k * 3;
    mt_real * rbm2 = rbm[2].data() + k * 3;
    mt_real * rbm3 = rbm[3].data() + k * 3;
    mt_real * rbm4 = rbm[4].data() + k * 3;
    mt_real * rbm5 = rbm[5].data() + k * 3;
    const mt_real * xyz  = mesh.xyz.data() + k * 3;
    rbm0[0] = 1.0; rbm0[1] = 0.0; rbm0[2] = 0.0;
    rbm1[0] = 0.0; rbm1[1] = 1.0; rbm1[2] = 0.0;
    rbm2[0] = 0.0; rbm2[1] = 0.0; rbm2[2] = 1.0;
    rbm3[0] = xyz[1]; rbm3[1] = -xyz[0]; rbm3[2] = 0.0;
    rbm4[0] = 0.0; rbm4[1] = -xyz[2]; rbm4[2] = xyz[1];
    rbm5[0] = xyz[2]; rbm5[1] = 0.0; rbm5[2] = -xyz[0];
  }

  //orthogonalize the vectors
  mt_real norm = std::sqrt(calc_lumped_l2_scalar_product(mesh, rbm[0], rbm[0], 3));
  std::transform(rbm[0].begin(), rbm[0].end(), rbm[0].begin(),
                 std::bind(std::multiplies<mt_real>(), std::placeholders::_1, 1. / norm));
  
  for(size_t i=1; i < rbm.size(); i++) {
    mt_real dot = 0.0;
    for(size_t j=0; j < i; j++) {
      dot = calc_lumped_l2_scalar_product(mesh, rbm[i], rbm[j], 3);
      for(size_t k=0; k < rbm[j].size(); k++)
        rbm[i][k] -= dot * rbm[j][k];
    }
    norm = std::sqrt(calc_lumped_l2_scalar_product(mesh, rbm[i], rbm[i],3));
    std::transform(rbm[i].begin(), rbm[i].end(), rbm[i].begin(),
                   std::bind(std::multiplies<mt_real>(), std::placeholders::_1, 1. / norm));
  }
}

int main(int argc, char** argv)
{
  struct timeval t1, t2;
  struct deactivate_igb_options opts;
  igb_header igb_head_from;
  igb_header igb_head_to;
  int ret = deactivate_igb_parse_options(argc, argv, opts);
  if(ret != 0) return 1;
  struct mt_meshdata mesh;
  mt_filename mshfile(opts.msh_base, opts.ifmt);
  std::cout << "Reading mesh: " << mshfile.base << std::endl;
  gettimeofday(&t1, NULL);
  read_mesh_selected(mesh, mshfile.format, mshfile.base);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
  // we need full mesh connectivity as we might need to recover surface info
  std::cout << "Setting up n2e / n2n graphs .. " << std::endl;
  compute_full_mesh_connectivity(mesh, mshfile.base);
  
  //parse the igb header
  //initialize the header
  init_igb_header(opts.input_data, igb_head_from);
  read_igb_header(igb_head_from);
  //first insert the data into the reference mesh
  const mt_int dpn = static_cast<mt_int>(Num_Components[igb_head_from.v_type]);

  //set the new header
  igb_head_to          = igb_head_from; // does not copy the fileptr
  igb_head_to.filename = opts.output_data;
  igb_head_to.fileptr  = NULL;
  write_igb_header(igb_head_to);

  std::stringstream out;
  out << "Nullifying " << opts.input_data << " in tags: ";
  for(std::set<int>::iterator it = opts.tags.begin(); it != opts.tags.end(); ++it)
    out << *it << " ";
  out << " and output to " << opts.output_data;
  
  gettimeofday(&t1, NULL);
  const mt_int tsteps = igb_head_from.v_t;
  const mt_int bsize  = 1;
  const mt_int nnodes = static_cast<mt_int>(mesh.xyz.size() / 3);
  const bool nodebased = (nnodes == static_cast<mt_int>(igb_head_from.v_x));
  const mt_int datasize = nodebased ? nnodes : static_cast<mt_int>(mesh.etags.size());
  
  std::vector<std::vector<float> > igb_data_input(bsize);
  std::vector<std::vector<float> > igb_data_output;
  igb_data_output.assign(bsize, std::vector<float>(dpn * datasize, 0.0f));

  mt_vector<mt_vector<mt_real> > rbm;
  fill_rigid_body_vectors(mesh, rbm);
  mt_vector<mt_real> rbmnorms(6);
  mt_vector<mt_real> rbmweights(6);
  for(size_t i=0; i < rbmnorms.size(); i++)
    rbmnorms[i] = calc_lumped_l2_scalar_product(mesh, rbm[i], rbm[i], 3);

  PROGRESS<mt_int> progress(tsteps, out.str().c_str());
  for(mt_int i=0; i < tsteps; i++)
  {
    progress.next();
    read_igb_block(igb_data_input, bsize, igb_head_from);
    mt_vector<mt_real> disp;
    disp.assign(igb_data_input[0].begin(), igb_data_input[0].end());
    mt_vector<mt_real> disp_ortho;
    disp_ortho.assign(disp.size(), 0.0);
    for(size_t j=0; j < rbmweights.size(); j++) {
      rbmweights[j] = calc_lumped_l2_scalar_product(mesh, disp, rbm[j], dpn);
      std::cout << "RBMWEIGHT = " << rbmweights[j] << " NORM = " << rbmnorms[j] << std::endl;
      for(size_t k = 0; k < disp_ortho.size(); k++)
        disp_ortho[k] += (rbmweights[j] / rbmnorms[j]) * rbm[j][k];
    }
    // Loop over elements
    for(size_t eidx=0; eidx < mesh.etags.size(); eidx++) {
      if(nodebased) {
        const mt_int start = mesh.e2n_dsp[eidx], stop = start + mesh.e2n_cnt[eidx];
        for(mt_int j = start; j < stop; j++) {
          mt_int nidx = mesh.e2n_con[j];
          if(!opts.tags.count(mesh.etags[eidx])) {
            std::copy(igb_data_input[0].data() + nidx*dpn, igb_data_input[0].data() + nidx*dpn + dpn,
                      igb_data_output[0].data() + nidx * dpn);
          } else {
            std::copy(disp_ortho.data() + nidx*dpn, disp_ortho.data() + nidx*dpn + dpn,
                      igb_data_output[0].data() + nidx * dpn);
          }
        }
      } else {
        if(!opts.tags.count(mesh.etags[eidx]))
          std::copy(igb_data_input[0].data() + eidx*dpn, igb_data_input[0].data() + eidx*dpn + dpn,
                    igb_data_output[0].data() + eidx * dpn);
      }
    }
    write_igb_block(igb_data_output, igb_head_to);
  }
  progress.finish();

  fclose(igb_head_from.fileptr);
  fclose(igb_head_to.fileptr);

  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
  return 0;
}
