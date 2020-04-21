/**
* @file tara_msh_pp.cpp
* @brief Comparison/correction of tags between image stack and Tarantula mesh.
* @author Anton Prassl
* @version
* @date 2019-06-11
*/

#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "mt_modes_base.h"
#include "itk_utils.h"

struct retag_options{
  std::string msh_base;
  std::string img_base;
  std::string lut;
  std::string outmsh_base;
  std::string ifmt;
  std::string ofmt;
};

static const std::string img_par = "-img=";
static const std::string lut_par = "-lut=";


void print_retag_help()
{
  fprintf(stderr, "tara_msh_retag: identify and re-tag elements that were not correctly mapped during mesh generation.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the image stack.\n", img_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh.\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) tags as LUT between image stack and mesh.\n", lut_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh.\n", mesh_par.c_str());
  fprintf(stderr, "%s<format> (output) path to the basename of the output mesh: \n", outmesh_par.c_str());
  fprintf(stderr, "%s<path>\t (optional) mesh input format.\n", inp_format_par.c_str());
  fprintf(stderr, "The format of the LUT is:\n");
  fprintf(stderr, "itagA1,itagA2,..[:]etagB1,etagB2,..\n");
  fprintf(stderr, "Image and element tag regions are separated by \",\".\n"
                  "Image tags and element tags are opposed by \":\".\n"
                  "The first image tag relates to the first element tag.\n"
                  "The number of image tags must match the number of element tags.\n\n");
  fprintf(stderr, "\n");
}

int retag_parse_options(int argc, char** argv, struct retag_options & opts)
{
  if(argc < 2) {
    print_retag_help();
    return 1;
  }

  // parse all retag parameters -----------------------------------------------------------------
  for(int i=1; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, img_par, opts.img_base);
    if(!match) match = parse_param(param, mesh_par, opts.msh_base);
    if(!match) match = parse_param(param, lut_par, opts.lut);
    if(!match) match = parse_param(param, outmesh_par, opts.outmsh_base);
    if(!match) match = parse_param(param, inp_format_par, opts.ifmt);
    if(!match) match = parse_param(param, out_format_par, opts.ofmt);
    
    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 3;
    }
  }
  fixBasename(opts.msh_base);
  fixBasename(opts.outmsh_base);

  // check if all relevant parameters have been set ---------------------------------------------------
  bool mshok = opts.msh_base.size() > 0, outmshok = opts.outmsh_base.size() > 0;

  if( !(mshok && outmshok) )
  {
    std::cerr << "Error: Insufficient parameters provided." << std::endl;
    print_retag_help();
    return 2;
  }

  return 0;
}

/**
* @brief Extract the lookup table translating image labels to mesh labels.
*
* @param [in]  operstr  The string defining the operations.
* @param [out] lut      The row index corresponds with the image label.
*                       The first value with the incorrect mesh label.
*                                                          
*/
void extract_lut(const std::string & lutstr,
                 mt_vector<mt_int> & lut)
{
  mt_vector<std::string> taskslist,task, ietagstr, etagstr;
  
  split_string(lutstr, ';', taskslist);
  
  for(size_t i=0; i<taskslist.size(); i++) {
    
    split_string(taskslist[i], ':', task); // extract ith set of instructions
    split_string(task[0], ',', ietagstr);  // two entries: image tag + incorrect element tag
    split_string(task[1], ',', etagstr);   // one number:  supposed correct element tag
    
    if(ietagstr.size() != 2) {
      std::cerr << "tara_msh_retag error: set of instructions does not "
                   "have the structure 'itag,wrong_etag:corr_etag'" << std::endl;
      print_retag_help();
      exit(4);
    }

    size_t itag   = atoi(ietagstr[0].c_str());
    size_t etag   = atoi(ietagstr[1].c_str());
    size_t etag_c = atoi(etagstr[0].c_str());

    size_t rows = lut.size() /2;
    if (rows <= itag) {
        lut.resize((itag+1)*2, 0); // +1 because of the zero
    }
    lut[itag*2+0] = etag;
    lut[itag*2+1] = etag_c;
  }
}



int main(int argc, char** argv)
{
  struct timeval       t1, t2;
  struct retag_options opts;
  struct mt_meshdata   mesh;
  mt_vector<mt_int>    lut;
  itk_image            img;
  const int dtype = 1; // enforce data type "unsigned char" for any image stack
  mt_int updates  = 0; // pixel label change counter
  
  int ret = retag_parse_options(argc, argv, opts);
  if (ret != 0) return 1;

  std::cout << "Reading mesh: " << opts.msh_base << std::endl;
  gettimeofday(&t1, NULL);
  read_mesh_selected(mesh, opts.ifmt, opts.msh_base);
  compute_full_mesh_connectivity(mesh);
  
  std::cout << "Extracting lookup table: " << std::endl;
  extract_lut(opts.lut, lut);

#ifdef MT_DEBUG
  // TODO: this define does not work - don't know why
  for (size_t i=0; i<lut.size()/2; i+=2)
    std::cout << "lut[" << i/2 << "]=" << lut[i/2+0] << " " <<lut[i/2+1] <<std::endl;
#endif

  std::cout << "Reading image stack: " << opts.img_base << std::endl;
  img.read_file(opts.img_base.c_str());
  img.convert(dtype);


  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;


  /* Loop over all elements
     Ignore element tags not listed in the LUT
     Change current element tag to match the image tag if necessary
  */

  mt_vector<mt_int> tags = mesh.etags;
  itk_access<unsigned char> imdata(img);
  
  PROGRESS<size_t> prg(tags.size(), "Identifying mapping mismatches: ");
  
  gettimeofday(&t1, NULL);
  #ifdef OPENMP
  #pragma omp parallel for schedule(guided, 10)
  #endif
  for(size_t eidx = 0; eidx < mesh.e2n_cnt.size(); eidx++)
  {
    // identify current element center
    vec3r ectr = element_centerpoint(mesh, eidx);
    // retrieve corresponding ijk position based on element center position
    vec3i ijk = imdata.pixel<float>(ectr);

    // Note: It can happen that the reference mesh is larger than the image stack due
    //       to e.g. data reduction reasons.
    // assert search to bbox of image stack
    if (ijk.x < 0 || ijk.x >= imdata.dim.v1  ||
        ijk.y < 0 || ijk.y >= imdata.dim.v2  ||
        ijk.z < 0 || ijk.z >= imdata.dim.v3) {
      //std::cout << "Element " << eidx << " is not located inside image bbx" << std::endl;
      continue;
    }

    // query pixel value
    unsigned char _pxlabel = imdata(ijk.x, ijk.y, ijk.z);
    mt_int         pxlabel = mt_int(_pxlabel);

    // don't care about image labels the user did not specify
    if(!lut[pxlabel*2+0])
      continue;
    if(pxlabel >= static_cast<mt_int>((lut.size()/2+0)))
      continue;
    
    // find out about mesh element tag
    mt_int ctag = tags[eidx];

    // decide wether to change the element tag or not
    if (ctag == lut[pxlabel*2+0]) {
      //imdata(ijk.x, ijk.y, ijk.z) = lut[pxlabel];
      mesh.etags[eidx] = lut[pxlabel*2+1];
      updates++;
    }

#ifdef MT_DEBUG
    // TODO: this define does not work - don't know why
    std::cout << "Pixel value " << i << " " << j << " " << k << " " << "= " << +px;
    std::cout << " " << pos.x << "," << pos.y << "," << pos.z << std::endl << std::flush;
#endif
    #ifdef OPENMP
    #pragma omp critical
    #endif
    prg.next();
  }
  prg.finish();

  std::cout << "Changed " << updates << " element tags." << std::endl;
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;


  std::cout << "Writing mesh: " << opts.outmsh_base << std::endl;
  gettimeofday(&t1, NULL);
  write_mesh_selected(mesh, opts.ofmt, opts.outmsh_base);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;


  return 0;
}
