/**
* @file itk_mode.cpp
* @brief Manipulation of ITK structured data.
* @author Aurel Neic, Matthias Gsell
* @version
* @date 2017-02-16
*/

#include "mt_modes_base.h"
#include "itk_utils.h"
#include "kdtree.h"

/// The operation to perform
enum ITK_OP {ITK_SMOOTH, ITK_CLOSE, ITK_NORMALIZE, ITK_PADDING, ITK_CROP, ITK_FLIP,
             ITK_DTYPE, ITK_REFINE, ITK_SAMPLE, ITK_EXTRACT};

static const std::string ref_par = "-ref=";
static const std::string reclamp_par = "-reclamp=";
static const std::string axes_par = "-axes=";
static const std::string dtype_par = "-dtype=";

#define ITK_PADDING_SIZE_DEFAULT "1,1,1"
#define ITK_AXES_DEFAULT "xyz"
#define ITK_DTYPE_DEFAULT 4

/// The itk mode options struct
struct itk_options
{
  std::string msh;
  std::string outmsh;
  std::string tags;
  std::string smth;
  std::string iter;
  std::string thr;
  std::string ref;
  std::string reclamp;
  std::string size;
  std::string axes;
  std::string dtype;
  std::string min;
  std::string max;
  std::string surf;
  std::string idx;
  ITK_OP op;
};

/**
* @brief itk smooth help message.
*/
void print_itk_smooth_help()
{
  fprintf(stderr, "itk smooth: Smooth the voxel data.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t\t (input) Path of the mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<int>[,<int>,<int>]\t (optional) Discrete radius used for smoothing, single int or int triple (default 1).\n", thr_par.c_str());
  fprintf(stderr, "%s<int>[,<int>,<int>]\t (optional) Pixel refinement factor, single int or int triple (default 1).\n", ref_par.c_str());
  fprintf(stderr, "%s<int>\t\t\t (optional) Number of smoothing iterations per clamping step.\n"
                  "\t\t\t\t    0 means only clamp after smoothing. Default is 0.\n", reclamp_par.c_str());
  fprintf(stderr, "%s<int>\t\t\t (optional) Number of smoothing iter (default %d).\n", iter_par.c_str(), SMOOTH_ITER_DEFAULT);
  fprintf(stderr, "%s<float>\t\t\t (optional) Smoothing coefficient (default %.2f).\n", smooth_par.c_str(), SMOOTH_DEFAULT);
  fprintf(stderr, "%s<tag1,tag2/tag1..>\t (optional) List of tags.\n", tags_par.c_str());
  fprintf(stderr, "%s<path>\t\t\t (output) Path of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "\n");
}

/**
* @brief itk smooth help message.
*/
void print_itk_close_help()
{
  fprintf(stderr, "itk close: Apply closing (i.e. dilate-erode) algorithm to itk data.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t\t (input) Path of the mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<tag1,tag2,..>\t\t (input) List of tags. Each tag is processed individually.\n", tags_par.c_str());
  fprintf(stderr, "%s<int>,[<int>,<int>]\t (optional) Discrete radius used for smoothing, single int or int triple (default 5).\n", thr_par.c_str());
  fprintf(stderr, "%s<int>[,<int>,<int>]\t (optional) Pixel refinement factor, single int or int triple (default 1).\n", ref_par.c_str());
  fprintf(stderr, "%s<int>\t\t\t (optional) Number of iterations. Default is 1.\n", iter_par.c_str());
  fprintf(stderr, "%s<path>\t\t\t (output) Path of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "\n");
}
/**
 * @brief itk normalize help
 */
void print_itk_normalize_help()
{
  fprintf(stderr, "itk normalize: Normalize voxel spacing.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t\t (input) Path of the mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t\t\t (output) Path of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "%s<int>[,<int>,<int>]\t (optional) Pixel refinement factor, single int or int triple (default 1).\n", ref_par.c_str());
  fprintf(stderr, "\n");
}


/**
 * @brief itk padding help
 */
void print_itk_padding_help()
{
  fprintf(stderr, "itk padding: add padding to voxel data.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t\t (input) Path of the mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t\t\t (output) Path of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "%s<int>,<int>,<int>\t\t (input) Comma separated padding sizes for each axis (default '%s').\n", size_par.c_str(), ITK_PADDING_SIZE_DEFAULT);
  fprintf(stderr, "\n");
}

/**
 * @brief itk crop help
 */
void print_itk_crop_help()
{
  fprintf(stderr, "itk crop: remove surrounding whitespace.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t\t (input) Path of the mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t\t\t (output) Path of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "%s<int>[,<int>,<int>]\t (optional) Pixel refinement factor, single int or int triple (default 1).\n", ref_par.c_str());
  fprintf(stderr, "\n");
}

/**
 * @brief itk flip help
 */
void print_itk_flip_help()
{
  fprintf(stderr, "itk flip: flip the voxel data along given axes.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t\t (input) Path of the mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t\t\t (output) Path of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "%s<str>\t\t\t (input) Axes along which the image is flipped (default '%s')\n", axes_par.c_str(), ITK_AXES_DEFAULT);
  fprintf(stderr, "\t\t\t\t    the string consists of the following characters {'x', 'y', 'z'}\n");
  fprintf(stderr, "\n");
}

/**
 * @brief itk dtype help
 */
void print_itk_dtype_help()
{
  fprintf(stderr, "itk dtype: convert datatype.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t\t (input) Path of the mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t\t\t (output) Path of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "%s<int>\t\t\t (input) New datatype (default '%d')\n", dtype_par.c_str(), ITK_DTYPE_DEFAULT);
  for (int i = 1; i < 11; i++)
    fprintf(stderr, "\t\t\t\t    %d  %s\n", i, itk_get_datatype_str(i));
  fprintf(stderr, "\t\t\t\t    %d  color scalars\n", 11);
  fprintf(stderr, "\n");
}

/**
 * @brief itk dtype help
 */
void print_itk_refine_help()
{
  fprintf(stderr, "itk refine: refine voxel data.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t\t (input) Path of the mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t\t\t (output) Path of the output mesh\n", outmesh_par.c_str());
  fprintf(stderr, "%s<int>[,<int>,<int>]\t (optional) Pixel refinement factor, single int or int triple (default 1).\n", ref_par.c_str());
  fprintf(stderr, "\n");
}

/**
 * @brief itk new help
 */
void print_itk_sample_help()
{
  fprintf(stderr, "itk sample: create an itk image stack from sampeling surfaces.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t\t (input) Path of the image stack to sample on\n", mesh_par.c_str());
  fprintf(stderr, "%s<x,y,z>\t\t\t (input) Lower left point of ITK image space\n", min_par.c_str());
  fprintf(stderr, "%s<x,y,z>\t\t\t (input) Upper right point of ITK image space\n", max_par.c_str());
  fprintf(stderr, "%s<size>\t\t\t (input) ITK pixel spacing\n", size_par.c_str());
  fprintf(stderr, "%s<surf1,surf2>\t\t (input) List of surfaces to sample\n", surf_par.c_str());
  fprintf(stderr, "%s<path>\t\t\t (output) Path of the output file\n", outmesh_par.c_str());
  fprintf(stderr, "%s<string>\t\t\t (optional) ITK datatype (default '%s')\n", dtype_par.c_str(), itk_get_datatype_str(ITK_DTYPE_DEFAULT));
  fprintf(stderr, "\n");
  fprintf(stderr, "Note: Either use (%s) or (%s, %s, %s) input parameters to define the image stack\n"
                  "to sample on.\n", mesh_par.c_str(), min_par.c_str(), max_par.c_str(), size_par.c_str());
  fprintf(stderr, "\n");
}

/**
 * @brief itk new help
 */
void print_itk_extract_help()
{
  fprintf(stderr, "itk extract: extract slices of an itk image stack.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t\t (input) Path of the image stack to sample on\n", mesh_par.c_str());
  fprintf(stderr, "%s<plane>\t\t\t (input) Plane-axes from which to extract (xy, xz, yz)\n", axes_par.c_str());
  fprintf(stderr, "%s<int,>\t\t\t (input) Comma separated indices of slices to extract\n", idx_par.c_str());
  fprintf(stderr, "%s<path>\t\t\t (output) Path of the output file\n", outmesh_par.c_str());
  fprintf(stderr, "\n");
}
/**
* @brief itk mode options parser.
*
* @param [in]  argc Arguments count.
* @param [in]  argv Arguments string-array.
* @param [out] opts Options structure.
*
* @return
*/
int itk_parse_options(int argc, char** argv, struct itk_options & opts)
{
  if(argc < 3) {
    fprintf(stderr, "Please choose one of the following modes:\n");
    print_itk_smooth_help();
    print_itk_close_help();
    print_itk_normalize_help();
    print_itk_padding_help();
    print_itk_crop_help();
    print_itk_flip_help();
    print_itk_dtype_help();
    print_itk_refine_help();
    print_itk_sample_help();
    print_itk_extract_help();
    return 1;
  }

  std::string itkmode = argv[2];

  // parse parameters -----------------------------------------------------------------
  for(int i=3; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh_par, opts.msh);
    if(!match) match = parse_param(param, outmesh_par, opts.outmsh);
    if(!match) match = parse_param(param, tags_par, opts.tags);
    if(!match) match = parse_param(param, smooth_par, opts.smth);
    if(!match) match = parse_param(param, iter_par, opts.iter);
    if(!match) match = parse_param(param, thr_par, opts.thr);
    if(!match) match = parse_param(param, ref_par, opts.ref);
    if(!match) match = parse_param(param, reclamp_par, opts.reclamp);
    if(!match) match = parse_param(param, size_par, opts.size);
    if(!match) match = parse_param(param, axes_par, opts.axes);
    if(!match) match = parse_param(param, dtype_par, opts.dtype);
    if(!match) match = parse_param(param, min_par, opts.min);
    if(!match) match = parse_param(param, max_par, opts.max);
    if(!match) match = parse_param(param, surf_par, opts.surf);
    if(!match) match = parse_param(param, idx_par, opts.idx);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 2;
    }
  }

  if(itkmode.compare("smooth") == 0)
  {
    opts.op = ITK_SMOOTH;

    // check if all relevant parameters have been set ---------------------------------------------------
    if( ! (opts.msh.size() > 0 && opts.outmsh.size() > 0 ) )
    {
      std::cerr << "itk smooth: Insufficient parameters provided." << std::endl;
      print_itk_smooth_help();
      return 4;
    }
  }
  else if (itkmode.compare("close") == 0)
  {
    opts.op = ITK_CLOSE;

    // check if all relevant parameters have been set --------------------------------------------------
    if( ! (opts.msh.size() > 0 && opts.tags.size() && opts.outmsh.size() > 0 ) )
    {
      std::cerr << "itk close: Insufficient parameters provided." << std::endl;
      print_itk_close_help();
      return 4;
    }
  }
  else if (itkmode.compare("normalize") == 0)
  {
    opts.op = ITK_NORMALIZE;

    // check if all relevant parameters have been set --------------------------------------------------
    if( ! (opts.msh.size() > 0 && opts.outmsh.size() > 0 ) )
    {
      std::cerr << "itk normalize: Insufficient parameters provided." << std::endl;
      print_itk_normalize_help();
      return 4;
    }
  }
  else if (itkmode.compare("padding") == 0)
  {
    opts.op = ITK_PADDING;

    // check if all relevant parameters have been set --------------------------------------------------
    if( ! (opts.msh.size() > 0 && opts.outmsh.size() > 0 ) )
    {
      std::cerr << "itk padding: Insufficient parameters provided." << std::endl;
      print_itk_padding_help();
      return 4;
    }
  }
  else if (itkmode.compare("crop") == 0)
  {
    opts.op = ITK_CROP;

    // check if all relevant parameters have been set --------------------------------------------------
    if( ! (opts.msh.size() > 0 && opts.outmsh.size() > 0 ) )
    {
      std::cerr << "itk crop: Insufficient parameters provided." << std::endl;
      print_itk_crop_help();
      return 4;
    }
  }
  else if (itkmode.compare("flip") == 0)
  {
    opts.op = ITK_FLIP;

    // check if all relevant parameters have been set -------------------------------
    if( ! (opts.msh.size() > 0 && opts.outmsh.size() > 0 ) )
    {
      std::cerr << "itk flip: Insufficient parameters provided." << std::endl;
      print_itk_flip_help();
      return 4;
    }
  }
  else if (itkmode.compare("dtype") == 0)
  {
    opts.op = ITK_DTYPE;

    // check if all relevant parameters have been set -------------------------------
    if( ! (opts.msh.size() > 0 && opts.outmsh.size() > 0 ) )
    {
      std::cerr << "itk dtype: Insufficient parameters provided." << std::endl;
      print_itk_dtype_help();
      return 4;
    }
  }
  else if (itkmode.compare("refine") == 0)
  {
    opts.op = ITK_REFINE;

    // check if all relevant parameters have been set -------------------------------
    if( ! (opts.msh.size() > 0 && opts.outmsh.size() > 0 ) )
    {
      std::cerr << "itk dtype: Insufficient parameters provided." << std::endl;
      print_itk_refine_help();
      return 4;
    }
  }
  else if (itkmode.compare("sample") == 0)
  {
    opts.op = ITK_SAMPLE;

    bool img_stack_defined = opts.msh.size() ||
         (opts.min.size() && opts.max.size() && opts.size.size());

    // check if all relevant parameters have been set -------------------------------
    if( !(opts.surf.size() && img_stack_defined && opts.outmsh.size()) )
    {
      std::cerr << "itk sample: Insufficient parameters provided." << std::endl;
      print_itk_sample_help();
      return 4;
    }
  }
  else if (itkmode.compare("extract") == 0)
  {
    opts.op = ITK_EXTRACT;

    bool plane_defined = opts.axes.size() && ((opts.axes == "xy") || (opts.axes == "yz") || (opts.axes == "yz"));

    // check if all relevant parameters have been set -------------------------------
    if( !(opts.msh.size() && opts.idx.size() && plane_defined && opts.outmsh.size()) )
    {
      std::cerr << "itk extract: Insufficient parameters provided." << std::endl;
      print_itk_extract_help();
      return 4;
    }
  }
  else {
    print_usage(argv[0]);
    return 2;
  }
  return 0;
}



void itk_sample(itk_image & img, const mt_meshdata & surfmesh, const kdtree & tree)
{
  itk_access<short> data(img);
  // we assume that the surface has a uniform tag index and require that it is nonzero
  short curtag = surfmesh.etags[0];
  curtag = curtag != 0 ? curtag : 1;

  PROGRESS<unsigned int> prg(data.dim.v1, "Sampling progress: ");

  for(unsigned int i=0; i<data.dim.v1; i++) {
    #ifdef OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for(unsigned int j=0; j<data.dim.v2; j++)
      for(unsigned int k=0; k<data.dim.v3; k++) {
        vec3f pos = data.position<float>(i,j,k);
        if(inside_closed_surface(tree, pos))
          data(i,j,k) = curtag;
      }

    prg.next();
  }

  prg.finish();
}



/**
* @brief Smoothing mode function.
*
* @param [in]  argc Arguments count.
* @param [in]  argv Arguments string-array.
*/
void itk_mode(int argc, char** argv)
{
  struct itk_options opts;
  int ret = itk_parse_options(argc, argv, opts);
  struct timeval t1, t2;

  if(ret != 0) return;

  switch(opts.op)
  {
    case ITK_SMOOTH:
    {
      double smooth  = opts.smth.size()    > 0 ? atof(opts.smth.c_str())    : SMOOTH_DEFAULT;
      int    iter    = opts.iter.size()    > 0 ? atoi(opts.iter.c_str())    : SMOOTH_ITER_DEFAULT;
      int    reclamp = opts.reclamp.size() > 0 ? atoi(opts.reclamp.c_str()) : 0;

      triple<short> rad = {1, 1, 1};
      if (opts.thr.size() > 0)
      {
        mt_vector<std::string> radlist;
        split_string(opts.thr.c_str(), ',', radlist);
        if (radlist.size() == 1)
        {      
          rad.v1 = rad.v2 = rad.v3 = atoi(radlist[0].c_str());
        }
        else if (radlist.size() == 3)
        {
          rad.v1 = atoi(radlist[0].c_str());
          rad.v2 = atoi(radlist[1].c_str());
          rad.v3 = atoi(radlist[2].c_str());
        }
        else 
        {
          std::cerr << std::endl << "argument error: '" << thr_par.c_str() << "' expects one single int or a int triple!" << std::endl;
          print_itk_smooth_help();
          return;
        }
      }

      triple<short> ref = {1, 1, 1};
      if (opts.ref.size() > 0)
      { 
        mt_vector<std::string> reflist;
        split_string(opts.ref.c_str(), ',', reflist);
        if (reflist.size() == 1)
        {      
          ref.v1 = ref.v2 = ref.v3 = atoi(reflist[0].c_str());
        }
        else if (reflist.size() == 3)
        {
          ref.v1 = atoi(reflist[0].c_str());
          ref.v2 = atoi(reflist[1].c_str());
          ref.v3 = atoi(reflist[2].c_str());
        }
        else 
        {
          std::cerr << std::endl << "argument error: '" << ref_par.c_str() << "' expects one single int or a int triple!" << std::endl;
          print_itk_smooth_help();
          return;
        }
      }

      itk_image input_img, comp_img;

      std::cout << "Reading image: " << opts.msh << std::endl;
      gettimeofday(&t1, NULL);
      input_img.read_file(opts.msh.c_str());
      if (((ref.v1 > 0) && (ref.v2 > 0) && (ref.v3 > 0)) &&
          ((ref.v1 > 1) || (ref.v2 > 1) || (ref.v3 > 1)))      
      input_img.refine(ref);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Converting image data to float .. " << std::endl;
      gettimeofday(&t1, NULL);
      comp_img.assign(input_img, "float");
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      if(opts.tags.size()) {
        mt_vector<std::string> operlist;
        split_string(opts.tags, '/' , operlist);

        itk_image timg;
        timg.assign(comp_img.dim, comp_img.orig, comp_img.pixspc, comp_img.ncomp, "float");

        for(size_t oper = 0; oper < operlist.size(); oper++)
        {
          mt_vector<std::string> taglist;
          std::set<float> tg;

          split_string(operlist[oper], ',' , taglist);
          for(std::string & s : taglist) tg.insert(atof(s.c_str()));

          timg.databuff.zero();
          itk_extract_vals(comp_img, timg, tg);

          printf("Smoothing tags %s\n", operlist[oper].c_str());
          gettimeofday(&t1, NULL);

          smooth_itk(timg, rad, smooth, iter);

          gettimeofday(&t2, NULL);
          std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

          std::cout << "Clamping .. " << std::endl;
          gettimeofday(&t1, NULL);
          clamp_itk(timg, input_img, rad);
          gettimeofday(&t2, NULL);
          std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

          itk_insert_vals(timg, comp_img, tg);
        }
      }
      else if (reclamp == 0) {
        std::cout << "Smoothing .. " << std::endl;
        gettimeofday(&t1, NULL);
        smooth_itk(comp_img, rad, smooth, iter);
        // comp_img.write_file("debug.smth.vtk");
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

        std::cout << "Clamping .. " << std::endl;
        gettimeofday(&t1, NULL);
        clamp_itk(comp_img, input_img, rad);
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      }
      else {
        std::cout << "Smoothing .. " << std::endl;
        gettimeofday(&t1, NULL);
        int it_idx = iter % reclamp;

        smooth_itk(comp_img, rad, smooth, it_idx);
        clamp_itk(comp_img, input_img, rad);

        while(it_idx < iter)
        {
          smooth_itk(comp_img, rad, smooth, reclamp);
          clamp_itk(comp_img, input_img, rad);
          it_idx += reclamp;
        }
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      }

      std::cout << "Writing output .. " << std::endl;
      gettimeofday(&t1, NULL);
      input_img.assign(comp_img, itk_get_datatype_str(input_img.comptypeid));
      input_img.write_file(opts.outmsh.c_str());
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      break;
    }

    case ITK_CLOSE:
    {
      int    iter    = opts.iter.size()    > 0 ? atoi(opts.iter.c_str())    : 1;

      triple<short> rad = {1, 1, 1};
      if (opts.thr.size() > 0)
      {      
        mt_vector<std::string> radlist;
        split_string(opts.thr.c_str(), ',', radlist);
        if (radlist.size() == 1)
        {      
          rad.v1 = rad.v2 = rad.v3 = atoi(radlist[0].c_str());
        }
        else if (radlist.size() == 3)
        {
          rad.v1 = atoi(radlist[0].c_str());
          rad.v2 = atoi(radlist[1].c_str());
          rad.v3 = atoi(radlist[2].c_str());
        }
        else 
        {
          std::cerr << std::endl << "argument error: '" << thr_par.c_str() << "' expects one single int or a int triple!" << std::endl;
          print_itk_close_help();
          return;
        }
      }
  
      triple<short> ref = {1, 1, 1};
      if (opts.ref.size() > 0)
      {
        mt_vector<std::string> reflist;
        split_string(opts.ref.c_str(), ',', reflist);
        if (reflist.size() == 1)
        {      
          ref.v1 = ref.v2 = ref.v3 = atoi(reflist[0].c_str());
        }
        else if (reflist.size() == 3)
        {
          ref.v1 = atoi(reflist[0].c_str());
          ref.v2 = atoi(reflist[1].c_str());
          ref.v3 = atoi(reflist[2].c_str());
        }
        else 
        {
          std::cerr << std::endl << "argument error: '" << ref_par.c_str() << "' expects one single int or a int triple!" << std::endl;
          print_itk_close_help();
          return;
        }
      }

      itk_image input_img, comp_img;

      std::cout << "Reading image: " << opts.msh << std::endl;
      gettimeofday(&t1, NULL);
      input_img.read_file(opts.msh.c_str());      
      if (((ref.v1 > 0) && (ref.v2 > 0) && (ref.v3 > 0)) &&
          ((ref.v1 > 1) || (ref.v2 > 1) || (ref.v3 > 1)))      
      input_img.refine(ref);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Converting image data to float .. " << std::endl;
      gettimeofday(&t1, NULL);
      comp_img.assign(input_img, "float");
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      mt_vector<std::string> taglist;
      split_string(opts.tags, ',' , taglist);

      itk_image timg;
      timg.assign(comp_img.dim, comp_img.orig, comp_img.pixspc, comp_img.ncomp, "float");

      for(size_t tidx = 0; tidx < taglist.size(); tidx++)
      {
        std::set<float> tg;
        tg.insert(atof(taglist[tidx].c_str()));

        timg.databuff.zero();
        itk_extract_vals(comp_img, timg, tg);

        printf("Closing tag %.2f\n", *tg.begin());
        gettimeofday(&t1, NULL);

        for(int it=0; it < iter; it++)
        {
          itk_dilate(timg, rad, *tg.begin());
          itk_erode(timg, rad);
        }
        itk_insert_vals(timg, comp_img, tg);

        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      }

      std::cout << "Writing output .. " << std::endl;
      gettimeofday(&t1, NULL);
      input_img.assign(comp_img, "short");
      input_img.write_file(opts.outmsh.c_str());
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      break;
    }

    case ITK_NORMALIZE:
    {
      mt_vector<std::string> reflist;
      split_string(opts.ref.c_str(), ',', reflist);
      triple<short> ref;
      if (reflist.size() == 1)
      {      
        ref.v1 = ref.v2 = ref.v3 = atoi(reflist[0].c_str());
      }
      else if (reflist.size() == 3)
      {
        ref.v1 = atoi(reflist[0].c_str());
        ref.v2 = atoi(reflist[1].c_str());
        ref.v3 = atoi(reflist[2].c_str());
      }
      else 
      {
        std::cerr << std::endl << "argument error: '" << ref_par.c_str() << "' expects one single int or a int triple!" << std::endl;
        print_itk_refine_help();
        return;
      }

      itk_image img;

      std::cout << "Reading image: " << opts.msh << std::endl;
      gettimeofday(&t1, NULL);
      img.read_file(opts.msh.c_str());
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Normalizing voxel spacing .. " << std::endl;
      gettimeofday(&t1, NULL);
      img.normalize_spacing();
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      
      if (((ref.v1 > 0) && (ref.v2 > 0) && (ref.v3 > 0)) && 
          ((ref.v1 > 1) || (ref.v2 > 1) || (ref.v3 > 1))) 
      {
        std::cout << "Refine " << ref.v1 << "-" << ref.v2 << "-" << ref.v3 << " x image data .." << std::endl;
        gettimeofday(&t1, NULL);
        img.refine(ref);
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      }

      std::cout << "Writing output .. " << std::endl;
      gettimeofday(&t1, NULL);
      img.write_file(opts.outmsh.c_str());
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      break;
    }

    case ITK_PADDING:
    {
      itk_image img;

      mt_vector<std::string> list;
      split_string(opts.size, ',', list);

      if (list.size() == 3)
      {
        char * endptr;
        triple<unsigned int> size;
        size.v1 = strtoul(list[0].c_str(), &endptr, 10);
        size.v2 = strtoul(list[1].c_str(), &endptr, 10);
        size.v3 = strtoul(list[2].c_str(), &endptr, 10);

        std::cout << "Reading image: " << opts.msh << std::endl;
        gettimeofday(&t1, NULL);
        img.read_file(opts.msh.c_str());
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

        std::cout << "Add padding to voxel data .. " << std::endl;
        gettimeofday(&t1, NULL);
        img.padding(size);
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

        std::cout << "Writing output .. " << std::endl;
        gettimeofday(&t1, NULL);
        img.write_file(opts.outmsh.c_str());
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      }
      else
        std::cerr << "Error, check padding size format, ( 'i,j,k' expected )" << std::endl;

      break;
    }

    case ITK_CROP:
    {
      triple<short> ref = {1, 1, 1};
      if (opts.ref.size() > 0)
      {
        mt_vector<std::string> reflist;
        split_string(opts.ref.c_str(), ',', reflist);
        if (reflist.size() == 1)
        {      
          ref.v1 = ref.v2 = ref.v3 = atoi(reflist[0].c_str());
        }
        else if (reflist.size() == 3)
        {
          ref.v1 = atoi(reflist[0].c_str());
          ref.v2 = atoi(reflist[1].c_str());
          ref.v3 = atoi(reflist[2].c_str());
        }
        else 
        {
          std::cerr << std::endl << "argument error: '" << ref_par.c_str() << "' expects one single int or a int triple!" << std::endl;
          print_itk_crop_help();
          return;
        }
      }

      itk_image img;

      std::cout << "Reading image: " << opts.msh << std::endl;
      gettimeofday(&t1, NULL);
      img.read_file(opts.msh.c_str());
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Crop voxel data .. " << std::endl;
      gettimeofday(&t1, NULL);
      img.crop();
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      if (((ref.v1 > 0) && (ref.v2 > 0) && (ref.v3 > 0)) && 
          ((ref.v1 > 1) || (ref.v2 > 1) || (ref.v3 > 1))) 
      {
        std::cout << "Refine " << ref.v1 << "-" << ref.v2 << "-" << ref.v3 << " x image data .." << std::endl;
        gettimeofday(&t1, NULL);
        img.refine(ref);
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      }

      std::cout << "Writing output .. " << std::endl;
      gettimeofday(&t1, NULL);
      img.write_file(opts.outmsh.c_str());
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      break;
    }

    case ITK_FLIP:
    {
      itk_image img;

      unsigned int axes = ITK_AXIS_0;
      for (unsigned int i = 0; i < opts.axes.length(); i++)
      {
        switch (opts.axes[i])
        {
          case 'x': axes |= ITK_AXIS_X; break;
          case 'y': axes |= ITK_AXIS_Y; break;
          case 'z': axes |= ITK_AXIS_Z; break;
          default : break;
        }
      }

      std::cout << "Reading image: " << opts.msh << std::endl;
      gettimeofday(&t1, NULL);
      img.read_file(opts.msh.c_str());
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Flip voxel data .. " << std::endl;
      gettimeofday(&t1, NULL);
      img.flip(axes);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Writing output .. " << std::endl;
      gettimeofday(&t1, NULL);
      img.write_file(opts.outmsh.c_str());
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      break;
    }

    case ITK_DTYPE:
    {
      itk_image inp_img, out_img;
      int dtype = atoi(opts.dtype.c_str());

      if ((dtype > 0) && (dtype < 11))
      {
        std::cout << "Reading image: " << opts.msh << std::endl;
        gettimeofday(&t1, NULL);
        inp_img.read_file(opts.msh.c_str());
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

        std::cout << "Converting image data .. " << std::endl;
        gettimeofday(&t1, NULL);
        out_img.assign(inp_img, itk_get_datatype_str(dtype));
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

        std::cout << "Writing output .. " << std::endl;
        gettimeofday(&t1, NULL);
        out_img.write_file(opts.outmsh.c_str());
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      }
      else if (dtype == 11)
      {
        std::cout << "Reading image: " << opts.msh << std::endl;
        gettimeofday(&t1, NULL);
        inp_img.read_file(opts.msh.c_str());
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

        std::cout << "Converting image data .. " << std::endl;
        gettimeofday(&t1, NULL);
        out_img.assign(inp_img, itk_get_datatype_str(1));
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

        std::cout << "Writing output .. " << std::endl;
        gettimeofday(&t1, NULL);
        out_img.write_file_color_scalars(opts.outmsh.c_str());
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      }
      else
        std::cerr << "Error, invalid dtype " << dtype << std::endl;

      break;
    }

    case ITK_REFINE:
    {
      triple<short> ref = {1, 1, 1};
      if (opts.ref.size() > 0)
      {
        mt_vector<std::string> reflist;
        split_string(opts.ref.c_str(), ',', reflist);
        if (reflist.size() == 1)
        {      
          ref.v1 = ref.v2 = ref.v3 = atoi(reflist[0].c_str());
        }
        else if (reflist.size() == 3)
        {
          ref.v1 = atoi(reflist[0].c_str());
          ref.v2 = atoi(reflist[1].c_str());
          ref.v3 = atoi(reflist[2].c_str());
        }
        else 
        {
          std::cerr << std::endl << "argument error: '" << ref_par.c_str() << "' expects one single int or a int triple!" << std::endl;
          print_itk_refine_help();
          return;
        }
      }

      itk_image img;

      std::cout << "Reading image: " << opts.msh << std::endl;
      gettimeofday(&t1, NULL);
      img.read_file(opts.msh.c_str());
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      if (((ref.v1 > 0) && (ref.v2 > 0) && (ref.v3 > 0)) &&
          ((ref.v1 > 1) || (ref.v2 > 1) || (ref.v3 > 1)))      
      {
        std::cout << "Refine " << ref.v1 << "-" << ref.v2 << "-" << ref.v3 << " x image data .." << std::endl;
        gettimeofday(&t1, NULL);
        img.refine(ref);
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      }

      std::cout << "Writing output .. " << std::endl;
      gettimeofday(&t1, NULL);
      img.write_file(opts.outmsh.c_str());
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      break;
    }

    case ITK_SAMPLE:
    {
      // set up itk image space from params
      itk_image img;

      if(opts.msh.size()) {
        std::cout << "Reading image stack " << opts.msh << " .." << std::endl;
        gettimeofday(&t1, NULL);
        itk_image rimg;
        rimg.read_file(opts.msh.c_str());
        img.assign(rimg, "short");
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      }
      else {
        triple<double> orig, max, spacing;
        std::cout << "Setting up empty image stack .." << std::endl;
        gettimeofday(&t1, NULL);
        mt_vector<std::string> minstr, maxstr;
        split_string(opts.min, ',', minstr);
        split_string(opts.max, ',', maxstr);

        assert(minstr.size() == 3);
        assert(maxstr.size() == 3);

        orig.v1 = atof(minstr[0].c_str());
        orig.v2 = atof(minstr[1].c_str());
        orig.v3 = atof(minstr[2].c_str());

        max.v1 = atof(maxstr[0].c_str());
        max.v2 = atof(maxstr[1].c_str());
        max.v3 = atof(maxstr[2].c_str());

        float size = atof(opts.size.c_str());
        triple<unsigned int> dim;
        dim.v1 = fabs((max.v1 - orig.v1) / size) + 0.5;
        dim.v2 = fabs((max.v2 - orig.v2) / size) + 0.5;
        dim.v3 = fabs((max.v3 - orig.v3) / size) + 0.5;

        spacing.v1 = fabs((max.v1 - orig.v1) / dim.v1);
        spacing.v2 = fabs((max.v2 - orig.v2) / dim.v2);
        spacing.v3 = fabs((max.v3 - orig.v3) / dim.v3);

        img.assign(dim, orig, spacing, 1, "short");
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      }

      // read surface into mesh and set up kdtree
      std::cout << "Reading surfaces " << opts.surf << " and setting up kdtrees .." << std::endl;
      gettimeofday(&t1, NULL);

      mt_vector<std::string> surf_files;
      split_string(opts.surf, ',' , surf_files);

      int num_surfs = surf_files.size();
      mt_vector<mt_meshdata> meshes(num_surfs);
      mt_vector<kdtree>      trees(num_surfs);
      const int min_tris_per_leaf = 10;

      for(int s=0; s<num_surfs; s++)
      {
        mt_filename surffile(surf_files[s], "");
        read_mesh_selected(meshes[s], surffile.format, surffile.base, CRP_READ_PTS | CRP_READ_ELEM);
        trees[s].items_per_leaf = min_tris_per_leaf;
        trees[s].build_tree(meshes[s]);
      }
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Sampling image stack .." << std::endl;
      gettimeofday(&t1, NULL);
      for(int s=0; s<num_surfs; s++)
      {
        itk_sample(img, meshes[s], trees[s]);
      }
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Writing image stack " << opts.outmsh << " .." << std::endl;
      img.write_file(opts.outmsh.c_str());

      break;
    }

    case ITK_EXTRACT:
    {
      itk_image img;

      unsigned int plane = ITK_AXIS_0;
      if (opts.axes == "xy") plane = ITK_AXIS_X | ITK_AXIS_Y;
      else if (opts.axes == "xz") plane = ITK_AXIS_X | ITK_AXIS_Z;
      else if (opts.axes == "yz") plane = ITK_AXIS_Y | ITK_AXIS_Z;
      
      if (plane == ITK_AXIS_0)
      {
        std::cerr << "Error, invalid plane " << opts.axes << " !" << std::endl;
        return;
      }

      mt_vector<std::string> str_idx;
      split_string(opts.idx, ',' , str_idx); 
      std::set<unsigned int> idx;
      for (size_t i = 0; i < str_idx.size(); i++)
        idx.insert(std::stoul(str_idx[i]));

      std::cout << "Reading image: " << opts.msh << std::endl;
      gettimeofday(&t1, NULL);
      img.read_file(opts.msh.c_str());
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Extract voxel data .. " << std::endl;
      gettimeofday(&t1, NULL);
      img.extract_slices(plane, idx);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Writing output .. " << std::endl;
      gettimeofday(&t1, NULL);
      img.write_file(opts.outmsh.c_str());
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      break;
    }
  }
}
