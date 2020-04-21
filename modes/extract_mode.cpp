/**
* @file extract_mode.h
* @brief Meshtool extract mode.
*
* This mode allows to extract many entities (submeshes, surfaces, etc.) from a mesh.
*
* @author Aurel Neic
* @version
* @date 2016-12-13
*/

#include "mt_modes_base.h"

#ifdef MT_ADDONS
#include "addons_utils.h"
#endif

static const std::string hyp_par = "-hybrid=";
static const std::string ang_thr_par = "-ang_thr=";
static const std::string tagfile_par = "-tag_file=";

/**
* @brief Extract mode options.
*/
struct extract_options {
  enum op_type type;
  std::string msh_base;
  std::string msh1_base;
  std::string msh2_base;
  std::string submsh_base;
  std::string msh_dat_file;
  std::string submsh_dat_file;
  std::string oper;
  std::string surf;
  std::string coord;
  std::string ifmt;
  std::string ofmt;
  std::string edge_thr;
  std::string rad;
  std::string lower_rad;
  std::string idat;
  std::string odat;
  std::string mode;
  std::string hybrid;
  std::string thr;
  std::string ang_thr;
  std::string tag_file;
  std::set<int> tags;
};

/**
* @brief Display extract mesh help message.
*/
void print_extract_mesh_help()
{
  fprintf(stderr, "extract mesh: a submesh is extracted from a given mesh based on given element tags\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t (input) path to basename of the mesh to extract from\n", mesh_par.c_str());
  fprintf(stderr, "%stag1%ctag2\t\t (input) \"%c\"-seperated list of tags\n", tags_par.c_str(), tag_separator, tag_separator);
  fprintf(stderr, "%s<path>\t (optional) path to an alternative tag file {*.tags, *.btags}.\n", tagfile_par.c_str());
  fprintf(stderr, "%s<path>\t\t (output) path to basename of the submesh to extract to\n", submesh_par.c_str());
  fprintf(stderr, "%s<format>\t\t (optional) mesh input format.\n", inp_format_par.c_str());
  fprintf(stderr, "%s<format>\t\t (optional) mesh output format.\n\n", out_format_par.c_str());
  fprintf(stderr, "The supported input formats are:\n%s\n", input_formats.c_str());
  fprintf(stderr, "The supported output formats are:\n%s\n", output_formats.c_str());
  fprintf(stderr, "\n");
}

/**
* @brief Display extract data help message.
*/
void print_extract_data_help()
{
  fprintf(stderr, "extract data: data defined on a mesh is extracted for a given submesh\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t (input) path to basename of the submesh\n", submesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) file the data is extracted from.\n", mesh_data_par.c_str());
  fprintf(stderr, "%s<path>\t (output) file the data is extracted into\n", submesh_data_par.c_str());
  fprintf(stderr, "%s<int>\t\t (optional) Data mode. 0 = nodal, 1 = element. Default is 0.\n", mode_par.c_str());
  fprintf(stderr, "\n");
  fprintf(stderr, "Note that the files defining the submesh must include a *.eidx and a *.nod file.\n"
                  "This files define how elements and nodes of the submesh map back into the original mesh.\n"
                  "The *.eidx and *.nod files are generated when using the \"extract mesh\" mode.\n");
  fprintf(stderr, "\n");
}

/**
* @brief Display extract surface help message.
*/
void print_extract_surf_help()
{
  fprintf(stderr, "extract surface: extract a sequence of surfaces defined by set operations on element tags\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t (input) path to basename of the mesh\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t\t (output) list of names associated to the given operations\n", surf_par.c_str());
  fprintf(stderr, "%soperations\t\t (optional) list of operations to perform. By default, the\n"
                  "\t\t\t surface of the full mesh is computed.\n", oper_par.c_str());
  fprintf(stderr, "%s<path>\t (optional) path to an alternative tag file {*.tags, *.btags}.\n", tagfile_par.c_str());
  fprintf(stderr, "%s<deg. angle>\t (optional) surface elements connected to sharp edges will\n"
                  "\t\t\t be removed. A sharp edge is defined by the nodes which connect elements\n"
                  "\t\t\t with normal vectors at angles above the given threshold.\n", edge_par.c_str());
  fprintf(stderr, "%s<deg. angle>\t (optional) if set, surface traversal stops when angle between\n"
                  "\t\t\t current and starting normal vectors exceeds threshold.\n", ang_thr_par.c_str());
  fprintf(stderr, "%s<xyz>:<xyz>:..\t (optional) restrict surfaces to those elements reachable by\n"
                  "\t\t\t surface edge-traversal from the surface vertices closest to the given\n"
                  "\t\t\t coordinates. If %s is also provided, sharp edges will block\n"
                  "\t\t\t traversal, thus limit what is reachable.\n", coord_par.c_str(), edge_par.c_str());
  fprintf(stderr, "%s<float>\t\t (optional) surface edge-traversal is limited to the given\n"
                  "\t\t\t radius from the initial index.\n", size_par.c_str());
  fprintf(stderr, "%s<float>\t (optional) surface edge-traversal limitation lower size (for extracting bands).\n", lower_size_par.c_str());
  fprintf(stderr, "%s<int>\t\t (optional) Write hybrid quad + tri surfaces. 1 == on, 0 == off. 0 is default.\n", hyp_par.c_str());
  fprintf(stderr, "%s<format>\t\t (optional) mesh input format.\n", inp_format_par.c_str());
  fprintf(stderr, "%s<format>\t\t (optional) mesh output format. If set, the surfaces will also\n"
                  "\t\t\t be written as surface meshes.\n\n", out_format_par.c_str());
  fprintf(stderr, "The supported input formats are:\n%s\n", input_formats.c_str());
  fprintf(stderr, "The supported output formats are:\n%s\n", output_formats.c_str());
  fprintf(stderr, "\n");
  fprintf(stderr, "The format of the operations is:\n");
  fprintf(stderr, "tagA1,tagA2,[surfA1,surfA2..]..[+-:]tagB1,tagB2,[surfB1..]..;tagA1,..[+-:]tagB1..;..\n");
  fprintf(stderr, "Tag regions separated by \",\" will be unified into submeshes and their surface computed.\n"
                  "Alternatively, surfaces can be provided directly by .surf surface files (only basename, no extension).\n"
                  "If two surfaces are separated by \"-\", the rhs surface will be removed from the\n"
                  "lhs surface (set difference). Similarly, using \"+\" will compute the surface union.\n"
                  "If the submeshes are separated by \":\", the set intersection of the two submesh surfaces will be computed.\n"
                  "Individual operations are separated by \";\".\n\n");
  fprintf(stderr, "The number of names provided with \"%s\" must match the number of operations. If no operations are provided,\n"
                  "the surface of the whole geometry will be extracted. Individual names are separated by \",\".\n\n", surf_par.c_str());
  fprintf(stderr, "Further restrictions can be added to the surface extraction with the %s , %s , %s options.\n", edge_par.c_str(), coord_par.c_str(), size_par.c_str());
  fprintf(stderr, "\n");
}

/**
* @brief Display extract myocard help message.
*/
void print_extract_myo_help()
{
  fprintf(stderr, "extract myocard: the myocardium is extracted from a given mesh. The myocard is identified based on non-zero fibers.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh to extract from\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to basename of the submesh to extract to\n", submesh_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh input format.\n", inp_format_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh output format.\n\n", out_format_par.c_str());
  fprintf(stderr, "The supported input formats are:\n%s\n", input_formats.c_str());
  fprintf(stderr, "The supported output formats are:\n%s\n", output_formats.c_str());
  fprintf(stderr, "\n");
}

/**
* @brief Display extract volume help message.
*/
void print_extract_volume_help()
{
  fprintf(stderr, "extract volume: extract elements inside a given box volume.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t\t (input) path to basename of the mesh to extract from\n", mesh_par.c_str());
  fprintf(stderr, "%sx,xd,y,yd,z,zd\t (input) volume definition\n", coord_par.c_str());
  fprintf(stderr, "%s<path>\t\t (output) path to basename of the submesh to extract to\n", submesh_par.c_str());
  fprintf(stderr, "%s<int>\t\t (optional) mesh output mode. 0 = submesh, 1 = vtx file. Default is 0.\n", mode_par.c_str());
  fprintf(stderr, "%s<format>\t\t (optional) mesh input format.\n", inp_format_par.c_str());
  fprintf(stderr, "%s<format>\t\t (optional) mesh output format.\n\n", out_format_par.c_str());
  fprintf(stderr, "The supported input formats are:\n%s\n", input_formats.c_str());
  fprintf(stderr, "The supported output formats are:\n%s\n", output_formats.c_str());
  fprintf(stderr, "\n");
}

/**
* @brief Display extract unreachable help message.
*/
void print_extract_unreachable_help()
{
  fprintf(stderr, "extract unreachable: extract elements that are unreachable when traversing the mesh connectivity from a given start vertex.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh to extract from\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to basename of the submesh to extract to\n", submesh_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh input format.\n", inp_format_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh output format.\n\n", out_format_par.c_str());
  fprintf(stderr, "The supported input formats are:\n%s\n", input_formats.c_str());
  fprintf(stderr, "The supported output formats are:\n%s\n", output_formats.c_str());
  fprintf(stderr, "\n");
}

void print_extract_gradient_help()
{
  fprintf(stderr, "extract gradient: compute gradient and gradient magnitude of a scalar function on a mesh.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh to extract from\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to input data. Also specify file extension.\n", idat_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to output data.\n", odat_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh input format.\n", inp_format_par.c_str());
  fprintf(stderr, "%s<int>\t (optional) output mode. 0 == nodal putput, 1 == element output. 0 is default.\n", mode_par.c_str());
  fprintf(stderr, "The supported input formats are:\n%s\n", input_formats.c_str());
  fprintf(stderr, "\n");
}

void print_extract_overlap_help()
{
  fprintf(stderr, "extract overlap: extract the elements of mesh1 that overlap with mesh2.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of mesh1.\n", mesh1_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to basename of mesh2.\n", mesh2_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to submesh (the overlaping elements).\n", submesh_par.c_str());
  fprintf(stderr, "%s<int>\t (optional) output mode. 0 == only overlap, 1 == overlap and complement.\n"
                  "\t\t 0 is default.\n", mode_par.c_str());
  fprintf(stderr, "%s<float>\t (optional) overlap region radius, in multiples of the average mesh1 edge length.\n"
                  "\t\t default is 2.0.\n", size_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh input format.\n", inp_format_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh output format.\n\n", out_format_par.c_str());
  fprintf(stderr, "The supported input formats are:\n%s\n", input_formats.c_str());
  fprintf(stderr, "The supported output formats are:\n%s\n", output_formats.c_str());
  fprintf(stderr, "\n");
}

void print_extract_tags_help()
{
  fprintf(stderr, "extract tags: extract the elements tags into an element data vector.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the mesh.\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) .dat file the tags are extracted to.\n", odat_par.c_str());
  fprintf(stderr, "%s<format>\t (optional) mesh input format.\n", inp_format_par.c_str());
  fprintf(stderr, "The supported input formats are:\n%s\n", input_formats.c_str());
  fprintf(stderr, "\n");
}

#ifdef MT_ADDONS
void print_extract_isosurf_help()
{
  fprintf(stderr, "extract isosurf: extract an isosurface.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of mesh.\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (input) path to real data.\n", idat_par.c_str());
  fprintf(stderr, "%s<float>\t (input) value to compute isosurf for.\n", thr_par.c_str());
  fprintf(stderr, "%s<name>\t (output) path to basename of output surface.\n", surf_par.c_str());
  fprintf(stderr, "\n");
}
#endif

/**
* @brief Parse extract mode options.
*
* @param [in]  argc Arguments count.
* @param [in]  argv Arguments string-array.
* @param [out] opts Options structure.
*
* @return 0 for success, >0 otherwise.
*/
int extract_parse_options(int argc, char** argv, struct extract_options & opts)
{
  if(argc < 3) {
    std::cerr << "Please choose one of the following extract modes: " << std::endl << std::endl;
    print_extract_mesh_help();
    print_extract_data_help();
    print_extract_surf_help();
    print_extract_myo_help();
    print_extract_volume_help();
    print_extract_unreachable_help();
    print_extract_gradient_help();
    print_extract_overlap_help();
    print_extract_tags_help();
#ifdef MT_ADDONS
    print_extract_isosurf_help();
#endif
    return 1;
  }

  std::string extract_type = argv[2];

  // parse all extract parameters -----------------------------------------------------------------
  for(int i=3; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh_par, opts.msh_base);
    if(!match) match = parse_param(param, mesh1_par, opts.msh1_base);
    if(!match) match = parse_param(param, mesh2_par, opts.msh2_base);
    if(!match) match = parse_param(param, submesh_par, opts.submsh_base);
    if(!match) match = parse_param(param, mesh_data_par, opts.msh_dat_file, chk_fexists);
    if(!match) match = parse_param(param, submesh_data_par, opts.submsh_dat_file);
    if(!match) match = parse_param(param, oper_par, opts.oper);
    if(!match) match = parse_param(param, surf_par, opts.surf);
    if(!match) match = parse_param(param, inp_format_par, opts.ifmt);
    if(!match) match = parse_param(param, out_format_par, opts.ofmt);
    if(!match) match = parse_param(param, coord_par, opts.coord);
    if(!match) match = parse_param(param, edge_par, opts.edge_thr);
    if(!match) match = parse_param(param, ang_thr_par, opts.ang_thr);
    if(!match) match = parse_param(param, size_par, opts.rad);
    if(!match) match = parse_param(param, lower_size_par, opts.lower_rad);
    if(!match) match = parse_param(param, idat_par, opts.idat, chk_fexists);
    if(!match) match = parse_param(param, odat_par, opts.odat);
    if(!match) match = parse_param(param, mode_par, opts.mode);
    if(!match) match = parse_param(param, thr_par, opts.thr);
    if(!match) match = parse_param(param, hyp_par, opts.hybrid);
    if(!match) match = parse_param(param, tagfile_par, opts.tag_file, chk_fexists);


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
      return 2;
    }
  }
  fixBasename(opts.surf);

  // check if all relevant parameters have been set ---------------------------------------------------
  if(extract_type.compare("mesh") == 0)
  {
    opts.type = OP_MESH;

    if( !(opts.msh_base.size() > 0 && opts.submsh_base.size() > 0 && opts.tags.size() > 0) )
    {
      std::cerr << "Mesh extract error: Insufficient parameters provided." << std::endl;
      print_extract_mesh_help();
      return 4;
    }
  }
  else if(extract_type.compare("data") == 0)
  {
    bool submsh_set = opts.submsh_base.size() > 0, submsh_data_set = opts.submsh_dat_file.size() > 0,
         msh_data_set = opts.submsh_dat_file.size() > 0;

    opts.type = OP_DATA;

    if( !(submsh_set && msh_data_set && submsh_data_set) )
    {
      std::cerr << "Data extract error: Insufficient parameters provided." << std::endl;
      print_extract_data_help();
      return 4;
    }
  }
  else if(extract_type.compare("surface") == 0)
  {
    opts.type = OP_SURF;
    if(! (opts.msh_base.size() > 0 && opts.surf.size() > 0))
    {
      std::cerr << "Surface extract error: Insufficient parameters provided." << std::endl;
      print_extract_surf_help();
      return 4;
    }
  }
  else if(extract_type.compare("myocard") == 0)
  {
    opts.type = OP_MYO;
    if(! (opts.msh_base.size() > 0 && opts.submsh_base.size() > 0))
    {
      std::cerr << "Myocard extract error: Insufficient parameters provided." << std::endl;
      print_extract_myo_help();
      return 4;
    }
  }
  else if(extract_type.compare("volume") == 0)
  {
    opts.type = OP_VOL;
    if(! (opts.msh_base.size() > 0 && opts.submsh_base.size() > 0 && opts.coord.size() > 0) )
    {
      std::cerr << "Volume extract error: Insufficient parameters provided." << std::endl;
      print_extract_volume_help();
      return 4;
    }
  }
  else if(extract_type.compare("unreachable") == 0)
  {
    opts.type = OP_UR;
    if(! (opts.msh_base.size() > 0 && opts.submsh_base.size() > 0 ) )
    {
      std::cerr << "Unreachables extract error: Insufficient parameters provided." << std::endl;
      print_extract_unreachable_help();
      return 4;
    }
  }
  else if(extract_type.compare("gradient") == 0)
  {
    opts.type = OP_GRAD;
    if(opts.msh_base.size() == 0 || opts.idat.size() == 0 || opts.odat.size() == 0 )
    {
      std::cerr << "Gradient extract error: Insufficient parameters provided." << std::endl;
      print_extract_gradient_help();
      return 4;
    }
  }
  else if(extract_type.compare("overlap") == 0)
  {
    opts.type = OP_OVERLP;
    if(opts.msh1_base.size() == 0 || opts.msh2_base.size() == 0 || opts.submsh_base.size() == 0 )
    {
      std::cerr << "Overlap extract error: Insufficient parameters provided." << std::endl;
      print_extract_overlap_help();
      return 4;
    }
  }
  else if(extract_type.compare("tags") == 0)
  {
    opts.type = OP_TAGS;
    if(opts.msh_base.size() == 0 || opts.odat.size() == 0)
    {
      std::cerr << "Tags extract error: Insufficient parameters provided." << std::endl;
      print_extract_tags_help();
      return 4;
    }
  }
#ifdef MT_ADDONS
  else if(extract_type.compare("isosurf") == 0)
  {
    opts.type = OP_ISOSURF;
    if(opts.msh_base.size() == 0 || opts.idat.size() == 0 || opts.surf.size() == 0 ||
       opts.thr.size() == 0)
    {
      std::cerr << "Isosurface extract error: Insufficient parameters provided." << std::endl;
      print_extract_isosurf_help();
      return 4;
    }
  }
#endif
  else {
    print_usage(argv[0]);
    return 2;
  }

  return 0;
}

void tokenize_operstr(std::string & operstr, const char token_indicator, mt_vector<std::string> & file_names)
{
  int stringlen = operstr.length();

  file_names.resize(0);
  int file_names_idx = 0, pos = 0, start_pos;

  while (pos < stringlen) {

    if (operstr[pos] == '"') {
      // if we find a '"' char, the user is passing a protected filename. thus we search
      // the closing '"' and store the string between into file_names.
      start_pos = pos++;
      while ((pos < stringlen) && (operstr[pos] != '"'))
        pos++;

      // we found start and end pos characters, we do the tokenization
      if ((operstr[pos] == '"') && (pos-start_pos > 1)) {
        // extract filename and add to vector
        file_names.push_back(operstr.substr(start_pos+1, pos-start_pos-1));

        // create string token and replace in `operstr`
        std::string token = token_indicator + std::to_string(file_names_idx++);
        operstr.replace(start_pos, pos-start_pos+1, token);

        // update string length
        stringlen = operstr.size();

        // correct position
        pos = start_pos + token.length(); 
      }
      // we found start and end pos characters, but the string is empty
      else if ((operstr[pos] == '"') && (pos-start_pos == 1)) {
        fprintf(stderr, "%s error: Cannot tokenize string: %s \nEmpty string! Aborting!\n",
                __func__, operstr.c_str());
        exit(EXIT_FAILURE);
      }
      else {
        fprintf(stderr, "%s error: Cannot tokenize string: %s \nNon-matching \" characters! Aborting!\n",
                __func__, operstr.c_str());
        exit(EXIT_FAILURE);
      }

    }
    pos++;
  }
}

void detokenize_surfaces(mt_vector<mt_vector<std::string>> & surf,
                         const char token_indicator,
                         const mt_vector<std::string> & file_names)
{

  for (size_t i=0; i<surf.size(); i++)
    for (size_t j=0; j<surf[i].size(); j++)
      if ((surf[i][j].empty() == false) && (surf[i][j][0] == token_indicator)) {
        const int idx = std::stoi(surf[i][j]);
        surf[i][j] = file_names[idx];
      }
}


/**
* @brief Extract the tag-sets for the different surface extraction operations.
*
* @param [in]  operstr  The string defining the operations.
* @param [out] operflg  The encoded operation types.
* @param [out] tagsA    The lhs tags associated with each operation.
* @param [out] tagsB    The rhs tags associated with each operation.
*/
void extract_surf_fill_taglists(const std::string & operstr,
                                std::string & operflg,
                                mt_vector<mt_vector<std::string> > & surfA,
                                mt_vector<mt_vector<std::string> > & surfB)
{
   mt_vector<std::string> operlist;

   split_string(operstr, ';', operlist);
   int numoper = operlist.size();
   surfA.resize(numoper), surfB.resize(numoper);
   operflg.resize(numoper);

   for(int i=0; i<numoper; i++) {
     mt_vector<std::string> tlist;
     split_string(operlist[i], ':', tlist);
     if(tlist.size() == 2) {
       operflg[i] = ':';
       split_string(tlist[0], ',', surfA[i]);
       split_string(tlist[1], ',', surfB[i]);
     }
     else {
       split_string(operlist[i], '-', tlist);
       if(tlist.size() == 2) {
         operflg[i] = '-';
         split_string(tlist[0], ',', surfA[i]);
         split_string(tlist[1], ',', surfB[i]);
       }
       else {
         split_string(operlist[i], '+', tlist);
         if(tlist.size() == 2) {
           operflg[i] = '+';
           split_string(tlist[0], ',', surfA[i]);
           split_string(tlist[1], ',', surfB[i]);
         }
         else {
           operflg[i] = '0';
           split_string(tlist[0], ',', surfA[i]);
         }
       }
     }
   }
}

void unified_surf_from_mixed_list(const mt_meshdata & mesh,
                                  const mt_vector<std::string> & list,
                                  MT_MAP<triple<mt_int>,    tri_sele>  & tri_surf,
                                  MT_MAP<quadruple<mt_int>, quad_sele> & quad_surf)
{
  MT_USET<mt_int> tags;

  for(size_t i=0; i<list.size(); i++) {
    if(is_integer(list[i])) {
      tags.insert(atoi(list[i].c_str()));
    }
    else {
      mt_meshdata curr_surf;
      mt_vector<mt_int> eidx;

      read_surf_info(mesh, curr_surf, eidx, list[i]);
      surfmesh_to_surfmap(curr_surf, eidx, tri_surf, quad_surf);
    }
  }

  if(tags.size())
    compute_surface(mesh, tags, tri_surf, quad_surf);
}


template<class T>
void extract_1dpn(const mt_vector<mt_int> & nod, mt_vector<T> & in, mt_vector<T> & out)
{
  #ifdef OPENMP
  #pragma omp parallel for schedule(dynamic, 256)
  #endif
  for(size_t i=0; i<nod.size(); i++)
    out[i] = in[nod[i]];
}
template<class T>
void extract_3dpn(const mt_vector<mt_int> & nod, mt_vector<T> & in, mt_vector<T> & out)
{
  #ifdef OPENMP
  #pragma omp parallel for schedule(dynamic, 128)
  #endif
  for(size_t i=0; i<nod.size(); i++) {
    out[i*3+0] = in[nod[i]*3+0];
    out[i*3+1] = in[nod[i]*3+1];
    out[i*3+2] = in[nod[i]*3+2];
  }
}


/**
* @brief Extract mode function.
*
* @param [in]  argc Arguments count.
* @param [in]  argv Arguments string-array.
*/
void extract_mode(int argc, char** argv)
{
  struct extract_options opts;
  int ret = extract_parse_options(argc, argv, opts);
  struct timeval t1, t2;

  if(ret != 0) return;

  // we "centralize" the mesh reading to reduce code length ------------------------------
  struct mt_meshdata mesh;
  mt_filename mshfile(opts.msh_base, opts.ifmt);

  if(opts.type != OP_DATA && opts.type != OP_OVERLP) {
    std::cout << "Reading mesh: " << mshfile.base << std::endl;
    gettimeofday(&t1, NULL);
    read_mesh_selected(mesh, mshfile.format, mshfile.base);
    gettimeofday(&t2, NULL);
    std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
  }

  switch(opts.type) {

    case OP_MESH: {
      // extract mesh  ------------------------------------------------------------
      if (opts.tag_file.size())
        readElementTags_general(mesh.etags, opts.tag_file);

      std::cout << "Extracting tags: ";
      for(std::set<int>::iterator it = opts.tags.begin(); it != opts.tags.end(); ++it)
        std::cout << *it << " ";
      std::cout << std::endl;
      gettimeofday(&t1, NULL);

      mt_vector<mt_int> nod, eidx;

      // generate vector of elements to keep
      mt_vector<bool> keep(mesh.etags.size(), false);
      for(size_t i=0; i<keep.size(); i++)
        if(opts.tags.count(mesh.etags[i])) keep[i] = true;

      restrict_meshdata(keep, mesh, nod, eidx);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      mt_filename submshfile(opts.submsh_base, opts.ofmt);
      std::cout << "Writing mesh: " << submshfile.base << std::endl;
      gettimeofday(&t1, NULL);
      write_mesh_selected(mesh, submshfile.format, submshfile.base);
      binary_write(nod.begin(), nod.end(), submshfile.base + NOD_EXT);
      binary_write(eidx.begin(), eidx.end(), submshfile.base + EIDX_EXT);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      break;
    }

    case OP_DATA: {
      // extract data  ------------------------------------------------------------
      bool nodal_data = true;
      if(opts.mode.size()) nodal_data = atoi(opts.mode.c_str()) != 1;

      mt_vector<mt_int> indices;
      std::string indices_file;
      if(nodal_data)
        indices_file = opts.submsh_base + NOD_EXT;
      else
        indices_file = opts.submsh_base + EIDX_EXT;

      std::cout << "Reading " << indices_file << std::endl;
      binary_read(indices, indices_file);

      // read and insert data
      short data_idx = -1;
      mt_vector<mt_real> idat, odat;
      igb_header igb_msh, igb_submsh;

      setup_data_format(opts.msh_dat_file, opts.submsh_dat_file, data_idx, igb_msh, igb_submsh);

      gettimeofday(&t1, NULL);
      std::cout << "Extracting data .. " << std::endl;
      switch(data_idx) {
        case 0:
        case 2:
          read_vector_ascii(idat, opts.msh_dat_file, true);

          if(data_idx == 0) {
            odat.resize(indices.size());
            extract_1dpn(indices, idat, odat);
            write_vector_ascii(odat, opts.submsh_dat_file, 1);
          }
          else {
            odat.resize(indices.size()*3);
            extract_3dpn(indices, idat, odat);
            write_vector_ascii(odat, opts.submsh_dat_file, 3);
          }
          break;

        case 1:
        case 3:
        {
          if(igb_msh.v_t != igb_submsh.v_t) {
            fprintf(stderr, "Number of time-slices in mesh and submesh do not match. Aborting!\n");
            exit(1);
          }

          // set new spatial size and write output igb header
          igb_submsh.v_x = data_idx == 1 ? indices.size() : indices.size()*3;
          write_igb_header(igb_submsh);
          odat.resize(igb_submsh.v_x);

          std::vector<mt_real> rbuff;
          printf("Extracting igb time-slices %s into %s. \n", opts.msh_dat_file.c_str(),
                 opts.submsh_dat_file.c_str());

          for(int t=0; t<igb_submsh.v_t; t++) {
            printf("\rcurrent time-slice %d / %d .. ", t+1, int(igb_submsh.v_t));
            fflush(stdout);

            read_igb_slice(rbuff, igb_msh);
            idat.assign(rbuff.begin(), rbuff.end());

            if(data_idx == 1) extract_1dpn(indices, idat, odat);
            else              extract_3dpn(indices, idat, odat);

            rbuff.assign(odat.begin(), odat.end());
            write_igb_slice(rbuff, igb_submsh);
          }
          printf("\n");

          fclose(igb_msh.fileptr);
          fclose(igb_submsh.fileptr);
          break;
        }
        default: break;
      }

      break;
    }

    case OP_SURF: {
      // extract surface --------------------------------------------------------------
      if (opts.tag_file.size())
        readElementTags_general(mesh.etags, opts.tag_file);

      // check for sharp angle detection, we do this first to not waste the users time
      // when the input is actually bad
      float edge_thr = opts.edge_thr.size() > 0 ? atof(opts.edge_thr.c_str()) : 0.0f;
      float ang_thr  = opts.ang_thr.size() > 0  ? atof(opts.ang_thr.c_str())  : 0.0f;

      mt_vector<std::string> coords;
      if(opts.coord.size())
        split_string(opts.coord, ':' , coords);

      // some checks regarding edge_thr and ang_thr
      if(edge_thr && ang_thr) {
        fprintf(stderr, "%s error: %s and %s cannot be set at the same time! Aborting!\n",
                __func__, edge_par.c_str(), ang_thr_par.c_str());
        exit(EXIT_FAILURE);
      }

      if(ang_thr && coords.size() != 1) {
        fprintf(stderr, "%s error: %s requires exactly one starting coordinate! Aborting!\n",
                __func__, ang_thr_par.c_str());
        exit(EXIT_FAILURE);
      }

      int npts = *(std::max_element(mesh.e2n_con.begin(), mesh.e2n_con.end())) + 1;

      mt_vector<mt_vector<std::string> > surfA, surfB;
      std::string operflg;

      // we need full mesh connectivity as we might need to recover surface info
      std::cout << "Setting up n2e / n2n graphs .. " << std::endl;
      compute_full_mesh_connectivity(mesh, mshfile.base);

      // define operations
      if(opts.oper.size() > 0) {
        const char token_indicator = ' ';
        mt_vector<std::string> file_names;
        tokenize_operstr(opts.oper, token_indicator, file_names);

        extract_surf_fill_taglists(opts.oper, operflg, surfA, surfB);

        detokenize_surfaces(surfA, token_indicator, file_names);
        detokenize_surfaces(surfB, token_indicator, file_names);
      }
      else {
        operflg = "0";
        std::string list;

        surfA.resize(1);
        mt_vector<mt_int> alltags = mesh.etags;
        binary_sort(alltags); unique_resize(alltags);
        surfA[0].resize(alltags.size());

        for(size_t i=0; i<alltags.size(); i++)
          surfA[0][i] = std::to_string(alltags[i]);
      }

      // extract surface names
      mt_vector<std::string> names;
      split_string(opts.surf, ',', names);

      if(names.size() != operflg.size()) {
          fprintf(stderr, "Error: Number of provided names does not match number of operations!\n");
          return;
      }

      bool hybrid_output = false;
      if(opts.hybrid.size())
        hybrid_output = atoi(opts.hybrid.c_str()) == 1;

      int numoper = operflg.size();
      // mt_vector< mt_vector<mt_int> > surf_con(numoper);
      mt_vector<mt_meshdata> surf_mesh(numoper);
      mt_vector<nbc_data> nbc(numoper);

      gettimeofday(&t1, NULL);
      for(int i=0; i<numoper; i++) {
        std::cout << "Computing surface " << i+1 << "/" << numoper << std::endl;

        MT_MAP<triple<mt_int>,    tri_sele>  tri_surf_a,  tri_surf_b;
        MT_MAP<quadruple<mt_int>, quad_sele> quad_surf_a, quad_surf_b;

        unified_surf_from_mixed_list(mesh, surfA[i], tri_surf_a, quad_surf_a);
        if(operflg[i] != '0')
          unified_surf_from_mixed_list(mesh, surfB[i], tri_surf_b, quad_surf_b);

        switch(operflg[i]) {
          default: break;

          case '-':
            surface_difference(tri_surf_a, quad_surf_a, tri_surf_b, quad_surf_b);
            break;
          case '+':
            surface_union(tri_surf_a, quad_surf_a, tri_surf_b, quad_surf_b);
            break;
          case ':':
            surface_intersection(tri_surf_a, quad_surf_a, tri_surf_b, quad_surf_b);
            break;
        }

        mt_vector<mt_int> elem_origin;

        // Later we want to do hybrid only, but for now we distinguish
        if(hybrid_output) {
          // we generate a surfmesh of quads and tris
          surfmap_to_surfmesh(tri_surf_a, quad_surf_a, surf_mesh[i], elem_origin);
        }
        else {
          // we convert all surface elements to triangles
          surfmap_to_vector(tri_surf_a, quad_surf_a, surf_mesh[i].e2n_con, elem_origin);
          surf_mesh[i].e2n_cnt.assign(surf_mesh[i].e2n_con.size() / 3, 3);
          surf_mesh[i].etype.assign(surf_mesh[i].e2n_con.size() / 3, Tri);
        }

        bucket_sort_offset(surf_mesh[i].e2n_cnt, surf_mesh[i].e2n_dsp);
        generate_nbc_data(mesh, surf_mesh[i], elem_origin, nbc[i]);
      }
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      int extr_error = 0;
      gettimeofday(&t1, NULL);
      for(int i=0; i<numoper; i++) {
        if(surf_mesh[i].e2n_con.size() == 0) {
          std::cout << "Warning: surface extraction " << i
                    << " yielded an empty surface! Skipping!" << std::endl;
          extr_error = 1;
          continue;
        }

        if(coords.size() > 0)
        {
          compute_full_mesh_connectivity(surf_mesh[i]);

          MT_USET<mt_int> surf_nod;
          surf_nod.insert(surf_mesh[i].e2n_con.begin(), surf_mesh[i].e2n_con.end());

          if(edge_thr || ang_thr ) {
            // we block connectiviy at sharp edges. therefore we have to traverse along
            // the surface elements and check for blocking edges (not blocking nodes)
            mt_mapping<mt_int>            ele2edge;
            MT_MAP<tuple<mt_int>, mt_int> edges;
            compute_edges(surf_mesh[i], ele2edge, edges);
            ele2edge.transpose(); ele2edge.setup_dsp();

            MT_USET<mt_int> sharp_edges;
            mt_vector<mt_real> elem_nrmls;

            if(edge_thr)
              identify_sharp_edges(surf_mesh[i], mesh.xyz, ele2edge, edges, edge_thr, sharp_edges);
            else if(ang_thr)
              compute_element_surface_normals(surf_mesh[i], mesh.xyz, elem_nrmls);

            mt_vector<bool> reached(surf_mesh[i].e2n_cnt.size(), false);
            vec3r ref_nrml;

            // we insert elements connected to input coords
            for(size_t cidx=0; cidx < coords.size(); cidx++)
            {
              mt_vector<std::string> b;
              split_string(coords[cidx], ',' , b);
              if(b.size() == 3) {
                vec3r p(atof(b[0].c_str()), atof(b[1].c_str()), atof(b[2].c_str()));

                mt_int inp_idx; mt_real dist;
                linear_search_vtx(mesh.xyz, surf_nod, p, inp_idx, dist);

                // we insert the first element connected to the identified vertex into reached
                mt_int eidx = surf_mesh[i].n2e_con[surf_mesh[i].n2e_dsp[inp_idx]];
                reached[eidx] = true;

                if(ang_thr) ref_nrml.get(elem_nrmls.data() + eidx*3);
              }
              else {
                fprintf(stderr, "Error: Failed to parse coord %s ! Aborting!\n", coords[cidx].c_str());
                exit(1);
              }
            }

            if(edge_thr)
              traverse_surfelem_connectivity(surf_mesh[i], ele2edge, sharp_edges, reached);
            else if(ang_thr)
              traverse_surfelem_connectivity(surf_mesh[i], ele2edge, elem_nrmls, ref_nrml,
                                             ang_thr, reached);

            remove_elems_from_surf(reached, surf_mesh[i].e2n_cnt, surf_mesh[i].e2n_con, nbc[i]);
          }
          else {
            // we do not have to care about sharp edges, therefore we can use the cheaper
            // traverse_nodal_connectivity functionality
            MT_USET<mt_int> reached_set, remove_set;
            size_t numnodes = surf_mesh[i].n2n_cnt.size();
            mt_vector<bool> reached      (numnodes, false);
            mt_vector<bool> lower_reached(numnodes, false);

            // we loop over coords and traverse
            for(size_t cidx=0; cidx<coords.size(); cidx++)
            {
              mt_vector<std::string> b;
              split_string(coords[cidx], ',' , b);
              if(b.size() == 3) {
                vec3r   p(atof(b[0].c_str()), atof(b[1].c_str()), atof(b[2].c_str()));

                mt_int inp_idx; mt_real dist;
                linear_search_vtx(mesh.xyz, surf_nod, p, inp_idx, dist);
                reached[inp_idx] = true, lower_reached[inp_idx] = true;

                if(opts.rad.size() > 0) {
                  // we traverse only a given distance
                  mt_real rad = atof(opts.rad.c_str());
                  vec3r   refpt(mesh.xyz.data() + inp_idx*3);
                  traverse_nodal_connectivity(surf_mesh[i].n2n_cnt, surf_mesh[i].n2n_dsp,
                                        surf_mesh[i].n2n_con, mesh.xyz, refpt, rad,
                                        reached);

                  if(opts.lower_rad.size() > 0) {
                    mt_real lower_rad = atof(opts.lower_rad.c_str());
                    check_condition(lower_rad < rad, "Lower radius smaller than upper radius",
                                    __func__);

                    traverse_nodal_connectivity(surf_mesh[i].n2n_cnt, surf_mesh[i].n2n_dsp,
                                          surf_mesh[i].n2n_con, mesh.xyz, refpt, lower_rad,
                                          lower_reached);

                    for(size_t r=0; r < lower_reached.size(); r++)
                      if(lower_reached[r])
                        reached[r] = false;
                  }
                }
                else {
                  // we traverse whole connectivity
                  traverse_nodal_connectivity(surf_mesh[i].n2n_cnt, surf_mesh[i].n2n_dsp,
                                        surf_mesh[i].n2n_con, reached);
                }

                // store reached nodes and reset reached boolean vec
                for(size_t r=0; r<reached.size(); r++)
                  if(reached[r]) reached_set.insert(r);

                reached.assign(surf_mesh[i].n2n_cnt.size(), false);
                lower_reached.assign(surf_mesh[i].n2n_cnt.size(), false);
              }
              else {
                fprintf(stderr, "Error: Failed to parse coord %s ! Aborting!\n", coords[cidx].c_str());
                exit(1);
              }
            }

            for(mt_int r=0; r<mt_int(numnodes); r++)
              if(reached_set.count(r) == 0) remove_set.insert(r);

            // remove nodes
            remove_nodes_from_surf(remove_set, surf_mesh[i].e2n_cnt, surf_mesh[i].e2n_con, nbc[i]);
          }
        }

        // currently neubc files are only set up for triangle surf elems, thus we dont
        // write them for hybrid surfs
        bool write_nbc = !hybrid_output;
        write_surf_info(surf_mesh[i], write_nbc ? &nbc[i] : NULL, mesh.e2n_cnt.size(), npts, names[i]);

        if(opts.ofmt.size() > 0)
        {
          mt_vector<mt_int> vtx;
          vtx.assign(surf_mesh[i].e2n_con.begin(), surf_mesh[i].e2n_con.end());
          binary_sort(vtx); unique_resize(vtx);

          // nbc_data_into_surface changes the surface numbering! we must write all indexing
          // relative to the original mesh to disk before this point!!
          std::string surfmesh_file = names[i] + ".surfmesh";
          nbc_data_into_surface(mesh, nbc[i], surf_mesh[i]);

          std::cout << "Writing surface mesh " << surfmesh_file << std::endl;
          write_mesh_selected(surf_mesh[i], opts.ofmt, surfmesh_file);

          surfmesh_file += NOD_EXT;
          binary_write(vtx.begin(), vtx.end(), surfmesh_file);
        }
      }
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      if(extr_error) exit(extr_error);
      break;
    }

#ifdef MT_ADDONS
    case OP_ISOSURF:
    {
      compute_full_mesh_connectivity(mesh, mshfile.base);

      // read data
      mt_vector<mt_real> data;
      read_vector_ascii(data, opts.idat);

      // read threshold value
      mt_real thr = atof(opts.thr.c_str());

      mt_meshdata surf;
      compute_crinkled_isosurf(mesh, data, thr, surf);

      std::string surfname = opts.surf + SURF_EXT;

      std::cout << "Writing surface " << surfname << std::endl;
      write_surf(surf.e2n_cnt, surf.e2n_con, surfname);

      surf.xyz = mesh.xyz;
      reindex_nodes(surf);

      surfname = opts.surf + ".surfmesh";
      std::cout << "Writing mesh " << surfname << std::endl;
      write_mesh_selected(surf, "vtk_bin", surfname);
      break;
    }
#endif

    case OP_MYO:
    {
      // extract myocard --------------------------------------------------------------
      std::cout << "Extracting myocardium .." << std::endl;
      gettimeofday(&t1, NULL);
      mt_vector<mt_int> nod, eidx;
      bool twoFib = (mesh.lon.size() == mesh.e2n_cnt.size()*6);

      // generate vector of elements to keep
      mt_vector<bool> keep(mesh.etags.size(), false);
      if(twoFib) {
        #ifdef OPENMP
        #pragma omp parallel for schedule(guided, 100)
        #endif
        for(size_t i=0; i<keep.size(); i++) {
          float l1 = mesh.lon[i*6+0];
          float l2 = mesh.lon[i*6+1];
          float l3 = mesh.lon[i*6+2];
          if( !(l1 == 0.0f && l2 == 0.0f && l3 == 0.0f) ) keep[i] = true;
        }
      }
      else {
        #ifdef OPENMP
        #pragma omp parallel for schedule(guided, 100)
        #endif
        for(size_t i=0; i<keep.size(); i++) {
          float l1 = mesh.lon[i*3+0];
          float l2 = mesh.lon[i*3+1];
          float l3 = mesh.lon[i*3+2];
          if( !(l1 == 0.0f && l2 == 0.0f && l3 == 0.0f) ) keep[i] = true;
        }
      }
      restrict_meshdata(keep, mesh, nod, eidx);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      mt_filename submshfile(opts.submsh_base, opts.ofmt);
      std::cout << "Writing mesh: " << submshfile.base << std::endl;
      gettimeofday(&t1, NULL);
      write_mesh_selected(mesh, submshfile.format, submshfile.base);
      binary_write(nod.begin(), nod.end(), submshfile.base + NOD_EXT);
      binary_write(eidx.begin(), eidx.end(), submshfile.base + EIDX_EXT);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      break;
    }

    case OP_VOL:
    {
      // extract volume --------------------------------------------------------------
      mt_vector<std::string> coord_string_list;
      split_string(opts.coord, ',', coord_string_list);

      if(coord_string_list.size() != 6) {
        std::cerr << "Error! Wrong volume format! Aborting!" << std::endl;
        return;
      }

      int output_mode = 0;
      if(opts.mode.size())
        output_mode = atoi(opts.mode.c_str());
      check_condition(output_mode == 0 || output_mode == 1, "output_mode == [0|1]", "extract volume");

      float ref_x  = atof(coord_string_list[0].c_str());
      float ref_xd = atof(coord_string_list[1].c_str());
      float ref_y  = atof(coord_string_list[2].c_str());
      float ref_yd = atof(coord_string_list[3].c_str());
      float ref_z  = atof(coord_string_list[4].c_str());
      float ref_zd = atof(coord_string_list[5].c_str());

      mesh.e2n_dsp.resize(mesh.e2n_cnt.size());
      bucket_sort_offset(mesh.e2n_cnt, mesh.e2n_dsp);

      std::cout << "Extracting volume .." << std::endl;
      gettimeofday(&t1, NULL);
      mt_vector<mt_int> nod, eidx;

      // generate vector of elements to keep
      size_t nelem = mesh.e2n_cnt.size();
      mt_vector<bool> keep(nelem, false);
      #ifdef OPENMP
      #pragma omp parallel for schedule(guided, 100)
      #endif
      for(size_t idx=0; idx < nelem; idx++) {
        float x = 0, y = 0, z = 0;
        int esize = mesh.e2n_cnt[idx];

        // we use the center point of an element to check if
        // the element is inside the given box
        for(int i=0; i<esize; i++) {
          int vidx = mesh.e2n_con[mesh.e2n_dsp[idx] + i];
          x += mesh.xyz[vidx*3+0];
          y += mesh.xyz[vidx*3+1];
          z += mesh.xyz[vidx*3+2];
        }
        x /= esize;
        y /= esize;
        z /= esize;

        bool inside_x = ( (x >= ref_x) && (x <= (ref_x + ref_xd)) );
        bool inside_y = ( (y >= ref_y) && (y <= (ref_y + ref_yd)) );
        bool inside_z = ( (z >= ref_z) && (z <= (ref_z + ref_zd)) );

        keep[idx] = inside_x && inside_y && inside_z;
      }
      restrict_meshdata(keep, mesh, nod, eidx);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      mt_filename submshfile(opts.submsh_base, opts.ofmt);
      if(output_mode == 1) {
        std::string ofile = submshfile.base + VTX_EXT;
        std::cout << "Writing vertex file " << ofile << std::endl;
        write_vtx(nod, ofile, true);
      }
      else if(output_mode == 0){
        std::cout << "Writing mesh: " << submshfile.base << std::endl;
        gettimeofday(&t1, NULL);
        write_mesh_selected(mesh, submshfile.format, submshfile.base);
        binary_write(nod.begin(), nod.end(), submshfile.base + NOD_EXT);
        binary_write(eidx.begin(), eidx.end(), submshfile.base + EIDX_EXT);
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      }

      break;
    }
    case OP_UR:
    {
      // we extract unreachable elements -------------------------------------------------
      mt_vector<mt_meshdata> parts;
      mt_vector<mt_vector<mt_int> > dcmp_eidx;

      compute_full_mesh_connectivity(mesh, mshfile.base);
      nodal_connectivity_decomposition(mesh, parts, dcmp_eidx);

      std::cout << "Extracted " << parts.size() << " parts." << std::endl;
      if(parts.size() > 1) {
        for(size_t p = 0; p<parts.size(); p++)
        {
          mt_vector<mt_int> nod;
          reindex_nodes(parts[p], nod, false);

          mt_filename submshfile(opts.submsh_base, opts.ofmt);
          submshfile.base += ".part" + std::to_string(p);
          std::cout << "Writing mesh " << submshfile.base << " .." << std::endl;
          write_mesh_selected(parts[p], submshfile.format, submshfile.base);

          binary_write(nod.begin(), nod.end(), submshfile.base + NOD_EXT);
          binary_write(dcmp_eidx[p].begin(), dcmp_eidx[p].end(), submshfile.base + EIDX_EXT);
        }
      }
      else {
        std::cout << "Skipping mesh output .." << std::endl;
      }

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      break;
    }

    case OP_GRAD:
    {
      compute_full_mesh_connectivity(mesh, mshfile.base);

      // data_idx may be:
      // 0 : scalar dat file
      // 1 : scalar igb file
      // 2 : vec file
      // 3 : vec3 igb file
      short data_idx = -1;
      mt_vector<mt_real>            idat,     odat;
      mt_vector<mt_point<mt_real> > idat_vec, odat_vec;
      igb_header igb, igb_out, igb_out_vec;
      setup_data_format(opts.idat, opts.odat, data_idx, igb, igb_out);
      if(data_idx > 1) {
        std::cerr << "Error: Only scalar input data allowed! Aborting!" << std::endl;
        exit(1);
      }

      mt_filename outfile(opts.odat, "");
      bool nodal_input = true;
      bool nodal_output = true;

      if(opts.mode.size()) {
        int flg = atoi(opts.mode.c_str());
        if(flg == 0) nodal_output = true;
        if(flg == 1) nodal_output = false;
      }

      if(nodal_output)
        std::cout << "Output data is node based." << std::endl;
      else
        std::cout << "Output data is element based." << std::endl;

      switch(data_idx) {
        case 0:
        {
          read_vector_ascii(idat, opts.idat, true);

          if(idat.size() == mesh.e2n_cnt.size()) {
            std::cout << "Input data is element based." << std::endl;
            nodal_input = false;
          }
          else if(idat.size() == mesh.n2e_cnt.size()) {
            std::cout << "Input data is node based." << std::endl;
            nodal_input = true;
          }
          else {
            std::cerr << "Error: Data dimension does not fit either nodes or elems. Aborting!" << std::endl;
            exit(1);
          }

          compute_gradient(mesh, idat, nodal_input, nodal_output, odat_vec, odat);

          std::string outstr = outfile.base + ".grad" + VEC_EXT;
          write_vector_ascii(odat_vec, outstr);
          outstr = outfile.base + ".gradmag" + DAT_EXT;
          write_vector_ascii(odat, outstr);

          break;
        }

        case 1:
        {
          if(igb.v_x == int(mesh.e2n_cnt.size())) {
            std::cout << "Input data is element based." << std::endl;
            nodal_input = false;
          }
          else if(igb.v_x == int(mesh.n2e_cnt.size())) {
            std::cout << "Input data is node based." << std::endl;
            nodal_input = true;
          }
          else {
            std::cerr << "Error: Data dimension does not fit either nodes or elems. Aborting!" << std::endl;
            exit(1);
          }

          igb_out = igb; igb_out_vec = igb;
          igb_out.filename     = outfile.base + ".gradmag" + IGB_EXT;
          set_igb_header_datatype("double", igb_out);
          igb_out_vec.filename = outfile.base + ".grad" + IGB_EXT;
          set_igb_header_datatype("vec3d", igb_out_vec);
          write_igb_header(igb_out);
          write_igb_header(igb_out_vec);

          std::vector<double> rbuff;
          printf("Processing igb time-slices: %s to %s: \n", opts.idat.c_str(), opts.odat.c_str());

          for(int t=0; t<igb.v_t; t++) {
            printf("\rcurrent time-slice %d / %d .. ", t+1, int(igb.v_t));
            fflush(stdout);

            read_igb_slice(rbuff, igb);
            idat.assign(rbuff.begin(), rbuff.end());

            compute_gradient(mesh, idat, nodal_input, nodal_output, odat_vec, odat);

            rbuff.assign(odat.begin(), odat.end());
            write_igb_slice(rbuff, igb_out);

            rbuff.resize(odat_vec.size() * 3);
            for(size_t i=0; i<odat_vec.size(); i++) {
              rbuff[i*3+0] = odat_vec[i].x;
              rbuff[i*3+1] = odat_vec[i].y;
              rbuff[i*3+2] = odat_vec[i].z;
            }
            write_igb_slice(rbuff, igb_out_vec);
          }
          printf("\n");

          fclose(igb.fileptr);
          fclose(igb_out.fileptr);
          fclose(igb_out_vec.fileptr);
          break;
        }
        default: break;
      }
      break;
    }

    case OP_OVERLP:
    {
      mt_meshdata mesh1, mesh2;

      mt_filename msh1(opts.msh1_base, opts.ifmt);
      mt_filename msh2(opts.msh2_base, opts.ifmt);

      std::cout << "Reading mesh 1: " << msh1.base << std::endl;
      gettimeofday(&t1, NULL);
      read_mesh_selected(mesh1, msh1.format, msh1.base);
      bucket_sort_offset(mesh1.e2n_cnt, mesh1.e2n_dsp);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Reading mesh 2: " << msh2.base << std::endl;
      gettimeofday(&t1, NULL);
      read_mesh_selected(mesh2, msh2.format, msh2.base);
      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      std::cout << "Identifying overlaping elements .." << std::endl;
      gettimeofday(&t1, NULL);
      mt_vector<mt_real> m1_dat(mesh1.xyz.size() / 3, 0.0), m1_dat_ele;

      float radius = 2.0;
      if(opts.rad.size()) radius = atof(opts.rad.c_str());

      mt_real avrg_est = avrg_edgelength_estimate(mesh1, true);
      radius *= avrg_est;

      // we select mesh1 nodes in a radius around mesh2 nodes
      mt_vector<vec3r> pts;
      array_to_points(mesh1.xyz, pts);

      kdtree tree(10);
      tree.build_vertex_tree(pts);

      size_t npts_mesh2 = mesh2.xyz.size() / 3;
      MT_USET<mt_int> sel_nod;

      #ifdef OPENMP
      #pragma omp parallel
      #endif
      {
        MT_USET<mt_int> loc_set;
        mt_vector<mixed_tuple<mt_real, int>> found_vtx;

        #ifdef OPENMP
        #pragma omp for schedule(guided)
        #endif
        for(size_t nidx = 0; nidx < npts_mesh2; nidx++) {
          vec3r ref(mesh2.xyz.data() + nidx*3);
          tree.vertices_in_sphere(ref, radius, found_vtx);
          if(found_vtx.size()) {
            for(mixed_tuple<mt_real, int> & t : found_vtx)
              loc_set.insert(t.v2);
          }
        }

        #ifdef OPENMP
        #pragma omp critical
        #endif
        {
          sel_nod.insert(loc_set.begin(), loc_set.end());
        }
      }

      for(mt_int & idx : sel_nod) m1_dat[idx] = 100.0;
      nodeData_to_elemData(mesh1, m1_dat, m1_dat_ele);

      // now we mark all elements that have nodes inside mesh2
      mt_vector<bool> keep(mesh1.e2n_cnt.size());
      for(size_t eidx = 0; eidx < mesh1.e2n_cnt.size(); eidx++)
        keep[eidx] = m1_dat_ele[eidx] > 0;

      gettimeofday(&t2, NULL);
      std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

      bool output_complement = false;
      if(opts.mode.size()) {
        if(atoi(opts.mode.c_str()) == 0) output_complement = false;
        if(atoi(opts.mode.c_str()) == 1) output_complement = true;
      }

      if(output_complement) {
        mt_vector<mt_int> m1_eidx, comp_eidx;
        m1_eidx.reserve(mesh1.e2n_cnt.size() / 2);
        comp_eidx.reserve(mesh1.e2n_cnt.size() / 2);

        for (size_t eidx = 0; eidx < mesh1.e2n_cnt.size(); eidx++) {
          if (keep[eidx])
            m1_eidx.push_back(mt_int(eidx));
          else
            comp_eidx.push_back(mt_int(eidx));
        }

        mt_meshdata overlap, complement;

        mt_vector<mt_int> nod;
        extract_mesh(m1_eidx, mesh1, overlap);
        reindex_nodes(overlap, nod, true);

        // write first mesh
        mt_filename submesh(opts.submsh_base, opts.ofmt);
        std::cout << "Writing mesh: " << submesh.base << std::endl;
        gettimeofday(&t1, NULL);
        write_mesh_selected(overlap, submesh.format, submesh.base);
        binary_write(nod.begin(), nod.end(), submesh.base + NOD_EXT);
        binary_write(m1_eidx.begin(), m1_eidx.end(), submesh.base + EIDX_EXT);
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

        extract_mesh(comp_eidx, mesh1, complement);
        reindex_nodes(complement, nod, true);

        // write second mesh
        submesh.base += ".compl";
        std::cout << "Writing mesh: " << submesh.base << std::endl;
        gettimeofday(&t1, NULL);
        write_mesh_selected(complement, submesh.format, submesh.base);
        binary_write(nod.begin(), nod.end(), submesh.base + NOD_EXT);
        binary_write(comp_eidx.begin(), comp_eidx.end(), submesh.base + EIDX_EXT);
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      }
      else {
        mt_vector<mt_int> nod, eidx;
        restrict_meshdata(keep, mesh1, nod, eidx);

        // write mesh
        mt_filename submesh(opts.submsh_base, opts.ofmt);
        std::cout << "Writing mesh: " << submesh.base << std::endl;
        gettimeofday(&t1, NULL);
        write_mesh_selected(mesh1, submesh.format, submesh.base);
        binary_write(nod.begin(), nod.end(), submesh.base + NOD_EXT);
        binary_write(eidx.begin(), eidx.end(), submesh.base + EIDX_EXT);
        gettimeofday(&t2, NULL);
        std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;
      }

      break;
    }

    case OP_TAGS:
    {
      write_vector_ascii(mesh.etags, opts.odat);
      break;
    }

    default: break;
  }
}
