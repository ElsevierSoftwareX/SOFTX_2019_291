#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "mt_modes_base.h"


struct retag_options{
  mt_filename msh;
  mt_filename outmsh;
  std::string ifmt;
  std::string ofmt;
  std::string con;
  std::string neigh;
  std::string verbose;
};


#define CON_DFLT 50 // face connected neighbors within single tag region
static const std::string neigh_par = "-neigh=";
static const std::string verbose_par = "-verbose=";
// status of processed elements
static const mt_int IS_OTHER = -1;  // element belongs to different tag region
static const mt_int IS_PROCESS = 1; // element was source for neighbor search
static const mt_int IS_NEIGH = 4;   // element found to be a neighbor
static const mt_int IS_DONE = 5;    // combination of both above


void print_retag_help()
{
  fprintf(stderr, "retag_adv: identify and re-tag elements that are not connected to enough elements of the same tag.\n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t  (input) path to basename of the mesh.\n", mesh_par.c_str());
  fprintf(stderr, "%s<path>\t  (output) path to basename of the ouput mesh.\n", outmesh_par.c_str());
  fprintf(stderr, "%s<int>\t  (optional) Number of neighbours threshold. Elements with less neighbours are isolated. Default is %d.\n",
      neigh_par.c_str(), CON_DFLT);
  fprintf(stderr, "%s<format>\t  (optional) mesh input format. may be: %s\n", inp_format_par.c_str(), input_formats.c_str());
  fprintf(stderr, "%s<format>\t  (optional) mesh output format. may be: %s\n", out_format_par.c_str(), output_formats.c_str());
  fprintf(stderr, "%s<format> (optional) increate commandline verbosity. Default: %s\n", verbose_par.c_str(), "0");
  fprintf(stderr, "\n");
}

int retag_parse_options(int argc, char** argv, struct retag_options & opts)
{
  if(argc < 2) {
    print_retag_help();
    return 1;
  }

  std::string msh_base, outmsh_base, ifmt, ofmt;

  // parse all retag parameters ------------------------------------------------
  for(int i=1; i<argc; i++){
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, mesh_par, msh_base);
    if(!match) match = parse_param(param, outmesh_par, outmsh_base);
    if(!match) match = parse_param(param, out_format_par, ofmt);
    if(!match) match = parse_param(param, inp_format_par, ifmt);
    if(!match) match = parse_param(param, neigh_par, opts.neigh);
    if(!match) match = parse_param(param, verbose_par, opts.verbose);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 3;
    }
  }
  opts.msh.assign(msh_base, ifmt);
  opts.outmsh.assign(outmsh_base, ofmt);

  // check if all relevant parameters have been set ----------------------------
  if( !(opts.msh.isSet() && opts.outmsh.isSet()) )
  {
    std::cerr << "Error: Insufficient parameters provided." << std::endl;
    print_retag_help();
    return 2;
  }

  return 0;
}

void postprocess_loop(mt_meshdata &mesh, mt_meshgraph &submeshg, mt_vector<mt_int> &elem_task,
                      size_t NNEIGH, size_t num_con, MT_USET<triple<mt_int>> &logbook, short verbose)
{
  std::set<mt_int> geom;
  std::set<mt_int> face_neigh;
  std::set<mt_int> geom_tags; // tags of all elements within geom (island)
                              // must be only one, otherwise something is wrong
  std::vector<mt_int> tag_candidates;

  if( num_con == submeshg.eidx.size() ){
    // all elements are connected with faces
    // let's continue with next tag region
    for(size_t eidx = 0; eidx < elem_task.size(); eidx++) {
      if( elem_task[eidx] == IS_DONE ) {
        elem_task[eidx] = IS_OTHER;
      }
    }
    return;
  }

  // some unprocessed elements seem to be remaining
  // continue processing elements of current tag group

  if( num_con > NNEIGH ) {
    // shadow processed elements of current loop
    for(size_t eidx = 0; eidx < elem_task.size(); eidx++) {
      if( elem_task[eidx] == IS_DONE ) {
        elem_task[eidx] = IS_OTHER;
      }
    }
    return;
  }

  /****************************************************************************
    Found node and/or edge connected elements - treatment necessary
  ****************************************************************************/

  // At first, collect global element indices of this island.
  // Did not want to create this collection in every outer loop!
  geom_tags.clear();
  for(size_t eidx = 0; eidx < elem_task.size(); eidx++) {
    if( elem_task[eidx] == IS_DONE ) {
      geom.insert(eidx);

      // in addition shadow these elements for the next outer loop
      elem_task[eidx] = IS_OTHER;

      // collect element tag(s) of all island elements
      geom_tags.insert(mesh.etags[eidx]);
    }
  }

  // re-assign these elements to non-zero (tagwise) neighboring element regions
  // current strategy: Collect all neighbor tags different from current tag.
  //                   Highest count on the same tag will be taken.
  // Note: On interfaces, it could be more reasonable to query the tag<F9>
  //       neighborhood element-wise!
  // Note: Taking again only face neighbors into consideration.
  //       Querying any node-connected neighbor would be as good (or better).
  face_neigh.clear();
  for(auto eidx = geom.begin(); eidx != geom.end(); ++eidx)
  {
    // determine forming nodes of element eidx
    mt_int v0 = mesh.e2n_con[mesh.e2n_dsp[*eidx] +0];
    mt_int v1 = mesh.e2n_con[mesh.e2n_dsp[*eidx] +1];
    mt_int v2 = mesh.e2n_con[mesh.e2n_dsp[*eidx] +2];
    mt_int v3 = mesh.e2n_con[mesh.e2n_dsp[*eidx] +3];

    // add elements which share a face with eidx to the set face_neigh
    elements_with_face(mesh, v0, v1, v3, face_neigh);
    elements_with_face(mesh, v0, v1, v2, face_neigh);
    elements_with_face(mesh, v0, v3, v2, face_neigh);
    elements_with_face(mesh, v1, v2, v3, face_neigh);
  }

  // collect element tags of neighbors
  tag_candidates.clear();

  for(auto feidx = face_neigh.begin(); feidx != face_neigh.end(); ++feidx) {
    if (mesh.etags[*feidx] == *geom_tags.begin() )
      // exclude tag of island
      continue;
    /*
    if (mesh.etags[*feidx] == 0 )
      // refuse to set bath as new tag
      continue;
    */
    tag_candidates.push_back(mesh.etags[*feidx]);
  }

  // query for highest frequency of neighboring element labels
  if( tag_candidates.size() ) {
    mt_int ntag = 0;
    mt_int ntag_freq = 0;
    for(auto it = tag_candidates.begin(); it != tag_candidates.end(); ++it) {
      mt_int freq = std::count(tag_candidates.begin(), tag_candidates.end(), *it);
      if( freq > ntag_freq ) {
        ntag = *it;
        ntag_freq = freq;
      }
    }
    // apply new label
    for(auto eidx = geom.begin(); eidx != geom.end(); ++eidx) {
      triple<mt_int> trp = {*eidx, *geom_tags.begin(), ntag};
      logbook.insert(trp);

      mesh.etags[*eidx] = ntag;
    }
  } else if (verbose) {
    std::cout << std::endl << "Retagging following element(s) failed for unknown reasons: " << std::endl;
    for(auto eidx = geom.begin(); eidx != geom.end(); ++eidx) {
      std::cout << *eidx << std::endl << std::flush;
      // add them to the logbook increasing tag by 1000
      triple<mt_int> trp = {*eidx, *geom_tags.begin(), *geom_tags.begin()+1000};
      logbook.insert(trp);
    }
  }
}


void write_logbook(std::string file, MT_USET<triple<mt_int>> &logbook, size_t numele)
{
  mt_vector<mt_int> data; // full size data vector
  std::string fullfile = file + ".retaglog.dat";
  FILE* dat_file = fopen(fullfile.c_str(), MT_FOPEN_WRITE);

  if(dat_file == NULL) treat_file_open_error(fullfile, errno);

  // assemble full length element vector with "new tag" as value
  data.resize(numele,-1);
  for(auto it = logbook.begin(); it != logbook.end(); ++it) {
    triple<mt_int> trp = *it;
    data[trp.v1] = trp.v3;
  }

  // write header and data values
  //fprintf(dat_file, "%lu\n", numele);
  for(size_t i=0; i<numele; i++) {
    fprintf(dat_file, "%ld\n", data[i]);
  }
  fclose(dat_file);
}



int main(int argc, char** argv)
{
  struct timeval t1, t2;
  struct retag_options opts;
  struct mt_meshdata mesh;
  MT_USET<triple<mt_int>> logbook;

  int ret = retag_parse_options(argc, argv, opts);
  if (ret != 0) return 1;

  std::cout << "Reading mesh: " << opts.msh.base << std::endl;
  gettimeofday(&t1, NULL);
  read_mesh_selected(mesh, opts.msh.format, opts.msh.base);
  compute_full_mesh_connectivity(mesh);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;

  size_t NNEIGH = opts.neigh.size() ? atoi(opts.neigh.c_str()) : CON_DFLT;
  short  verbose = opts.verbose.size() ? atoi(opts.verbose.c_str()) : 0;

  // get unique tag list
  mt_vector<mt_int> tags = mesh.etags;
  binary_sort(tags); unique_resize(tags);


  for(size_t i=0; i<tags.size(); i++)
  {
    std::set<mt_int> face_neigh;
    MT_USET<mt_int>  tagset;
    mt_meshgraph     submeshg;


    mt_int etag = tags[i]; tagset.clear(); tagset.insert(etag);
    extract_tagged_meshgraph(tagset, mesh, submeshg);


    // keep track of elems done with face neighbor search
    mt_vector<mt_int> elem_task;
    elem_task.resize(mesh.e2n_cnt.size(), -1);

    // mark elements of current tag region as not yet processed
    for(size_t j = 0; j < submeshg.eidx.size(); j++) {
      elem_task[submeshg.eidx[j]] = 0;
    }

    //std::cout << "Procossing tag region " << etag << std::endl << std::flush;

    // collect face connected elements: limit search range
    size_t EIDX_BEG = submeshg.eidx[0];
    size_t loop_cnt = 0;

    // mark first element already as a neighbor
    // every next search must be done with an element already labelled as a neighbor
    elem_task[EIDX_BEG] |= IS_NEIGH;

    do {

      for(size_t eidx = EIDX_BEG; eidx < elem_task.size(); eidx++)
      {
        if( elem_task[eidx] == IS_OTHER )   continue;
        if( elem_task[eidx] == IS_PROCESS ) continue;
        if( elem_task[eidx] == IS_DONE )    continue;
        if(!elem_task[eidx] )               continue;

        // eidx is an neighbor to a previous element
        // but was not yet used as source in a neighbor search

        // determine forming nodes of element eidx
        mt_int v0 = mesh.e2n_con[mesh.e2n_dsp[eidx] +0];
        mt_int v1 = mesh.e2n_con[mesh.e2n_dsp[eidx] +1];
        mt_int v2 = mesh.e2n_con[mesh.e2n_dsp[eidx] +2];
        mt_int v3 = mesh.e2n_con[mesh.e2n_dsp[eidx] +3];

        // add elements which share a face with eidx to the set face_neigh
        face_neigh.clear();
        elements_with_face(mesh, v0, v1, v3, face_neigh);
        elements_with_face(mesh, v0, v1, v2, face_neigh);
        elements_with_face(mesh, v0, v3, v2, face_neigh);
        elements_with_face(mesh, v1, v2, v3, face_neigh);


        for(auto it = face_neigh.begin(); it != face_neigh.end(); ++it) {
          switch( elem_task[*it] ) {
            case IS_OTHER : break;
            case IS_DONE  : break;
            default : elem_task[*it] |= IS_NEIGH;
          }
        }
        // set status of current element to "processed"
        elem_task[eidx] |= IS_PROCESS;
      }
      loop_cnt++;


      // analyze results of current loop and reset variables
      size_t num_con = 0;
      size_t neigh_cnt = 0;
      size_t unprocessed_cnt = 0;

      for(size_t eidx = 0; eidx < elem_task.size(); eidx++) {
        switch( elem_task[eidx] ) {
        case 0 :        unprocessed_cnt++; break;
        case IS_OTHER : break;
        case IS_DONE :  num_con++; break; // collected amount of connected elements
        case IS_NEIGH : neigh_cnt++;
                        if( eidx < EIDX_BEG ) EIDX_BEG = eidx;
                        break;
        default : std::cout << "Unhandled condition found (analyze results)!" << std::endl;
        }
      }
      std::cout << "\r" << "Collecting face connected elements [tag: " << etag << "]: ";
      std::cout << num_con << std::flush;

      if( !neigh_cnt ) {
        postprocess_loop(mesh, submeshg, elem_task, NNEIGH, num_con, logbook, verbose);
        std::cout << std::endl << std::flush;

        if( unprocessed_cnt ) {
          // There is no neighbor left. Need to move on to the next unprocessed element.
          for(size_t eidx = 0; eidx < elem_task.size(); eidx++) {
            if(! elem_task[eidx]){
              EIDX_BEG = eidx;
              elem_task[eidx] |= IS_NEIGH;
              break;
            }
          }
        } else {
          // we are done with current tag
          break;
        }
      }
    } while( true );
  }

  if( logbook.size() ) {
    std::cout << std::string(40 + log10(logbook.size())+1,'=') << std::endl;
    std::cout << "=== Identified " << logbook.size() << " problematic elements ===" << std::endl;
    std::cout << std::string(40 + log10(logbook.size())+1,'=') << std::endl;

    std::cout << "Writing logbook: " << opts.outmsh.base+".retaglog.dat" << std::endl;
    write_logbook(opts.outmsh.base, logbook, mesh.etags.size());
  } else {
    std::cout << "No tag changes were made!" << std::endl;
  }
  std::cout << "Writing mesh: " << opts.outmsh.base << std::endl;
  gettimeofday(&t1, NULL);
  write_mesh_selected(mesh, opts.outmsh.format, opts.outmsh.base);
  gettimeofday(&t2, NULL);
  std::cout << "Done in " << (float)timediff_sec(t1, t2) << " sec" << std::endl;


  return 0;
}
