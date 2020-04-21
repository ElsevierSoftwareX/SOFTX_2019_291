/**
* @file main.cpp
* @brief File holding the meshtool main function.
* @author Aurel Neic
* @version 
* @date 2017-09-12
*/

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <thread>

#include "mt_modes.h"


int main(int argc, char** argv)
{
  #ifdef OPENMP
  unsigned int maxThreads = std::thread::hardware_concurrency();

  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    if (tid == 0)
    {
      int nthreads = omp_get_num_threads();
      printf("OpenMP thread utilization: %d / %d threads.\n", nthreads, maxThreads);
    }
  }
  #endif

  if(argc < 2) {
    print_usage(argv[0]);
    return 1;
  }

  std::string mode = argv[1];

  if(mode.compare("insert") == 0)
    insert_mode(argc, argv);
  else if(mode.compare("extract") == 0)
    extract_mode(argc, argv);
  else if(mode.compare("convert") == 0)
    convert_mode(argc, argv);
  else if(mode.compare("collect") == 0)
    collect_mode(argc, argv);
  else if(mode.compare("split") == 0)
    split_mode(argc, argv);
  else if(mode.compare("map") == 0)
    map_mode(argc, argv);
  else if(mode.compare("merge") == 0)
    merge_mode(argc, argv);
  else if(mode.compare("generate") == 0)
    generate_mode(argc, argv);
  else if(mode.compare("clean") == 0)
    clean_mode(argc, argv);
  else if(mode.compare("smooth") == 0)
    smooth_mode(argc, argv);
  else if(mode.compare("query") == 0)
    query_mode(argc, argv);
  else if(mode.compare("resample") == 0)
    resample_mode(argc, argv);
  else if(mode.compare("restore") == 0)
    restore_mode(argc, argv);
  else if(mode.compare("itk") == 0)
    itk_mode(argc, argv);
  else if(mode.compare("interpolate") == 0)
    interpolate_mode(argc, argv);
  else if(mode.compare("reindex") == 0)
    reindex_mode(argc, argv);
  else if(mode.compare("transform") == 0)
    transform_mode(argc, argv);
  else if(mode.compare("help") == 0)
    print_general_help(argv[0]);
  else
    print_usage(argv[0]);

  return 0;
}
