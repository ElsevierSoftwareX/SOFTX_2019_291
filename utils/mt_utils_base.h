/**
* @file mt_utils_base.h
* @brief The meshtool utilities base header.
* @author Aurel Neic
* @version
* @date 2017-09-12
*/
#ifndef _MT_UTILS_BASE
#define _MT_UTILS_BASE

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <sys/time.h>

#include <iostream>
#include <cstring>
#include <fstream>
#include <list>
#include <set>
#include <map>

#include "data_structs.h"
#include "sort.h"

#include "hashmap.hpp"
#define MT_MAP hashmap::unordered_map
#define MT_USET hashmap::unordered_set

#define MT_PI 3.14159265

#ifdef WINBUILD
#define drand48() (((double) rand())/RAND_MAX) 
#endif

#define POW2(A) ((A)*(A))
#define POW5(A) (POW2(A)*POW2(A)*(A))
#define POW10(A) (POW5(A)*POW5(A))


#include "progress.hpp"

#ifdef SILENT_PROGRESS
#define PROGRESS progress_silent
#else
#define PROGRESS progress_bar_eta
#endif 

#include "check_utils.h"
#endif
