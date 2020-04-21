/**
* @file mt_modes.h
* @brief The header file of the meshtool modes.
* @author Aurel Neic
* @version 
* @date 2017-09-12
*/

#ifndef _MT_MODES
#define _MT_MODES

#include "mt_modes_base.h"

void split_mode(int argc, char** argv);
void clean_mode(int argc, char** argv);
void convert_mode(int argc, char** argv);
void collect_mode(int argc, char** argv);
void extract_mode(int argc, char** argv);
void generate_mode(int argc, char** argv);
void insert_mode(int argc, char** argv);
void itk_mode(int argc, char** argv);
void map_mode(int argc, char** argv);
void merge_mode(int argc, char** argv);
void query_mode(int argc, char** argv);
void resample_mode(int argc, char** argv);
void restore_mode(int argc, char** argv);
void smooth_mode(int argc, char** argv);
void reindex_mode(int argc, char** argv);
void interpolate_mode(int argc, char** argv);
void transform_mode(int argc, char** argv);
#endif

