/**
* @file mmg_utils.h
* @brief MMG format specific IO utils.
* @author Aurel Neic
* @version
* @date 2016-12-13
*/


#ifndef _MMG_UTILS
#define _MMG_UTILS

void mmg_process_points(mt_meshdata & mesh,
                        char* buff, const int buffsize, FILE* fin);

void mmg_process_tet(mt_meshdata & mesh,
                     char* buff, const int buffsize, FILE* fin);
void mmg_process_hex(mt_meshdata & mesh,
                     char* buff, const int buffsize, FILE* fin);
void mmg_process_tri(mt_meshdata & mesh,
                     char* buff, const int buffsize, FILE* fin);
void mmg_process_quad(mt_meshdata & mesh,
                     char* buff, const int buffsize, FILE* fin);
void mmg_process(mt_meshdata & mesh,
                 char* buff, const int buffsize, FILE* fin);


void read_mmg_mesh(mt_meshdata & mesh, std::string filename);


void write_mmg_mesh(mt_meshdata & mesh, std::string filename);

#endif

