/**
* @file netgen_utils.h
* @brief Netgen neutral format specific IO utils.
* @author Elias Karabelas
* @version
* @date 2019-04-03
*/


#ifndef _NETGEN_UTILS
#define _NETGEN_UTILS

void netgen_process(mt_meshdata & mesh,
                 char* buff, const int buffsize, FILE* fin);


void read_netgen_mesh(mt_meshdata & mesh, std::string filename);


void write_netgen_mesh(mt_meshdata & mesh, std::string filename);

#endif

