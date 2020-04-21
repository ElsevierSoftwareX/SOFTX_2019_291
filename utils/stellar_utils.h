/**
* @file stellar_utils.h
* @brief IO functions for the Stellar mesh improvement tool.
* @author Aurel Neic
* @version
* @date 2017-02-13
*/

#ifndef _STELLAR_UTILS
#define _STELLAR_UTILS


void read_stellar_points(mt_vector<mt_real> & xyz, std::set<mt_int> & surf_vtx, bool & zero_indexing, std::string filename);

void write_stellar_points(mt_vector<mt_real> & xyz, std::set<mt_int> & surf_vtx, std::string filename);

void read_stellar_elems(mt_meshdata & mesh, std::string filename);

void write_stellar_elems(mt_meshdata & mesh, std::string filename);

void read_stellar_mesh(mt_meshdata & mesh, std::set<mt_int> & surf_vtx, std::string basename);

void write_stellar_mesh(mt_meshdata & mesh, std::set<mt_int> & surf_vtx, std::string basename);

#endif

