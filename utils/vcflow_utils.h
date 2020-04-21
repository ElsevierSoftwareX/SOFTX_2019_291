
/**
* @file vcflow_utils.h
* @brief VCFlow format specific IO utils.
* @author Elias Karabelas
* @version
* @date 2018-08-22
*/


#ifndef _VCFLOW_UTILS
#define _VCFLOW_UTILS

void read_vc_flow_points(mt_vector<mt_real> & xyz, std::string file);

void write_vc_flow_points(mt_vector<mt_real> & xyz, std::string file);

void read_vc_flow_elements(mt_meshdata & mesh, std::string file);

void write_vc_flow_elements(mt_meshdata & mesh, const std::string file);

void write_vc_flow_adjacency(mt_meshdata & mesh, const std::string file);

#endif