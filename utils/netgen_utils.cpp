/**
* @file netgen_utils.cpp
* @brief Netgen format specific IO utils.
* @author Elias Karabelas
* @version
* @date 2019-04-03
*/

#include "mt_utils_base.h"
#include "netgen_utils.h"
#include "io_utils.h"
#include "vtk_utils.h"

void netgen_process(mt_meshdata & mesh,
                     char* buff, const int buffsize, FILE* fin)
{

  unsigned long int nelem = 0;
  unsigned long int npts = 0;
  char* ptr = fgets(buff, buffsize, fin);

  //start with the points
  if(ptr) sscanf(buff, "%lu", &npts);
  mesh.xyz.resize(npts * 3);
  float pts[3];
  mt_real* p = mesh.xyz.data();
  for (size_t i=0; i<npts; i++)
  {
    ptr = fgets(buff, buffsize, fin);
    if(ptr) sscanf(buff, "%f %f %f", pts, pts+1, pts+2);
    p[0] = pts[0], p[1] = pts[1], p[2] = pts[2];
    p += 3;
  }
  //next the elements
  if(ptr) sscanf(buff, "%lu", &nelem);

  mesh.e2n_cnt.assign(nelem, 4);
  mesh.etype.assign  (nelem, Tetra);

  mesh.etags.resize  (nelem);
  mesh.e2n_con.resize(nelem*4);

  unsigned long int vtx[4];
  int tag;
  mt_int* con = mesh.e2n_con.data();

  for (size_t i=0; i<nelem; i++ )
  {
    ptr = fgets(buff, buffsize, fin);
    if(ptr) sscanf(buff, "%d %lu %lu %lu %lu", &tag, vtx, vtx+1, vtx+2, vtx+3);
    con[0] = vtx[0]-1;
    con[1] = vtx[1]-1;
    con[2] = vtx[2]-1;
    con[3] = vtx[3]-1;
    mesh.etags[i] = tag;
    con += 4;
  }
}

void read_netgen_mesh(mt_meshdata & mesh, std::string filename)
{
  FILE* fd = fopen(filename.c_str(), MT_FOPEN_READ);
  if(fd == NULL) treat_file_open_error(filename, errno);

  const int buffsize = 2048;
  char buffer[buffsize];

  // start with processing
  netgen_process(mesh, buffer, buffsize, fd);

  mt_int min = *(std::min_element(mesh.e2n_con.begin(), mesh.e2n_con.end()));

  if(min > 0) {
    for(size_t i=0; i < mesh.e2n_con.size(); i++)
      mesh.e2n_con[i] -= 1;
  }
  generate_default_fib(mesh.e2n_cnt.size(), mesh.lon);

  fclose(fd);
}


void write_netgen_mesh(mt_meshdata & mesh, std::string filename)
{
  FILE* fd = fopen(filename.c_str(), MT_FOPEN_WRITE);
  if(fd == NULL) treat_file_open_error(filename, errno);

  unsigned long int npts = mesh.xyz.size() / 3;
  fprintf(fd, "%lu\n", npts);

  for(size_t i=0; i<npts; i++)
    fprintf(fd, "%f %f %f\n", mesh.xyz[i*3+0], mesh.xyz[i*3+1], mesh.xyz[i*3+2]);

  fprintf(fd, "\n");

  unsigned long int nele = mesh.e2n_cnt.size();

  fprintf(fd, "%lu\n", nele);

  if(mesh.etype[0] == Tetra) {
    for(size_t i=0; i<nele; i++)
      fprintf(fd, "%d %lu %lu %lu %lu\n", (int) mesh.etags[i],
                                          (long unsigned int) mesh.e2n_con[i*4+0] + 1,
                                          (long unsigned int) mesh.e2n_con[i*4+1] + 1,
                                          (long unsigned int) mesh.e2n_con[i*4+2] + 1,
                                          (long unsigned int) mesh.e2n_con[i*4+3] + 1);
  }
  fclose(fd);
}
