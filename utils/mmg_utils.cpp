/**
* @file mmg_utils.cpp
* @brief MMG format specific IO utils.
* @author Aurel Neic
* @version
* @date 2017-08-04
*/

#include "mt_utils_base.h"
#include "mmg_utils.h"
#include "io_utils.h"
#include "vtk_utils.h"

void mmg_process_points(mt_meshdata & mesh,
                        char* buff, const int buffsize, FILE* fin)
{
  unsigned long int npts = 0;
  char* ptr = fgets(buff, buffsize, fin);
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
}

void mmg_process_tet(mt_meshdata & mesh,
                     char* buff, const int buffsize, FILE* fin)
{
  unsigned long int nelem = 0;

  char* ptr = fgets(buff, buffsize, fin);
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
    if(ptr) sscanf(buff, "%lu %lu %lu %lu %d", vtx, vtx+1, vtx+2, vtx+3, &tag);
    con[0] = vtx[0];
    con[1] = vtx[1];
    con[2] = vtx[2];
    con[3] = vtx[3];
    mesh.etags[i] = tag;
    con += 4;
  }
}

void mmg_process_hex(mt_meshdata & mesh,
                     char* buff, const int buffsize, FILE* fin)
{
  unsigned long int nelem = 0;

  char* ptr = fgets(buff, buffsize, fin);
  if(ptr) sscanf(buff, "%lu", &nelem);

  mesh.e2n_cnt.assign(nelem, 8);
  mesh.etype.assign  (nelem, Hexa);

  mesh.etags.resize  (nelem);
  mesh.e2n_con.resize(nelem*8);

  unsigned long int vtx[8];
  int tag;
  mt_int* con = mesh.e2n_con.data();

  for (size_t i=0; i<nelem; i++ )
  {
    ptr = fgets(buff, buffsize, fin);
    if(ptr) sscanf(buff, "%lu %lu %lu %lu %lu %lu %lu %lu %d", vtx, vtx+1, vtx+2, vtx+3, vtx+4,vtx+5,vtx+6,vtx+7,&tag);
    con[0] = vtx[7];
    con[1] = vtx[6];
    con[2] = vtx[5];
    con[3] = vtx[4];
    con[4] = vtx[3];
    con[5] = vtx[0];
    con[6] = vtx[1];
    con[7] = vtx[2];

    mesh.etags[i] = tag;
    con += 8;
  }
}

void mmg_process_tri(mt_meshdata & mesh,
                     char* buff, const int buffsize, FILE* fin)
{
  // dont read triangles if tets have already been stored
  if(mesh.e2n_cnt.size()) return;

  unsigned long int nelem = 0;

  char* ptr = fgets(buff, buffsize, fin);
  if(ptr) sscanf(buff, "%lu", &nelem);

  mesh.e2n_cnt.assign(nelem, 3);
  mesh.etype.assign  (nelem, Tri);

  mesh.etags.resize  (nelem);
  mesh.e2n_con.resize(nelem*3);

  unsigned long int vtx[3];
  int tag;
  mt_int* con = mesh.e2n_con.data();

  for (size_t i=0; i<nelem; i++ )
  {
    ptr = fgets(buff, buffsize, fin);
    if(ptr) sscanf(buff, "%lu %lu %lu %d", vtx, vtx+1, vtx+2, &tag);
    con[0] = vtx[0];
    con[1] = vtx[1];
    con[2] = vtx[2];
    mesh.etags[i] = tag;
    con += 3;
  }
}

void mmg_process_quad(mt_meshdata & mesh,
                     char* buff, const int buffsize, FILE* fin)
{
  // dont read triangles if tets have already been stored
  if(mesh.e2n_cnt.size()) return;

  unsigned long int nelem = 0;

  char* ptr = fgets(buff, buffsize, fin);
  if(ptr) sscanf(buff, "%lu", &nelem);

  mesh.e2n_cnt.assign(nelem, 4);
  mesh.etype.assign  (nelem, Quad);

  mesh.etags.resize  (nelem);
  mesh.e2n_con.resize(nelem*4);

  unsigned long int vtx[4];
  int tag;
  mt_int* con = mesh.e2n_con.data();

  for (size_t i=0; i<nelem; i++ )
  {
    ptr = fgets(buff, buffsize, fin);
    if(ptr) sscanf(buff, "%lu %lu %lu %lu %d", vtx, vtx+1, vtx+2, vtx+3,&tag);
    con[0] = vtx[0];
    con[1] = vtx[1];
    con[2] = vtx[2];
    con[3] = vtx[3];
    mesh.etags[i] = tag;
    con += 4;
  }
}

void mmg_process(mt_meshdata & mesh,
                 char* buff, const int buffsize, FILE* fin)
{
  char keyword[256];
  char rest[2048];

  int n = sscanf(buff, "%s %s", keyword, rest);
  if(n < 1) return;

  // process the keyword
  if(strcmp(keyword, "Vertices") == 0) {
    mmg_process_points(mesh, buff, buffsize, fin);
  }
  else if(strcmp(keyword, "Tetrahedra") == 0) {
    mmg_process_tet(mesh, buff, buffsize, fin);
  }
  else if(strcmp(keyword, "Triangles") == 0) {
    mmg_process_tri(mesh, buff, buffsize, fin);
  }
  else if(strcmp(keyword, "Quadrilaterals") == 0) {
    mmg_process_quad(mesh, buff, buffsize, fin);
  }
  else if(strcmp(keyword, "Hexahedra") == 0) {
    mmg_process_hex(mesh, buff, buffsize, fin);
  }
}


void read_mmg_mesh(mt_meshdata & mesh, std::string filename)
{
  FILE* fd = fopen(filename.c_str(), MT_FOPEN_READ);
  if(fd == NULL) treat_file_open_error(filename, errno);

  const int buffsize = 2048;
  char buffer[buffsize];

  // start with processing
  char* ptr = next_line(buffer, buffsize, fd);
  while(ptr != NULL) {
    mmg_process(mesh, buffer, buffsize, fd);
    ptr = next_line(buffer, buffsize, fd);
  }

  mt_int min = *(std::min_element(mesh.e2n_con.begin(), mesh.e2n_con.end()));

  if(min > 0) {
    for(size_t i=0; i < mesh.e2n_con.size(); i++)
      mesh.e2n_con[i] -= 1;
  }
  generate_default_fib(mesh.e2n_cnt.size(), mesh.lon);

  fclose(fd);
}


void write_mmg_mesh(mt_meshdata & mesh, std::string filename)
{
  FILE* fd = fopen(filename.c_str(), MT_FOPEN_WRITE);
  if(fd == NULL) treat_file_open_error(filename, errno);

  fprintf(fd, "MeshVersionFormatted 2\n\n");
  fprintf(fd, "Dimension 3\n\n");

  unsigned long int npts = mesh.xyz.size() / 3;
  fprintf(fd, "Vertices\n");
  fprintf(fd, "%lu\n", npts);

  for(size_t i=0; i<npts; i++)
    fprintf(fd, "%f %f %f 0\n", mesh.xyz[i*3+0], mesh.xyz[i*3+1], mesh.xyz[i*3+2]);

  fprintf(fd, "\n");

  unsigned long int nele = mesh.e2n_cnt.size();

  // we assume only one element type -- either Tet or Tri -- exists in the mesh
  const char* etype = (mesh.etype[0] == Tetra) ? "Tetrahedra" : "Triangles";

  fprintf(fd, "%s\n", etype);
  fprintf(fd, "%lu\n", nele);

  if(mesh.etype[0] == Tetra) {
    for(size_t i=0; i<nele; i++)
      fprintf(fd, "%lu %lu %lu %lu %d\n", (long unsigned int) mesh.e2n_con[i*4+0] + 1,
                                          (long unsigned int) mesh.e2n_con[i*4+1] + 1,
                                          (long unsigned int) mesh.e2n_con[i*4+2] + 1,
                                          (long unsigned int) mesh.e2n_con[i*4+3] + 1,
                                          (int) mesh.etags[i]);
  }
  else {
    for(size_t i=0; i<nele; i++)
      fprintf(fd, "%lu %lu %lu %d\n", (long unsigned int) mesh.e2n_con[i*3+0] + 1,
                                      (long unsigned int) mesh.e2n_con[i*3+1] + 1,
                                      (long unsigned int) mesh.e2n_con[i*3+2] + 1,
                                      (int) mesh.etags[i]);
  }

  fprintf(fd, "\nEnd\n\n");

  fclose(fd);
}
