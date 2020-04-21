/**
* @file stellar_utils.cpp
* @brief Stellar format specific IO utils.
* @author Aurel Neic
* @version 
* @date 2017-08-04
*/

#include "mt_utils_base.h"
#include "io_utils.h"
#include "stellar_utils.h"

void read_stellar_points(mt_vector<mt_real> & xyz, std::set<mt_int> & surf_vtx, bool & zero_indexing, std::string filename)
{
  FILE* fd = fopen(filename.c_str(), MT_FOPEN_READ);
  if(fd == NULL) treat_file_open_error(filename, errno);

  unsigned long int numpoints = 0, dummy;
  int dim, attr, boundary;

  const int bufsize = 2048;
  char buffer[bufsize], *ptr;
  double pt[3];

  ptr = fgets(buffer, bufsize, fd);
  if(ptr != NULL) sscanf(buffer, "%lu %d %d %d", &numpoints, &dim, &attr, &boundary);

  xyz.resize(numpoints*3);

  for(size_t i=0; i<numpoints; i++)
  {
    ptr = fgets(buffer, bufsize, fd);
    if(ptr != NULL) sscanf(buffer, "%lu %lf %lf %lf %d %d", &dummy, pt, pt+1, pt+2, &dim, &boundary);

    // test indexing
    if(!i)
      zero_indexing = (dummy == 0);

    xyz[i*3+0] = pt[0];
    xyz[i*3+1] = pt[1];
    xyz[i*3+2] = pt[2];
    // we assume attr is at max 1, usually 0
    // if attr is 0, then the boundary info was the forth value
    // in the line and we have to do a copy
    if(!attr) boundary = dim;
    if(boundary) surf_vtx.insert((mt_int)i);
  }

  fclose(fd);
}

void write_stellar_points(mt_vector<mt_real> & xyz, std::set<mt_int> & surf_vtx, std::string filename)
{
  FILE* fd = fopen(filename.c_str(), MT_FOPEN_WRITE);
  if(fd == NULL) treat_file_open_error(filename, errno);

  double pt[3];

  size_t numpts = xyz.size()/3;
  int boundary = surf_vtx.size() > 0 ? 1 : 0;

  fprintf(fd, "%lu 3 0 %d\n", numpts, boundary);

  for(size_t i=0; i<numpts; i++)
  {
    pt[0] = xyz[i*3+0];
    pt[1] = xyz[i*3+1];
    pt[2] = xyz[i*3+2];

    boundary = surf_vtx.count(i);

    fprintf(fd, "%lu %lf %lf %lf %d\n", i+1, pt[0], pt[1], pt[2], boundary);
  }

  fclose(fd);
}

void read_stellar_elems(mt_meshdata & mesh, std::string filename)
{
  FILE* fd = fopen(filename.c_str(), MT_FOPEN_READ);
  if(fd == NULL) treat_file_open_error(filename, errno);

  unsigned long int numelem = 0, dummy;
  int esize, attr;

  const int bufsize = 2048;
  char buffer[bufsize], *ptr;
  unsigned long int vtx[4];

  ptr = fgets(buffer, bufsize, fd);
  if(ptr != NULL) sscanf(buffer, "%lu %d %d", &numelem, &esize, &attr);

  mesh.e2n_cnt.assign(numelem, 4);
  mesh.e2n_con.resize(numelem *4);

  mesh.etags.assign(numelem, mt_int(0));
  mesh.etype.assign(numelem, Tetra);

  for(size_t i=0; i<numelem; i++)
  {
    ptr = fgets(buffer, bufsize, fd);
    if(ptr != NULL) sscanf(buffer, "%lu %lu %lu %lu %lu", &dummy, vtx, vtx+1, vtx+2, vtx+3);
    #if 0
    mesh.e2n_con[i*4+0] = vtx[0];
    mesh.e2n_con[i*4+1] = vtx[1];
    mesh.e2n_con[i*4+2] = vtx[2];
    mesh.e2n_con[i*4+3] = vtx[3];
    #else
    mesh.e2n_con[i*4+0] = vtx[1];
    mesh.e2n_con[i*4+1] = vtx[2];
    mesh.e2n_con[i*4+2] = vtx[3];
    mesh.e2n_con[i*4+3] = vtx[0];
    #endif
  }

  fclose(fd);
}

void write_stellar_elems(mt_meshdata & mesh, std::string filename)
{
  FILE* fd = fopen(filename.c_str(), MT_FOPEN_WRITE);
  if(fd == NULL) treat_file_open_error(filename, errno);

  size_t numelem = mesh.e2n_cnt.size();
  unsigned long int vtx[4];

  fprintf(fd, "%lu 4 0\n", numelem);

  for(size_t i=0; i<numelem; i++)
  {
    vtx[0] = mesh.e2n_con[i*4+3] + 1;
    vtx[1] = mesh.e2n_con[i*4+0] + 1;
    vtx[2] = mesh.e2n_con[i*4+1] + 1;
    vtx[3] = mesh.e2n_con[i*4+2] + 1;
    fprintf(fd, "%lu %lu %lu %lu %lu\n", i+1, vtx[0], vtx[1], vtx[2], vtx[3]);
  }

  fclose(fd);
}


void read_stellar_mesh(mt_meshdata & mesh, std::set<mt_int> & surf_vtx, std::string basename)
{
  std::string nodefile = basename + "node";
  std::string elemfile = basename + "ele";

  bool zero_indexing = false;

  read_stellar_points(mesh.xyz, surf_vtx, zero_indexing, nodefile);
  read_stellar_elems(mesh, elemfile);
  if(!zero_indexing) {
    for(size_t i=0; i<mesh.e2n_con.size(); i++) mesh.e2n_con[i] -= 1;
  }
  generate_default_fib(mesh.e2n_cnt.size(), mesh.lon);
}


void write_stellar_mesh(mt_meshdata & mesh, std::set<mt_int> & surf_vtx, std::string basename)
{
  std::string nodefile = basename + "node";
  std::string elemfile = basename + "ele";

  write_stellar_points(mesh.xyz, surf_vtx, nodefile);
  write_stellar_elems(mesh, elemfile);
}

