/**
* @file vcflow_utils.cpp
* @brief VCFlow format specific IO utils.
* @author Elias Karabelas
* @version
* @date 2018-08-22
*/

#include "mt_utils_base.h"
#include "io_utils.h"
#include "vcflow_utils.h"
#include "topology_utils.h"

void read_vc_flow_points(mt_vector<mt_real> & xyz, std::string file)
{
  FILE* pts_file = fopen(file.c_str(), MT_FOPEN_READ);
  if(pts_file == NULL) treat_file_open_error(file, errno);

  int numpts;
  fread(&numpts, sizeof(int), 1, pts_file);

  xyz.resize(numpts*3);
  mt_real* wp = xyz.data();
  double pts[3];
  PROGRESS<int> progress(numpts, "Reading points in VCFlow format: ");

  for(int i=0; i<numpts; i++)
  {
    progress.next();

    fread(pts, sizeof(double), 3, pts_file);
    wp[0] = (mt_real)(pts[0]); wp[1] = (mt_real)(pts[1]); wp[2] = (mt_real)(pts[2]);
    wp += 3;
  }
  progress.finish();

  fclose(pts_file);
}

void write_vc_flow_points(mt_vector<mt_real> & xyz, std::string file)
{
  FILE* pts_file = fopen(file.c_str(), MT_FOPEN_WRITE);
  if(pts_file == NULL) treat_file_open_error(file, errno);

  int numpts = (int)(xyz.size() / 3);
  fwrite(&numpts, sizeof(int), 1, pts_file);

  mt_real* wp = xyz.data();
  double pts[3];
  PROGRESS<int> progress(numpts, "Writing points in VCFlow format: ");

  for(int i=0; i<numpts; i++)
  {
    progress.next();

    pts[0] = wp[0]; pts[1] = wp[1]; pts[2] = wp[2];
    fwrite(pts, sizeof(double), 3, pts_file);
    wp += 3;
  }
  progress.finish();

  fclose(pts_file);
}

void read_vc_flow_elements(mt_meshdata & mesh, std::string file)
{
  FILE* ele_file = fopen(file.c_str(), MT_FOPEN_READ);
  if(ele_file == NULL) treat_file_open_error(file, errno);

  int max_elem_size = 4;
  int numele;
  size_t numentries = 0;
  int n[max_elem_size]; memset(n, 0, max_elem_size*sizeof(int));

  fread(&numele, sizeof(int), 1, ele_file);

  mesh.e2n_cnt.resize(numele);
  mesh.etags.resize(numele);
  mesh.etype.resize(numele);
  mesh.e2n_con.resize(numele*max_elem_size, -1);
  mt_int* elem = mesh.e2n_con.data();
  PROGRESS<int> progress(numele, "Reading elements in VCFlow format: ");

  for(int i=0; i<numele; i++)
  {
    progress.next();
    mesh.etype[i] = Tetra;
    mt_int nodes = 4;
    mesh.e2n_cnt[i] = nodes;

    fread(n, sizeof(int), nodes, ele_file);
    // copy the element connectivity
    for(int j=0; j<nodes; j++) elem[j] = n[j];
    elem += nodes;
    numentries += nodes;
    mesh.etags[i] = 999;
  }
  progress.finish();

  mesh.e2n_con.resize(numentries);
  mesh.e2n_con.reallocate();

  fclose(ele_file);
}

void write_vc_flow_elements(mt_meshdata & mesh, const std::string file)
{
  FILE* ele_file = fopen(file.c_str(), MT_FOPEN_WRITE);
  if(ele_file == NULL) treat_file_open_error(file, errno);

  int numele = (int)(mesh.e2n_cnt.size());
  int n[4];

  fwrite(&numele, sizeof(int), 1, ele_file);
  PROGRESS<int> progress(numele, "Writing elements in VCFlow format: ");

  mt_int* elem = mesh.e2n_con.data();
  for(int i=0; i<numele; i++)
  {
    progress.next();
    mt_int nodes = mesh.e2n_cnt[i];
    if(mesh.etype[i] != Tetra)
    {
      fprintf(stderr, "Error in write_vcflow_elements: Unsupported element type!\n");
      fclose(ele_file);
      exit(1);
    }
    // element connectivity
    for(mt_int j=0; j<nodes; j++) n[j] = elem[j];
    fwrite(n, sizeof(int), nodes, ele_file);
    elem += nodes;
  }
  progress.finish();
  fclose(ele_file);
}

void write_vc_flow_adjacency(mt_meshdata & mesh, const std::string file)
{
  FILE* ele_file = fopen(file.c_str(), MT_FOPEN_WRITE);
  if(ele_file == NULL) treat_file_open_error(file, errno);

  int numele = (int)(mesh.e2n_cnt.size());
  fwrite(&numele, sizeof(int), 1, ele_file);
  PROGRESS<int> progress(numele, "Writing adjacency file in VCFlow format: ");

  //compute faces
  mt_mapping<mt_int> ele2face;
  MT_MAP<triple<mt_int>, mt_int> face_map;
  compute_faces(mesh, ele2face, face_map);
  // we have set up the element -> face directed graph, by transposing we get
  // face -> element graph.
  ele2face.transpose();
  ele2face.setup_dsp();

  const int m1 = -1;

  //loop over element indices
  for(int i = 0; i < numele; i++)
  {
    progress.next();
    if(mesh.etype[i] != Tetra)
    {
      fprintf(stderr, "Error in write_vcflow_adjacency: Unsupported element type!\n");
      fclose(ele_file);
      exit(1);
    }
    const mt_int numfaces_of_tet = ele2face.fwd_cnt[i]; // how many faces are attached to element (for a tet this should be four)
    const mt_int * faceidx       = ele2face.fwd_con.data() + ele2face.fwd_dsp[i];
    //for each face we loop over the backward vector
    for(mt_int fidx = 0; fidx < numfaces_of_tet; fidx++)
    {
      mt_int start = ele2face.bwd_dsp[faceidx[fidx]], stop = start + ele2face.bwd_cnt[faceidx[fidx]]; // gives where we need to search
      mt_int my_cnt = 0;
      for(mt_int j=start; j < stop; j++)
      {
        mt_int myidx = ele2face.bwd_con[j];
        if(myidx != (mt_int)(i))
        {
          my_cnt++;
          fwrite(&myidx, sizeof(int), 1, ele_file);
        }
      }
      //write out -1 for the rest of the faces (since they are boundary faces)
      for(mt_int j=my_cnt; j < numfaces_of_tet; j++)
        fwrite(&m1, sizeof(int),1,ele_file);
    }
  }
  progress.finish();
  fclose(ele_file);
}
