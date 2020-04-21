/**
* @file vtk_utils.cpp
* @brief VTK format specific IO utils.
* @author Aurel Neic
* @version 
* @date 2017-08-04
*/

#include "mt_utils_base.h"
#include <strings.h>  // needed for strcasecmp()

#include "io_utils.h"
#include "mesh_utils.h"
#include "vtk_utils.h"

char* next_line(char* buff, const int buffsize, FILE* fin)
{
  char* ptr = fgets(buff, buffsize, fin);
  while(ptr && strcmp(buff, "\n") == 0)
    ptr = fgets(buff, buffsize, fin);

  return ptr;
}

void vtk_process_points(mt_meshdata & mesh,
                        vtk_parameters & par,
                        char* buff, const int bufsize, FILE* fin)
{
  sscanf(buff, "POINTS %lu float", &par.numpoints);
  mesh.xyz.resize(par.numpoints * 3);

  float pts[3];
  mt_real* p = mesh.xyz.data();

  if(par.binary_data) {
    for (size_t i=0; i<par.numpoints; i++ ) {
      fread(pts, sizeof(float), 3, fin);
      p[0] = betoh(pts[0]);
      p[1] = betoh(pts[1]);
      p[2] = betoh(pts[2]);
      p += 3;
    }
  }
  else {
    for (unsigned long int i=0; i<par.numpoints; i++ )
    {
      fscanf(fin, "%f %f %f", pts, pts+1, pts+2);
      p[0] = pts[0], p[1] = pts[1], p[2] = pts[2];
      p += 3;
    }
  }
}

void vtk_process_cells(mt_meshdata & mesh,
                       vtk_parameters & par,
                       char* buff, const int bufsize, FILE* fin)
{
  unsigned long int numentries;
  sscanf(buff, "CELLS %lu %lu\n", &par.numelems, &numentries);

  mesh.e2n_cnt.resize(par.numelems);
  mesh.e2n_con.resize(numentries - par.numelems);
  mesh.etags.resize(par.numelems);
  int n[9];

  if(par.binary_data) {
    mt_int* elem = mesh.e2n_con.data();
    for(unsigned long int i=0; i<par.numelems; i++)
    {
      fread(n, sizeof(int), 1, fin);
      int nodes = betoh(n[0]);
      mesh.e2n_cnt[i] = nodes;
      fread(n, sizeof(int), nodes, fin);
      for(int j=0; j<nodes; j++) {
        elem[j] = betoh(n[j]);
      }
      elem += nodes;
    }
  }
  else {
    mt_int* elem = mesh.e2n_con.data();
    for(unsigned long int i=0; i<par.numelems; i++)
    {
      char* ptr = fgets( buff, bufsize, fin);
      if(ptr != NULL) sscanf(buff, "%d %d %d %d %d %d %d %d %d", n, n+1, n+2, n+3, n+4, n+5, n+6, n+7, n+8);
      else break;

      int nodes = n[0];
      mesh.e2n_cnt[i] = nodes;
      for(int j=0; j<nodes; j++) elem[j] = n[j+1];
      elem += nodes;
    }
  }
}

void vtk_process_polygons(mt_meshdata & mesh,
                          vtk_parameters & par,
                          char* buff, const int bufsize, FILE* fin)
{
  unsigned long int numentries;
  sscanf(buff, "POLYGONS %lu %lu\n", &par.numelems, &numentries);

  mesh.e2n_cnt.resize(par.numelems);
  mesh.e2n_con.resize(numentries - par.numelems);
  mesh.etags.resize(par.numelems);
  int n[9];

  if(par.binary_data) {
    mt_int* elem = mesh.e2n_con.data();
    for(unsigned long int i=0; i<par.numelems; i++)
    {
      fread(n, sizeof(int), 1, fin);
      int nodes = betoh(n[0]);
      mesh.e2n_cnt[i] = nodes;
      fread(n, sizeof(int), nodes, fin);
      for(int j=0; j<nodes; j++) {
        elem[j] = betoh(n[j]);
      }
      elem += nodes;
    }
  }
  else {
    mt_int* elem = mesh.e2n_con.data();
    for(unsigned long int i=0; i<par.numelems; i++)
    {
      char* ptr = fgets( buff, bufsize, fin);
      if(ptr != NULL) sscanf(buff, "%d %d %d %d %d %d %d %d %d", n, n+1, n+2, n+3, n+4, n+5, n+6, n+7, n+8);
      else break;

      int nodes = n[0];
      mesh.e2n_cnt[i] = nodes;
      for(int j=0; j<nodes; j++) elem[j] = n[j+1];
      elem += nodes;
    }
  }
}

void vtk_process_cellTypes(mt_meshdata & mesh,
                           vtk_parameters & par,
                           char* buff, const int bufsize, FILE* fin)
{
  sscanf(buff, "CELL_TYPES %lu\n", &par.numelems);
  mesh.etype.resize(par.numelems);

  if(par.binary_data) {
    mt_vector<int> etype(par.numelems);
    fread(etype.data(), sizeof(int), par.numelems, fin);

    for(unsigned long int i = 0; i<par.numelems; i++) {
      switch(betoh(int(etype[i]))) {
        case 3: // Lines are encoded as index 3
          mesh.etype[i] = Line;
          break;
        case 5: // Triangles are encoded as index 5
          mesh.etype[i] = Tri;
          break;
        case 9: // Quads are encoded as index 9
          mesh.etype[i] = Quad;
          break;
        case 10: // Tetras are encoded as index 10
          mesh.etype[i] = Tetra;
          break;
        case 13: // Prisms are encoded as index 13
          mesh.etype[i] = Prism;
          break;
        case 14: // Pyramids are encoded as index 14
          mesh.etype[i] = Pyramid;
          break;
        case 12: // Hexahedras are encoded as index 12
          mesh.etype[i] = Hexa;
          break;
      }
    }
  }
  else {
    int etype;
    for(unsigned long int i = 0; i<par.numelems; i++) {
      char* ptr = fgets(buff, bufsize, fin);
      if(ptr != NULL) sscanf(buff, "%d", &etype);
      else break;

      switch(etype) {
        case 3: // Lines are encoded as index 3
          mesh.etype[i] = Line;
          break;
        case 5: // Triangles are encoded as index 5
          mesh.etype[i] = Tri;
          break;
        case 9: // Quads are encoded as index 9
          mesh.etype[i] = Quad;
          break;
        case 10: // Tetras are encoded as index 10
          mesh.etype[i] = Tetra;
          break;
        case 13: // Prisms are encoded as index 13
          mesh.etype[i] = Prism;
          break;
        case 14: // Pyramids are encoded as index 14
          mesh.etype[i] = Pyramid;
          break;
        case 12: // Hexahedras are encoded as index 12
          mesh.etype[i] = Hexa;
          break;
      }
    }
  }
}

void vtk_process_cellData(mt_meshdata & mesh,
                          vtk_parameters & par,
                          char* buff, const int bufsize, FILE* fin)
{
  char keyword[32], rest[2048];
  mt_vector<mt_real> fiber, sheet;

  sscanf(buff, "CELL_DATA %lu", &par.numelemdata);

  while(true) {
    char* ptr = fgets(buff, bufsize, fin);
    int numread;
    if(ptr != NULL) numread = sscanf(buff, "%s %s", keyword, rest);
    else break;

    if(numread != 2) continue;

    if(strcasecmp(keyword, "SCALARS") == 0)
    {
      sscanf(buff, "SCALARS %s %s", keyword, rest);
      if(strcmp(keyword, "elemTag") == 0)
      {
        fgets(buff, bufsize, fin); // jump over lookup table line
        mesh.etags.resize(par.numelemdata);
        int ibuff;
        if(par.binary_data) {
          for(unsigned long int i=0; i<par.numelemdata; i++) {
            fread(&ibuff, sizeof(int), 1, fin);
            mesh.etags[i] = betoh(ibuff);
          }
        }
        else {
          for(unsigned long int i=0; i<par.numelemdata; i++) {
            fscanf(fin, "%d", &ibuff);
            mesh.etags[i] = ibuff;
          }
        }
      }
    }
    else if(strcasecmp(keyword, "VECTORS") == 0) {
      sscanf(buff, "VECTORS %s %s", keyword, rest);
      if(strcmp(keyword, "fiber") == 0)
      {
        fiber.resize(par.numelemdata*3);
        mt_real* fib = fiber.data();
        float fbuff[3];

        for(unsigned long int i=0; i<par.numelemdata; i++) {
          if(par.binary_data) {
            fread(fbuff, sizeof(float), 3, fin);
            fbuff[0] = betoh(fbuff[0]);
            fbuff[1] = betoh(fbuff[1]);
            fbuff[2] = betoh(fbuff[2]);
          }
          else
            fscanf(fin, "%f %f %f", fbuff, fbuff+1, fbuff+2);

          fib[0] = fbuff[0], fib[1] = fbuff[1], fib[2] = fbuff[2];
          fib += 3;
        }
      }
      if(strcmp(keyword, "sheet") == 0) {
        sheet.resize(par.numelemdata*3);
        mt_real* she = sheet.data();
        float fbuff[3];

        for(unsigned long int i=0; i<par.numelemdata; i++) {
          if(par.binary_data) {
            fread(fbuff, sizeof(float), 3, fin);
            fbuff[0] = betoh(fbuff[0]);
            fbuff[1] = betoh(fbuff[1]);
            fbuff[2] = betoh(fbuff[2]);
          }
          else
            fscanf(fin, "%f %f %f", fbuff, fbuff+1, fbuff+2);

          she[0] = fbuff[0], she[1] = fbuff[1], she[2] = fbuff[2];
          she += 3;
        }
      }
    }
    else break;
  }

  if(fiber.size()) {
    // at least fiber must have been written
    bool twofibers = fiber.size() == sheet.size();
    if(twofibers) {
      mesh.lon.resize(par.numelemdata * 6);
      mt_real *lon = mesh.lon.data(), *fib = fiber.data(), *she = sheet.data();

      for(unsigned long int i=0; i<par.numelemdata; i++) {
        lon[0] = fib[0];
        lon[1] = fib[1];
        lon[2] = fib[2];
        lon[3] = she[0];
        lon[4] = she[1];
        lon[5] = she[2];
        lon += 6;
        fib += 3;
        she += 3;
      }
    }
    else
      mesh.lon.assign(fiber.begin(), fiber.end());
  }
  else {
    std::cerr << "Warning: Could not read mesh fibers!!" << std::endl;
    std::cerr << "Generating default fibers." << std::endl;
    generate_default_fib(par.numelemdata, mesh.lon);
  }
}

void vtk_process(mt_meshdata & mesh,
                 vtk_parameters & par,
                 char* buff, const int bufsize, FILE* fin,
                 progress<int> & progress)
{
  char keyword[256];
  char rest[2048];

  int n = sscanf(buff, "%s %s", keyword, rest);
  if(n != 2) return;

  // process the keyword
  if(strcmp(keyword, "POINTS") == 0) {
    progress.next();
    par.points_read = true;
    vtk_process_points(mesh, par, buff, bufsize, fin);
  }
  else if(strcmp(keyword, "CELLS") == 0) {
    progress.next();
    par.elems_read = true;
    vtk_process_cells(mesh, par, buff, bufsize, fin);
  }
  else if(strcmp(keyword, "CELL_TYPES") == 0) {
    progress.next();
    par.elem_types_read = true;
    vtk_process_cellTypes(mesh, par, buff, bufsize, fin);
  }
  else if(strcmp(keyword, "CELL_DATA") == 0) {
    progress.next();
    vtk_process_cellData(mesh, par, buff, bufsize, fin);
  }
  else if(strcmp(keyword, "POLYGONS") == 0) {
    progress.next();
    par.elems_read = true;
    par.is_polydata = true;
    vtk_process_polygons(mesh, par, buff, bufsize, fin);
  }
}

void readVTKmesh(mt_meshdata & mesh, std::string file)
{
  FILE* vtk_file = fopen(file.c_str(), MT_FOPEN_READ);
  if(vtk_file == NULL) treat_file_open_error(file, errno);

  const int bufsize = 2048;
  char buffer[bufsize], *ptr;
  vtk_parameters par;

  // skip unused lines
  ptr = next_line(buffer, bufsize, vtk_file);
  ptr = next_line(buffer, bufsize, vtk_file);
  // here we can treat binary / ascii
  ptr = next_line(buffer, bufsize, vtk_file);
  if(ptr != NULL) {
    char format[32];
    sscanf(buffer, "%s", format);
    if(strcasecmp(format, "binary") == 0)
      par.binary_data = true;
  }
  ptr = next_line(buffer, bufsize, vtk_file);

  // start with processing
  PROGRESS<int> progress(5, "Reading vtk file: ");

  ptr = next_line(buffer, bufsize, vtk_file);
  while(ptr != NULL && strlen(ptr) ) {
    vtk_process(mesh, par, buffer, bufsize, vtk_file, progress);
    ptr = next_line(buffer, bufsize, vtk_file);
  }
  progress.finish();
  fclose(vtk_file);

  if(par.points_read && par.elems_read) {
    if(par.is_polydata) {
      mesh.etype.resize(mesh.e2n_cnt.size());
      for(size_t i=0; i<mesh.etype.size(); i++) {
        switch(mesh.e2n_cnt[i]) {
          case 3: mesh.etype[i] = Tri; break;
          case 4: mesh.etype[i] = Quad; break;

          default:
            fprintf(stderr, "%s error: meshtool support only Tri and Quad Polygons! Aborting!\n",
                    __func__);
            exit(1);
        }
      }
    }
    mt_vector<mt_int> e2n_con(mesh.e2n_con);
    perm_con_vtk_to_carp(mesh, e2n_con);
  }
  else {
    fprintf(stderr, "%s error: File %s did not contain a complete mesh (points and elements)!\n"
                    "Aborting!\n", __func__, file.c_str());
    exit(1);
  }
}

void writeVTKmesh(const mt_meshdata & mesh, std::string file)
{
  FILE* vtk_file = fopen(file.c_str(), MT_FOPEN_WRITE);
  if(vtk_file == NULL) treat_file_open_error(file, errno);

  mt_vector<mt_int> e2n_con;
  perm_con_carp_to_vtk(mesh, e2n_con);

  unsigned long int numnodes = mesh.xyz.size() / 3;
  unsigned long int numelems = mesh.e2n_cnt.size();
  bool twofibers = mesh.lon.size() == (numelems*6);

  size_t num_next = twofibers ? (numnodes + 6*numelems) : (numnodes + 5*numelems);
  PROGRESS<size_t> progress(num_next, "Writing vtk file in text format: ");

  fprintf (vtk_file, "# vtk DataFile Version 3.0\n");
  fprintf (vtk_file, "vtk output\n");
  fprintf (vtk_file, "ASCII\n");
  fprintf (vtk_file, "DATASET UNSTRUCTURED_GRID\n\n");
  fprintf (vtk_file, "POINTS %lu float\n", numnodes);

  float pts[3];
  const mt_real* p = mesh.xyz.data();

  for (unsigned long int i=0; i<numnodes; i++ ) {
    progress.next();

    pts[0] = p[0], pts[1] = p[1], pts[2] = p[2];
    fprintf(vtk_file, "%.3f %.3f %.3f\n", pts[0], pts[1], pts[2]);
    p += 3;
  }

  fprintf (vtk_file, "CELL_TYPES %lu\n", numelems);
  unsigned long int valcount = 0;
  for(unsigned long int i=0; i<numelems; i++) {
    progress.next();

    mt_int actnodes = mesh.e2n_cnt[i];
    switch(mesh.etype[i]) {
      // Line for Purkinje
      case Line:
        fprintf(vtk_file, "3\n"); // Lines are encoded as index 3
        break;
      case Tri:
        fprintf(vtk_file, "5\n"); // Triangles are encoded as index 5
        break;
      case Quad:
        fprintf(vtk_file, "9\n"); // Quads are encoded as index 9
        break;
      case Tetra:
        fprintf(vtk_file, "10\n"); // Tetras are encoded as index 10
        break;
      case Pyramid:
        fprintf(vtk_file, "14\n"); // Pyramids are encoded as index 14
        break;
      case Prism:
        fprintf(vtk_file, "13\n"); // Prisms are encoded as index 13
        break;
      case Hexa:
        fprintf(vtk_file, "12\n"); // Hexahedras are encoded as index 12
        break;
    }
    valcount += actnodes+1;
  }

  fprintf(vtk_file, "CELLS %lu %lu\n", numelems, valcount);
  const mt_int* elem = e2n_con.data();
  for(unsigned long int i=0; i<numelems; i++)
  {
    progress.next();

    int nodes = mesh.e2n_cnt[i];
    fprintf(vtk_file, "%d ", nodes);

    for(int j=0; j<nodes; j++)
      fprintf(vtk_file, "%ld ", (long int)elem[j]);

    fprintf(vtk_file, "\n");
    elem += nodes;
  }

  fprintf (vtk_file, "CELL_DATA %lu \n", numelems);
  fprintf (vtk_file, "SCALARS elemTag int 1\n");
  fprintf (vtk_file, "LOOKUP_TABLE default\n");
  for (unsigned long int i=0; i<numelems; i++ )
  {
    progress.next();
    fprintf(vtk_file, "%d \n", (int)mesh.etags[i]);
  }

  // write fiber data
  if(twofibers) {
    fprintf (vtk_file, "VECTORS fiber float\n");
    for (unsigned long int i=0; i<numelems; i++ ) {
      progress.next();
      fprintf(vtk_file, "%f %f %f\n", mesh.lon[i*6], mesh.lon[i*6+1], mesh.lon[i*6+2]);
    }
    fprintf (vtk_file, "VECTORS sheet float\n");
    for (unsigned long int i=0; i<numelems; i++ ) {
      progress.next();
      fprintf(vtk_file, "%f %f %f\n", mesh.lon[i*6+3], mesh.lon[i*6+4], mesh.lon[i*6+5]);
    }
  }
  else {
    fprintf (vtk_file, "VECTORS fiber float\n");
    for (unsigned long int i=0; i<numelems; i++ ) {
      progress.next();
      fprintf(vtk_file, "%f %f %f\n", mesh.lon[i*3], mesh.lon[i*3+1], mesh.lon[i*3+2]);
    }
  }

  progress.finish();
  fclose(vtk_file);
}

void writeVTKmesh_binary(const mt_meshdata & mesh, std::string file)
{
  FILE* vtk_file = fopen(file.c_str(), MT_FOPEN_WRITE);
  if(vtk_file == NULL) treat_file_open_error(file, errno);

  mt_vector<mt_int> e2n_con;
  perm_con_carp_to_vtk(mesh, e2n_con);

  unsigned long int numnodes = mesh.xyz.size() / 3;
  unsigned long int numelems = mesh.e2n_cnt.size();
  bool twofibers = mesh.lon.size() == (numelems*6);

  size_t num_next = twofibers ? (numnodes + 6*numelems) : (numnodes + 5*numelems);
  PROGRESS<size_t> progress(num_next, "Writing vtk file in binary format: ");

  fprintf (vtk_file, "# vtk DataFile Version 3.0\n");
  fprintf (vtk_file, "vtk output\n");
  fprintf (vtk_file, "binary\n");
  fprintf (vtk_file, "DATASET UNSTRUCTURED_GRID\n\n");
  fprintf (vtk_file, "POINTS %lu float\n", numnodes);

  mt_vector<float> pts(mesh.xyz.size());
  const mt_real* p = mesh.xyz.data();

  for (unsigned long int i=0; i<numnodes; i++ ) {
    progress.next();

    pts[i*3+0] = htobe(float(p[0]));
    pts[i*3+1] = htobe(float(p[1]));
    pts[i*3+2] = htobe(float(p[2]));
    p += 3;
  }
  fwrite(pts.data(), sizeof(float), pts.size(), vtk_file);
  fprintf(vtk_file, "\n");

  mt_vector<int> write_buff(numelems);

  fprintf (vtk_file, "CELL_TYPES %lu\n", numelems);
  unsigned long int valcount = 0;
  int vtk_type = 0;
  for(unsigned long int i=0; i<numelems; i++) {
    progress.next();

    elem_t etype = mesh.etype[i];
    switch(etype) {
      // Line for Purkinje
      case Line:
        vtk_type =  3; // Lines are encoded as index 3
        break;
      case Tri:
        vtk_type = 5; // Triangles are encoded as index 5
        break;
      case Tetra:
        vtk_type = 10; // Tetras are encoded as index 10
        break;
      case Quad:
        vtk_type = 9; // Quads are encoded as index 9
        break;
      case Pyramid:
        vtk_type = 14; // Pyramids are encoded as index 14
        break;
      case Prism:
        vtk_type = 13; // Prisms are encoded as index 13
        break;
      case Hexa:
        vtk_type = 12; // Hexahedras are encoded as index 12
        break;
    }
    write_buff[i] = htobe(vtk_type);
    valcount += mesh.e2n_cnt[i]+1;
  }

  fwrite(write_buff.data(), sizeof(int), numelems, vtk_file);
  fprintf(vtk_file, "\n");

  fprintf(vtk_file, "CELLS %lu %lu\n", numelems, valcount);
  write_buff.resize(numelems + mesh.e2n_con.size());
  size_t widx = 0;

  const mt_int* elem = e2n_con.data();
  for(unsigned long int i=0; i<numelems; i++)
  {
    progress.next();
    int nodes = mesh.e2n_cnt[i];
    write_buff[widx++] = htobe(nodes);

    for(int j=0; j<mesh.e2n_cnt[i]; j++)
      write_buff[widx++] = htobe(int(elem[j]));

    elem += mesh.e2n_cnt[i];
  }
  fwrite(write_buff.data(), sizeof(int), write_buff.size(), vtk_file);
  fprintf(vtk_file, "\n");

  fprintf (vtk_file, "CELL_DATA %lu \n", numelems);
  fprintf (vtk_file, "SCALARS elemTag int 1\n");
  fprintf (vtk_file, "LOOKUP_TABLE default\n");

  write_buff.resize(numelems);
  for (unsigned long int i=0; i<numelems; i++ ) {
    progress.next();
    write_buff[i] = htobe(int(mesh.etags[i]));
  }
  fwrite(write_buff.data(), sizeof(int), numelems, vtk_file);

  pts.resize(numelems*3);
  // write fiber data
  if(twofibers) {
    fprintf (vtk_file, "VECTORS fiber float\n");
    for (unsigned long int i=0; i<numelems; i++ ) {
      progress.next();
      pts[i*3+0] = htobe(float(mesh.lon[i*6+0])), pts[i*3+1] = htobe(float(mesh.lon[i*6+1])), pts[i*3+2] = htobe(float(mesh.lon[i*6+2]));
    }
    fwrite(pts.data(), sizeof(float), pts.size(), vtk_file);
    fprintf(vtk_file, "\n");

    fprintf (vtk_file, "VECTORS sheet float\n");
    for (unsigned long int i=0; i<numelems; i++ ) {
      progress.next();
      pts[i*3+0] = htobe(float(mesh.lon[i*6+3])), pts[i*3+1] = htobe(float(mesh.lon[i*6+4])), pts[i*3+2] = htobe(float(mesh.lon[i*6+5]));
    }
    fwrite(pts.data(), sizeof(float), pts.size(), vtk_file);
    fprintf(vtk_file, "\n");
  }
  else {
    fprintf (vtk_file, "VECTORS fiber float\n");
    for (unsigned long int i=0; i<numelems; i++ ) {
      progress.next();
      pts[i*3+0] = htobe(float(mesh.lon[i*3+0])), pts[i*3+1] = htobe(float(mesh.lon[i*3+1])), pts[i*3+2] = htobe(float(mesh.lon[i*3+2]));
    }
    fwrite(pts.data(), sizeof(float), pts.size(), vtk_file);
    fprintf(vtk_file, "\n");
  }

  progress.finish();
  fclose(vtk_file);
}

void writeVTKPolydataMesh(const mt_meshdata & mesh, std::string file)
{
  FILE* vtk_file = fopen(file.c_str(), MT_FOPEN_WRITE);
  if(vtk_file == NULL) treat_file_open_error(file, errno);

  size_t numnodes = mesh.xyz.size() / 3;
  size_t numelems = mesh.e2n_cnt.size();

  fprintf (vtk_file, "# vtk DataFile Version 2.3\n");
  fprintf (vtk_file, "vtk output\n");
  fprintf (vtk_file, "ASCII\n");
  fprintf (vtk_file, "DATASET POLYDATA\n\n");
  fprintf (vtk_file, "POINTS %lu float\n", numnodes);

  float pts[3];
  const mt_real* p = mesh.xyz.data();

  for (size_t i=0; i<numnodes; i++ ) {
    pts[0] = p[0], pts[1] = p[1], pts[2] = p[2];
    fprintf(vtk_file, "%f %f %f\n", pts[0], pts[1], pts[2]);
    p += 3;
  }

  /*fprintf (vtk_file, "CELL_TYPES %d\n", numelems);*/
  size_t valcount = 0;
  for (size_t i=0; i<numelems; i++) {
    mt_int actnodes = mesh.e2n_cnt[i];
    valcount += actnodes+1;
  }

  fprintf(vtk_file, "POLYGONS %lu %lu\n", numelems, valcount);
  const mt_int* elem = mesh.e2n_con.data();
  for(size_t i=0; i<numelems; i++)
  {
    mt_int nodes = mesh.e2n_cnt[i];
    fprintf(vtk_file, "%d ", int(nodes));

    for(mt_int j=0; j<nodes; j++)
      fprintf(vtk_file, "%ld ", (long int)elem[j]);

    fprintf(vtk_file, "\n");
    elem += nodes;
  }


  fprintf (vtk_file, "CELL_DATA %lu \n", numelems);
  fprintf (vtk_file, "SCALARS %s int %d\n", "elemTag", 1);
  fprintf (vtk_file, "LOOKUP_TABLE default\n");
  for (size_t i=0; i<numelems; i++ ) {
    fprintf(vtk_file, "%d \n", (int)mesh.etags[i]);
  }

  fclose(vtk_file);
}

void write_VTU_points(const mt_vector<mt_real> & xyz, bool binary, FILE* output_file)
{
  size_t npoints = xyz.size() / 3;

  if(binary == false)
  {
    fprintf(output_file, "      <Points>\n");
    fprintf(output_file, "        <DataArray type=\"Float32\" Name=\"Points\"");
    fprintf(output_file, " NumberOfComponents=\"3\" format=\"ascii\">");
    for(size_t i=0; i < npoints; i++)
        fprintf(output_file, "\n         %f %f %f", xyz[i*3+0], xyz[i*3+1], xyz[i*3+2]);
    fprintf(output_file, "\n        </DataArray>\n");
    fprintf(output_file, "      </Points>\n");
  }
  else {
    mt_vector<unsigned char> data_buf, encode_buf;
    mt_vector<float> wbuff;
    wbuff.assign(xyz.begin(), xyz.end());

    fprintf(output_file, "      <Points>\n");
    fprintf(output_file, "        <DataArray type=\"Float32\" Name=\"Points\"");
    fprintf(output_file, " NumberOfComponents=\"3\" format=\"binary\">\n");
    fprintf(output_file, "          ");

    write_vtu_data(output_file, data_buf, encode_buf, wbuff);

    fprintf(output_file, "        </DataArray>\n");
    fprintf(output_file, "      </Points>\n");
  }
}

void write_VTU_types(const mt_vector<elem_t> & types, bool binary, FILE* output_file)
{
  size_t nelems = types.size();

  if(binary == false) {
    fprintf(output_file, "        <DataArray type=\"UInt8\" Name=\"types\"");
    fprintf(output_file, " format=\"ascii\">");
    int line_count = 0;

    for(size_t j=0; j < nelems; j++) {
      if(line_count % 6 == 0)
        fprintf(output_file, "\n         ");

      line_count++;

      switch(types[j]) {
        case Line:    fprintf(output_file, " 3");  break;
        case Tri:     fprintf(output_file, " 5");  break;
        case Quad:    fprintf(output_file, " 9");  break;
        case Tetra:   fprintf(output_file, " 10"); break;
        case Hexa:    fprintf(output_file, " 12"); break;
        case Prism:   fprintf(output_file, " 13"); break;
        case Pyramid: fprintf(output_file, " 14"); break;
        default:
          fprintf( stderr, "ERROR in %s: Element type not supported\n", __func__);
          exit(EXIT_FAILURE);
      }
    }
  }
  else {
    fprintf(output_file, "        <DataArray type=\"UInt8\" Name=\"types\"");
    fprintf(output_file, " format=\"binary\">\n");
    mt_vector<unsigned char> num_array(types.size());

    // Write out element connectivity
    for(size_t j=0; j<nelems; j++) {
      switch(types[j]) {
        case Line:    num_array[j] = 3;  break;
        case Tri:     num_array[j] = 5;  break;
        case Quad:    num_array[j] = 9;  break;
        case Tetra:   num_array[j] = 10; break;
        case Hexa:    num_array[j] = 12; break;
        case Prism:   num_array[j] = 13; break;
        case Pyramid: num_array[j] = 14; break;
        default:
          fprintf( stderr, "ERROR in %s: Element type not supported\n", __func__);
          exit(EXIT_FAILURE);
        }
    }

    mt_vector<unsigned char> data_buf, encode_buf;
    write_vtu_data(output_file, data_buf, encode_buf, num_array);
  }
  fprintf(output_file, "        </DataArray>\n");
}

void write_VTU_cell_offsets(const mt_meshdata & mesh, bool binary, FILE* output_file)
{
  int node_count = 0;
  int line_count = 0;
  size_t nelems = mesh.e2n_cnt.size();

  if(binary == false) {
    fprintf(output_file, "        <DataArray type=\"Int64\" Name=\"offsets\"");
    fprintf(output_file, " format=\"ascii\">");
    for(size_t i=0; i < nelems; i++) {
      node_count += mesh.e2n_cnt[i];
      if(line_count % 6 == 0)
        fprintf(output_file, "\n         ");
      fprintf(output_file, " %d", node_count);
      line_count++;
    }
  }
  else {
    fprintf(output_file, "        <DataArray type=\"Int64\" Name=\"offsets\"");
    fprintf(output_file, " format=\"binary\">\n");

    mt_vector<long int> num_array(nelems);
    for(size_t i=0; i < nelems; i++) {
      node_count += mesh.e2n_cnt[i];
      num_array[i] = node_count;
    }

    mt_vector<unsigned char> data_buf, encode_buf;
    write_vtu_data(output_file, data_buf, encode_buf, num_array);
  }

  fprintf(output_file, "        </DataArray>\n");
}

void write_VTU_cells(const mt_meshdata & mesh, bool binary, FILE* output_file)
{
  size_t   nelems   = mesh.e2n_cnt.size();

  fprintf(output_file, "      <Cells>\n");

  write_VTU_types(mesh.etype, binary, output_file);

  write_VTU_cell_offsets(mesh, binary, output_file);

  mt_vector<mt_int> e2n_con;
  perm_con_carp_to_vtk(mesh, e2n_con);

  if(binary == false) {
    fprintf(output_file, "        <DataArray type=\"Int64\" Name=\"connectivity\"");
    fprintf(output_file, " format=\"ascii\"");
    fprintf(output_file, " RangeMin=\"0\" RangeMax=\"%ld\">", (long int)nelems);
    for(size_t k=0, ridx=0; k < nelems; k++) {
      if(!(k%6))
        fprintf(output_file, "\n         ");

      for(mt_int j=0; j < mesh.e2n_cnt[k]; j++, ridx++)
        fprintf(output_file, " %ld", (long int) e2n_con[ridx]);
    }
  }
  else {
    fprintf(output_file, "        <DataArray type=\"Int64\" Name=\"connectivity\"");
    fprintf(output_file, " format=\"binary\"");
    fprintf(output_file, " RangeMin=\"0\" RangeMax=\"%ld\">\n", (long int) nelems);
    mt_vector<long int> num_array;
    num_array.assign(e2n_con.begin(), e2n_con.end());

    mt_vector<unsigned char> data_buf, encode_buf;
    write_vtu_data(output_file, data_buf, encode_buf, num_array);
  }

  fprintf(output_file, "        </DataArray>\n");
  fprintf(output_file, "      </Cells>\n");
}

void finish_VTU(FILE* output_file)
{
  fprintf(output_file, "    </Piece>\n");
  fprintf(output_file, "  </UnstructuredGrid>\n");
  fprintf(output_file, "</VTKFile>\n");
  fclose (output_file);
}

void writeXMLVTKmesh(mt_meshdata & mesh, std::string file)
{
  FILE* vtk_file = fopen(file.c_str(), MT_FOPEN_WRITE);
  if(vtk_file == NULL) treat_file_open_error(file, errno);

  mt_vector<mt_int> e2n_con;
  perm_con_carp_to_vtk(mesh, e2n_con);

  unsigned long int numnodes = mesh.xyz.size() / 3;
  unsigned long int numelems = mesh.e2n_cnt.size();
  short numFib = mesh.lon.size() == (numelems*6) ? 2 : 1;

  #if MT_ENDIANNESS == 0
  const char* endian = "LittleEndian";
  #else
  const char* endian = "BigEndian";
  #endif

  const char *istr, *fstr;

  if(sizeof(mt_int) == sizeof(int))
    istr = "Int32";
  else if(sizeof(mt_int) == sizeof(long))
      istr = "Int64";
    else {
      std::cerr << "writeXMLVTKmesh: Unknown integer datatype. Aborting." << std::endl;
      exit(1);
    }

  if(sizeof(mt_real) == sizeof(float))
    fstr = "Float32";
  else if(sizeof(mt_real) == sizeof(double))
      fstr = "Float64";
    else {
      std::cerr << "writeXMLVTKmesh: Unknown float datatype. Aborting." << std::endl;
      exit(1);
    }

  if(mesh.e2n_dsp.size() == 0) {
    mesh.e2n_dsp.resize(mesh.e2n_cnt.size());
    bucket_sort_offset(mesh.e2n_cnt, mesh.e2n_dsp);
  }

  PROGRESS<short> prg(8, "Writing binary vtu file: ");

  mt_vector<unsigned char>      data_buf;
  mt_vector<unsigned char>      encode_buf;

  fprintf(vtk_file, "<?xml version=\"1.0\"?>\n");
  fprintf(vtk_file, "<VTKFile type=\"UnstructuredGrid\" byte_order=\"%s\" header_type=\"UInt64\">\n", endian);
  fprintf(vtk_file, "  <UnstructuredGrid>\n");
  fprintf(vtk_file, "    <Piece NumberOfPoints=\"%ld\" NumberOfCells=\"%ld\">\n", numnodes, numelems);

  fprintf(vtk_file, "      <Points>\n");
  fprintf(vtk_file, "        <DataArray type=\"%s\" format=\"binary\" NumberOfComponents=\"3\">\n", fstr);

  write_vtu_data(vtk_file, data_buf, encode_buf, mesh.xyz);
  prg.next();

  fprintf(vtk_file, "        </DataArray>\n");
  fprintf(vtk_file, "      </Points>\n");
  fprintf(vtk_file, "      <Cells>\n");
  fprintf(vtk_file, "        <DataArray type=\"%s\" format=\"binary\" Name=\"connectivity\">\n", istr);

  write_vtu_data(vtk_file, data_buf, encode_buf, e2n_con);
  prg.next();

  fprintf(vtk_file, "        </DataArray>\n");
  fprintf(vtk_file, "        <DataArray type=\"%s\" format=\"binary\" Name=\"offsets\">\n", istr);

  mt_vector<mt_int> dsp(mesh.e2n_dsp.size());
  for(size_t i=0; i<mesh.e2n_dsp.size(); i++)
    dsp[i] = mesh.e2n_cnt[i] + mesh.e2n_dsp[i];

  write_vtu_data(vtk_file, data_buf, encode_buf, dsp);
  prg.next();

  fprintf(vtk_file, "        </DataArray>\n");
  fprintf(vtk_file, "        <DataArray type=\"UInt8\" format=\"binary\" Name=\"types\">\n");
  mt_vector<unsigned char> tp(mesh.etype.size());
  for(size_t i=0; i<tp.size(); i++) {
    switch(mesh.etype[i]) {
      case Line:
        tp[i] =  3; // Lines are encoded as index 3
        break;
      case Tri:
        tp[i] = 5; // Triangles are encoded as index 5
        break;
      case Tetra:
        tp[i] = 10; // Tetras are encoded as index 10
        break;
      case Quad:
        tp[i] = 9; // Quads are encoded as index 9
        break;
      case Pyramid:
        tp[i] = 14; // Pyramids are encoded as index 14
        break;
      case Prism:
        tp[i] = 13; // Prisms are encoded as index 13
        break;
      case Hexa:
        tp[i] = 12; // Hexahedras are encoded as index 12
        break;
    }
  }

  write_vtu_data(vtk_file, data_buf, encode_buf, tp);
  prg.next();

  fprintf(vtk_file, "        </DataArray>\n");
  fprintf(vtk_file, "      </Cells>\n");
  fprintf(vtk_file, "      <CellData Scalars=\"elemTag\" Vectors=\"fiber\">\n");
  fprintf(vtk_file, "        <DataArray type=\"%s\" format=\"binary\" Name=\"elemTag\">\n", istr);

  write_vtu_data(vtk_file, data_buf, encode_buf, mesh.etags);
  prg.next();

  fprintf(vtk_file, "        </DataArray>\n");
  fprintf(vtk_file, "        <DataArray type=\"%s\" format=\"binary\" NumberOfComponents=\"3\" Name=\"fiber\">\n", fstr);
  mt_vector<mt_real> fib(numelems*3);
  for(size_t i=0; i<numelems; i++) {
    fib[i*3+0] = mesh.lon[i*3*numFib+0];
    fib[i*3+1] = mesh.lon[i*3*numFib+1];
    fib[i*3+2] = mesh.lon[i*3*numFib+2];
  }
  write_vtu_data(vtk_file, data_buf, encode_buf, fib);
  prg.next();

  fprintf(vtk_file, "        </DataArray>\n");
  if(numFib == 2) {
    fprintf(vtk_file, "      <DataArray type=\"%s\" format=\"binary\" NumberOfComponents=\"3\" Name=\"sheet\">\n", fstr);
    for(size_t i=0; i<numelems; i++) {
      fib[i*3+0] = mesh.lon[i*6+3];
      fib[i*3+1] = mesh.lon[i*6+4];
      fib[i*3+2] = mesh.lon[i*6+5];
    }
    write_vtu_data(vtk_file, data_buf, encode_buf, fib);
    fprintf(vtk_file, "      </DataArray>\n");
  }
  prg.next();

  fprintf(vtk_file, "      </CellData>\n");
  fprintf(vtk_file, "    </Piece>\n");
  fprintf(vtk_file, "  </UnstructuredGrid>\n");
  fprintf(vtk_file, "</VTKFile>\n");
  fclose(vtk_file);

  prg.finish();
}

void write_VTU_elemtags (const mt_vector<mt_int> & tags, const char *name, bool binary,
                        FILE* output_file)
{
  std::vector<int> wbuff;
  wbuff.assign(tags.begin(), tags.end());

  bool is_float = false;
  write_VTU_dataset(wbuff.data(), wbuff.size(), 1, name, is_float, binary, output_file);
}

void write_VTU_londata(const mt_vector<mt_real> & lon, const size_t nelem,
                       const int num_axes, const char* prefix, bool binary,
                       FILE* output_file)
{
  if (num_axes != 1 && num_axes != 2) {
    fprintf(stderr, "Currently only 2 fiber directions supported.\n");
    exit(EXIT_FAILURE);
  }

  std::string name = prefix;
  mt_vector<float> wbuff(3 * nelem);
  const int  fibvals  = num_axes * 3;
  const bool is_float = true;

  name += "fibres";
  for(size_t k =0; k < nelem; k++) {
    wbuff[k*3+0] = lon[k*fibvals+0];
    wbuff[k*3+1] = lon[k*fibvals+1];
    wbuff[k*3+2] = lon[k*fibvals+2];
  }

  write_VTU_dataset(wbuff.data(), nelem, 3, name.c_str(), is_float, binary,
                    output_file);

  if(num_axes == 2) {
    name = prefix;
    name += "sheets";
    for(size_t k =0; k < nelem; k++) {
      wbuff[k*3+0] = lon[k*fibvals+3];
      wbuff[k*3+1] = lon[k*fibvals+4];
      wbuff[k*3+2] = lon[k*fibvals+5];
    }

    write_VTU_dataset(wbuff.data(), nelem, 3, name.c_str(), is_float, binary,
                      output_file);
  }
}

void write_VTU_mesh(const mt_meshdata & mesh, bool binary, bool finish_file, FILE* vtk_file)
{
  long int npts   = mesh.xyz.size() / 3;
  long int nelems = mesh.e2n_cnt.size();

  // write header
#if MT_ENDIANNESS == 0
  const char* endian = "LittleEndian";
#else
  const char* endian = "BigEndian";
#endif

  fprintf(vtk_file, "<VTKFile type=\"UnstructuredGrid\"");
  fprintf(vtk_file, " version=\"1.0\" byte_order=\"%s\" \n", endian);
  fprintf(vtk_file, " header_type=\"UInt64\">\n");
  fprintf(vtk_file, "  <UnstructuredGrid>\n");
  fprintf(vtk_file, "    <Piece NumberOfPoints=\"%ld\" ", npts);
  fprintf(vtk_file, " NumberOfCells=\"%ld\">\n", nelems);
  // Write content
  write_VTU_points(mesh.xyz, binary, vtk_file);
  write_VTU_cells (mesh, binary, vtk_file);

  if(finish_file) {
    int num_axes = mesh.lon.size() == mesh.e2n_cnt.size() * 6 ? 2 : 1;

    fprintf(vtk_file, "      <CellData>\n");
    write_VTU_elemtags (mesh.etags, "elemTags", binary, vtk_file);
    write_VTU_londata(mesh.lon, mesh.e2n_cnt.size(), num_axes, "", binary, vtk_file);
    fprintf(vtk_file, "      </CellData>\n");

    finish_VTU(vtk_file);
  }
}

void write_VTU_mesh(const mt_meshdata & mesh, std::string file, bool binary)
{
  FILE* vtk_file = fopen(file.c_str(), MT_FOPEN_WRITE);
  if(vtk_file == NULL) treat_file_open_error(file, errno);

  write_VTU_mesh(mesh, binary, true, vtk_file);
}

unsigned long int read_DataArray_header(FILE* vtk_file,
                                        const std::string & header_type)
{
  short header_len = 0;

  if(header_type.compare("UInt64") == 0)
    header_len = 12;
  else if(header_type.compare("UInt32") == 0)
    header_len = 8;
  else
    treat_file_error("vtu read error: unsupported DataArray header type.", vtk_file);

  mt_vector<unsigned char>     data_buf(header_len);
  mt_vector<unsigned long int> data_len;

  seek_over_char(vtk_file, ' ');
  long arr_start = ftell(vtk_file);
  fread(data_buf.data(), header_len, 1, vtk_file);
  fseek(vtk_file, arr_start, SEEK_SET);
  decode_base64(data_buf, data_len);

  return data_len[0];
}

template<class T>
void read_DataArray_data(FILE* vtk_file,
                         const std::string & header_type,
                         const std::string & data_type,
                         const unsigned long int data_len,
                         mt_vector<unsigned char> & data_buf,
                         mt_vector<unsigned char> & decode_buf,
                         mt_vector<T> & data)
{
  short header_len = 0;

  if(header_type.compare("UInt64") == 0)
    header_len = sizeof(unsigned long int);
  else if(header_type.compare("UInt32") == 0)
    header_len = sizeof(unsigned int);
  else
    treat_file_error("vtu read error: unsupported DataArray header type.", vtk_file);

  data_buf.resize((data_len + header_len + 2) / 3 * 4 + 2);

  fgets((char*)data_buf.data(), data_buf.size(), vtk_file);
  decode_base64(data_buf, decode_buf);

  if(data_type.compare("Float32") == 0) {
    float* start = (float*)(decode_buf.begin() + header_len);
    float* end   = (float*)decode_buf.end();
    data.assign(start, end);
  }
  else if(data_type.compare("Float64") == 0) {
    double* start = (double*)(decode_buf.begin() + header_len);
    double* end   = (double*)decode_buf.end();
    data.assign(start, end);
  }
  else if(data_type.compare("UInt8") == 0) {
    unsigned char* start = (unsigned char*)(decode_buf.begin() + header_len);
    unsigned char* end   = (unsigned char*)decode_buf.end();
    data.assign(start, end);
  }
  else if(data_type.compare("Int32") == 0) {
    int* start = (int*)(decode_buf.begin() + header_len);
    int* end   = (int*)decode_buf.end();
    data.assign(start, end);
  }
  else if(data_type.compare("Int64") == 0) {
    long int* start = (long int*)(decode_buf.begin() + header_len);
    long int* end   = (long int*)decode_buf.end();
    data.assign(start, end);
  }
  else
    treat_file_error("vtu read error: unsupported DataArray data type.", vtk_file);

}

unsigned short DataArray_type_size(const std::string & data_type)
{
  if(data_type.compare("Float32") == 0) {
    return sizeof(float);
  }
  else if(data_type.compare("Float64") == 0) {
    return sizeof(double);
  }
  else if(data_type.compare("UInt8") == 0) {
    return sizeof(unsigned char);
  }
  else if(data_type.compare("Int32") == 0) {
    return sizeof(int);
  }
  else if(data_type.compare("Int64") == 0) {
    return sizeof(long int);
  }
  return 0;
}


bool process_xmlvtk_piece(mt_meshdata & mesh,
                          FILE* vtk_file,
                          const std::string & header_type,
                          PROGRESS<short> & prg)
{
  const int bufsize = 2048;
  char buffer[bufsize], *ptr = NULL;
  MT_MAP<std::string, std::string> fields;
  bool finished = false;

  mt_vector<unsigned char> data_buf;
  mt_vector<unsigned char> decode_buf;
  unsigned long int data_len;

  ptr = next_line(buffer, bufsize, vtk_file);
  while(ptr && !finished) {
    parse_xml_line(ptr, fields);

    if(fields["tag"].compare("Points") == 0) {
      ptr = next_line(buffer, bufsize, vtk_file);
      parse_xml_line(ptr, fields);

      data_len = read_DataArray_header(vtk_file, header_type);
      if(data_len / DataArray_type_size(fields["type"]) != mesh.xyz.size())
        treat_file_error("vtu read error: size of Points DataArray"
                         "does not match number of points.", vtk_file);

      read_DataArray_data(vtk_file, header_type, fields["type"],
                          data_len, data_buf, decode_buf, mesh.xyz);
      prg.next();
    }
    else if(fields["tag"].compare("Cells") == 0)
    {
      while(ptr) {
        ptr = next_line(buffer, bufsize, vtk_file);
        parse_xml_line(ptr, fields);

        if((fields["tag"].compare("DataArray") == 0) &&
           (fields["Name"].compare("connectivity") == 0))
        {
          data_len = read_DataArray_header(vtk_file, header_type);
          read_DataArray_data(vtk_file, header_type, fields["type"],
                              data_len, data_buf, decode_buf, mesh.e2n_con);
          prg.next();
        }
        else if((fields["tag"].compare("DataArray") == 0) &&
                (fields["Name"].compare("offsets") == 0))
        {
          mt_vector<mt_int> offset;

          data_len = read_DataArray_header(vtk_file, header_type);
          if(data_len / DataArray_type_size(fields["type"]) != mesh.e2n_cnt.size())
            treat_file_error("vtu read error: size of offsets DataArray"
                             "does not match number of elements.", vtk_file);

          read_DataArray_data(vtk_file, header_type, fields["type"],
                              data_len, data_buf, decode_buf, offset);
          mesh.e2n_dsp.resize(mesh.e2n_cnt.size()); mesh.e2n_dsp[0] = 0;
          for(size_t i=1; i<mesh.e2n_dsp.size(); i++)
            mesh.e2n_dsp[i] = offset[i-1];
          for(size_t i=0; i<mesh.e2n_cnt.size(); i++)
            mesh.e2n_cnt[i] = offset[i] - mesh.e2n_dsp[i];
          prg.next();
        }
        else if((fields["tag"].compare("DataArray") == 0) &&
                (fields["Name"].compare("types") == 0))
        {
          mt_vector<unsigned char> types;

          data_len = read_DataArray_header(vtk_file, header_type);
          if(data_len / DataArray_type_size(fields["type"]) != mesh.e2n_cnt.size())
            treat_file_error("vtu read error: size of types DataArray"
                             "does not match number of elements.", vtk_file);

          read_DataArray_data(vtk_file, header_type, fields["type"],
                              data_len, data_buf, decode_buf, types);
          for(size_t i=0; i<types.size(); i++) {
            switch(types[i]) {
              case 3:
                mesh.etype[i] = Line; // Lines are encoded as index 3
                break;
              case 5:
                mesh.etype[i] = Tri; // Triangles are encoded as index 5
                break;
              case 10:
                mesh.etype[i] = Tetra; // Tetras are encoded as index 10
                break;
              case 9:
                mesh.etype[i] = Quad; // Quads are encoded as index 9
                break;
              case 14:
                mesh.etype[i] = Pyramid; // Pyramids are encoded as index 14
                break;
              case 13:
                mesh.etype[i] = Prism; // Prisms are encoded as index 13
                break;
              case 12:
                mesh.etype[i] = Hexa; // Hexahedras are encoded as index 12
                break;
            }
          }
          prg.next();
        }
        else if(fields["tag"].compare("/Cells") == 0)
          break;
      }
    }
    else if(fields["tag"].compare("CellData") == 0) {
      mt_vector<mt_real> fib, she;

      while(ptr) {
        ptr = next_line(buffer, bufsize, vtk_file);
        parse_xml_line(ptr, fields);

        if((fields["tag"].compare("DataArray") == 0) &&
           (fields["Name"].compare("elemTag") == 0))
        {
          data_len = read_DataArray_header(vtk_file, header_type);
          if(data_len / DataArray_type_size(fields["type"]) != mesh.etags.size())
            treat_file_error("vtu read error: size of tags DataArray "
                             "does not match number of elements.", vtk_file);

          read_DataArray_data(vtk_file, header_type, fields["type"],
                              data_len, data_buf, decode_buf, mesh.etags);
          prg.next();
        }
        else if((fields["tag"].compare("DataArray") == 0) &&
                (fields["Name"].compare("fiber") == 0))
        {
          data_len = read_DataArray_header(vtk_file, header_type);
          if(data_len / DataArray_type_size(fields["type"]) / 3 != mesh.etags.size())
            treat_file_error("vtu read error: size of fiber DataArray "
                             "does not match number of elements.", vtk_file);

          read_DataArray_data(vtk_file, header_type, fields["type"],
                              data_len, data_buf, decode_buf, fib);
        }
        else if((fields["tag"].compare("DataArray") == 0) &&
                (fields["Name"].compare("sheet") == 0))
        {
          data_len = read_DataArray_header(vtk_file, header_type);
          if(data_len / DataArray_type_size(fields["type"]) / 3 != mesh.etags.size())
            treat_file_error("vtu read error: size of sheet DataArray "
                             "does not match number of elements.", vtk_file);

          read_DataArray_data(vtk_file, header_type, fields["type"],
                              data_len, data_buf, decode_buf, she);
        }
        else if(fields["tag"].compare("/CellData") == 0)
          break;
      }

      size_t numele = mesh.e2n_cnt.size();
      if(fib.size() == numele*3 && she.size() == numele*3) {
        mesh.lon.resize(mesh.e2n_cnt.size()*6);
        for(size_t i=0; i<mesh.e2n_cnt.size(); i++) {
          mesh.lon[i*6+0] = fib[i*3+0];
          mesh.lon[i*6+1] = fib[i*3+1];
          mesh.lon[i*6+2] = fib[i*3+2];
          mesh.lon[i*6+3] = she[i*3+0];
          mesh.lon[i*6+4] = she[i*3+1];
          mesh.lon[i*6+5] = she[i*3+2];
        }
      }
      else if(fib.size() == numele*3) {
        mesh.lon.resize(mesh.e2n_cnt.size()*3);
        for(size_t i=0; i<mesh.e2n_cnt.size(); i++) {
          mesh.lon[i*3+0] = fib[i*3+0];
          mesh.lon[i*3+1] = fib[i*3+1];
          mesh.lon[i*3+2] = fib[i*3+2];
        }
      }
      else
        fprintf(stderr, "Warning, no fibers provided!\n");

      prg.next();
    }
    else if(fields["tag"].compare("/Piece") == 0)
      finished = true;

    ptr = next_line(buffer, bufsize, vtk_file);
  }

  return finished;
}

void readXMLVTKmesh(mt_meshdata & mesh, std::string file)
{
  const int bufsize = 2048;
  char buffer[bufsize], *ptr;

  FILE* vtk_file = fopen(file.c_str(), MT_FOPEN_READ);
  if(vtk_file == NULL) treat_file_open_error(file, errno);
  MT_MAP<std::string, std::string> fields;

  ptr = next_line(buffer, bufsize, vtk_file);
  parse_xml_line(ptr, fields);

  // if the first line was the xml specifier expression, read another line
  if(fields["tag"].compare("?xml") == 0) {
    ptr = next_line(buffer, bufsize, vtk_file);
    parse_xml_line(ptr, fields);
  }

  // we only read UnstructuredGrid files
  if(!((fields["tag"].compare("VTKFile") == 0) && 
       (fields["type"].compare("UnstructuredGrid") == 0)))
    treat_file_error("vtu read error: file is not an UnstructuredGrid vtk file.", vtk_file);

  std::string header_type = "UInt64";
  if(fields.count("header_type") == 1)
    header_type = fields["header_type"];

  PROGRESS<short> prg(7, "Reading binary vtu file: ");

  bool finished = false;
  ptr = next_line(buffer, bufsize, vtk_file);
  while(ptr && !finished)
  {
    parse_xml_line(ptr, fields);
    if(fields["tag"].compare("Piece") == 0) {
      int numnod = atoi(fields["NumberOfPoints"].c_str());
      int numele = atoi(fields["NumberOfCells"].c_str());
      mesh.e2n_cnt.resize(numele);
      mesh.etags.resize(numele);
      mesh.etype.resize(numele);
      mesh.xyz.resize(numnod*3);
      finished = process_xmlvtk_piece(mesh, vtk_file, header_type, prg);
    }
    ptr = next_line(buffer, bufsize, vtk_file);
  }

  fclose(vtk_file);
  prg.finish();

  mt_vector<mt_int> e2n_con(mesh.e2n_con);
  perm_con_vtk_to_carp(mesh, e2n_con);
}

void perm_con_carp_to_vtk(const mt_meshdata & mesh, mt_vector<mt_int> & e2n_con)
{
  e2n_con.resize(mesh.e2n_con.size());

  const mt_int * inp_elem = mesh.e2n_con.data();
  mt_int *       out_elem = e2n_con.data();

  for(size_t eidx=0; eidx < mesh.e2n_cnt.size(); eidx++) {
    int nnodes = mesh.e2n_cnt[eidx];

    switch(mesh.etype[eidx]) {
      case Hexa:
        out_elem[0] = inp_elem[5];
        out_elem[1] = inp_elem[6];
        out_elem[2] = inp_elem[7];
        out_elem[3] = inp_elem[4];
        out_elem[4] = inp_elem[3];
        out_elem[5] = inp_elem[2];
        out_elem[6] = inp_elem[1];
        out_elem[7] = inp_elem[0];
        break;

      case Prism:
        out_elem[0] = inp_elem[0];
        out_elem[1] = inp_elem[1];
        out_elem[2] = inp_elem[2];
        out_elem[3] = inp_elem[3];
        out_elem[4] = inp_elem[5];
        out_elem[5] = inp_elem[4];
        break;

      default:
        memcpy(out_elem, inp_elem, nnodes*sizeof(mt_int));
        break;
    }
    out_elem += nnodes;
    inp_elem += nnodes;
  }
}


void perm_con_vtk_to_carp(mt_meshdata & mesh, const mt_vector<mt_int> & e2n_con)
{
  const mt_int * inp_elem = e2n_con.data();
  mt_int *       out_elem = mesh.e2n_con.data();

  for(size_t eidx=0; eidx < mesh.e2n_cnt.size(); eidx++) {
    int nnodes = mesh.e2n_cnt[eidx];

    switch(mesh.etype[eidx]) {
      case Hexa:
        out_elem[0] = inp_elem[7];
        out_elem[1] = inp_elem[6];
        out_elem[2] = inp_elem[5];
        out_elem[3] = inp_elem[4];
        out_elem[4] = inp_elem[3];
        out_elem[5] = inp_elem[0];
        out_elem[6] = inp_elem[1];
        out_elem[7] = inp_elem[2];
        break;

      case Prism:
        out_elem[0] = inp_elem[0];
        out_elem[1] = inp_elem[1];
        out_elem[2] = inp_elem[2];
        out_elem[3] = inp_elem[3];
        out_elem[4] = inp_elem[5];
        out_elem[5] = inp_elem[4];
        break;

      default:
        memcpy(out_elem, inp_elem, nnodes*sizeof(mt_int));
        break;
    }
    out_elem += nnodes;
    inp_elem += nnodes;
  }
}

