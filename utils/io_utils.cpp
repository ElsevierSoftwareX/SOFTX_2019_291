/**
* @file io_utils.cpp
* @brief General IO utils.
* @author Aurel Neic
* @version
* @date 2017-08-04
*/

#include "mt_utils.h"
#include "asciireader.hpp"

#include <unistd.h>  // needed for access()

elem_t getElemTypeID(char *eletype)
{
  elem_t ret;

  if ( !strcmp( eletype, "Tt" ) ) {
    ret = Tetra;
  } else if ( !strcmp( eletype, "Hx" ) ) {
    ret = Hexa;
//  } else if ( !strcmp( eletype, "Oc" ) ) {
//    ret = Octa;
  } else if ( !strcmp( eletype, "Py" ) ) {
    ret = Pyramid;
  } else if ( !strcmp( eletype, "Pr" ) ) {
    ret = Prism;
  } else if ( !strcmp( eletype, "Qd" ) ) {
    ret = Quad;
  } else if ( !strcmp( eletype, "Tr" ) ) {
    ret = Tri;
  } else {
    ret = Line;
  }
  return ret;
}


void treat_file_open_error(std::string file, int errnum, bool do_exit)
{
  fprintf(stderr, "An IO error occured when opening file %s.\n%s\n\n", file.c_str(), strerror(errnum));
  if(do_exit)
    exit(1);
}

void treat_file_read_error(std::string file, FILE* fd)
{
  fprintf(stderr, "An error occured when reading file %s.\n\n", file.c_str());
  fclose(fd);
  exit(1);
}

void treat_file_error(const char* msg, FILE* fd)
{
  fprintf(stderr, "%s\n", msg);
  fclose(fd);
  exit(1);
}

bool file_exists(std::string filename)
{
  return (access(filename.c_str(), F_OK) == 0);
}

size_t file_size(FILE* fd)
{
  size_t oldpos = ftell(fd);
  fseek(fd, 0L, SEEK_END);
  size_t sz = ftell(fd);
  fseek(fd, oldpos, SEEK_SET);

  return sz;
}

size_t file_size(const char* file)
{
  FILE* fd = fopen(file, MT_FOPEN_READ);
  if(fd == NULL) return 0;

  size_t ret = file_size(fd);
  fclose(fd);

  return ret;
}


void generate_default_fib(size_t numelem, mt_vector<mt_real> & lon)
{
  lon.resize(numelem*3);
  for(size_t i=0; i<numelem; i++) {
    lon[i*3+0] = 1.0;
    lon[i*3+1] = 0.0;
    lon[i*3+2] = 0.0;
  }
}

void writeElements(mt_meshdata & mesh, std::string file)
{
  FILE* ele_file = fopen(file.c_str(), MT_FOPEN_WRITE);
  if(ele_file == NULL) treat_file_open_error(file, errno);

  size_t numele = mesh.e2n_cnt.size();
  fprintf(ele_file, "%lu\n", numele);
  PROGRESS<size_t> prg(numele, "Writing elements in text CARP format: ");

  for(size_t i=0, idx = 0; i<numele; i++)
  {
    prg.next();

    std::string etyp;
    switch(mesh.etype[i])
    {
      case Line:
        etyp = "Ln";
        break;

      case Tri:
        etyp = "Tr";
        break;

      case Quad:
        etyp = "Qd";
        break;

      case Tetra:
        etyp = "Tt";
        break;

      case Pyramid:
        etyp = "Py";
        break;

      case Prism:
        etyp = "Pr";
        break;

      case Hexa:
        etyp = "Hx";
        break;

      default:
        fprintf(stderr, "Error: Unsupported element type!\n");
        fclose(ele_file);
        exit(1);
    }
    fprintf(ele_file, "%s ", etyp.c_str());
    for(mt_int j=0; j<mesh.e2n_cnt[i]; j++, idx++)
      fprintf(ele_file, "%lu ", (unsigned long)mesh.e2n_con[idx]);

    fprintf(ele_file, "%d\n", (int)mesh.etags[i]);
  }
  prg.finish();

  fclose(ele_file);
}

void writeElementsBinary(mt_meshdata & mesh, const std::string file)
{
  FILE* ele_file = fopen(file.c_str(), MT_FOPEN_WRITE);
  if(ele_file == NULL) treat_file_open_error(file, errno);

  size_t numele = mesh.e2n_cnt.size();

  // we need to write the CARP header
  char header[BIN_HDR_SIZE]; memset(header, 0, sizeof(header));
  int checksum = 666;
  int n[10];

  // first elements, second endiannes, third checksum
  sprintf(header, "%lu %d %d", numele, MT_ENDIANNESS, checksum);
  fwrite(header, sizeof(char), BIN_HDR_SIZE, ele_file);
  PROGRESS<size_t> prg(numele, "Writing elements in binary CARP format: ");

  mt_int* elem = mesh.e2n_con.data();
  for(size_t i=0; i<numele; i++)
  {
    prg.next();

    mt_int nodes = mesh.e2n_cnt[i];

    // element type
    n[0] = mesh.etype[i];
    fwrite(n, sizeof(int), 1, ele_file);

    // element connectivity
    for(mt_int j=0; j<nodes; j++) n[j] = elem[j];
    fwrite(n, sizeof(int), nodes, ele_file);

    // element tag
    n[0] = mesh.etags[i];
    fwrite(n, sizeof(int), 1, ele_file);

    elem += nodes;
  }
  prg.finish();

  fclose(ele_file);
}

void readElements(mt_meshdata & mesh, std::string file)
{
  const int bufsize = 2048;
  char buffer[bufsize];
  const int max_elem_size = 9;

  unsigned long int numele = 0;

  asciireader reader;
  reader.read_file(file);

  if(reader.get_line(0, buffer, bufsize))
    sscanf(buffer, "%lu", &numele);

  mesh.e2n_cnt.resize(numele);
  mesh.etags.assign(size_t(numele), mt_int(0));
  mesh.etype.resize(numele);
  mesh.e2n_con.assign(size_t(numele*max_elem_size), mt_int(-1));
  PROGRESS<size_t> prg(numele, "Reading elements in text CARP format: ");

  #ifdef OPENMP
  #pragma omp parallel private(buffer)
  #endif
  {
    char etype_str[16];
    long int n[max_elem_size]; memset(n, 0, max_elem_size*sizeof(long int));
    int num_read = 0;
    const int check_intervall = 100;

    #ifdef OPENMP
    #pragma omp for schedule(dynamic, 100)
    #endif
    for(size_t i=0; i<numele; i++)
    {
      if(reader.get_line(i+1, buffer, bufsize))
      {
        // read the element
        //int numread =
        int nread = sscanf(buffer, "%s %ld %ld %ld %ld %ld %ld %ld %ld %ld", etype_str,
                                   n, n+1, n+2, n+3, n+4, n+5, n+6, n+7, n+8);
        mesh.etype[i] = getElemTypeID(etype_str);
        mt_int nodes;
        switch(mesh.etype[i])
        {
          case Line:
            nodes = 2;
            break;

          case Tri:
            nodes = 3;
            break;

          case Quad:
            nodes = 4;
            break;

          case Tetra:
            nodes = 4;
            break;

          case Pyramid:
            nodes = 5;
            break;

          case Prism:
            nodes = 6;
            break;

          case Hexa:
            nodes = 8;
            break;

          default:
            fprintf(stderr, "Error: Unsupported element type!\n");
            exit(1);
        }
        // set element size
        mesh.e2n_cnt[i] = nodes;

        // copy the element connectivity
        for(int j=0; j<nodes; j++) mesh.e2n_con[i*max_elem_size+j] = n[j];

        if(nread > nodes)
          mesh.etags[i] = n[nodes];

        num_read++;
      }

      if(num_read == check_intervall) {
        #ifdef OPENMP
        #pragma omp critical
        #endif
        {
          prg.increment(num_read);
        }
        num_read = 0;
      }
    }
  }

  prg.finish();

  reduce_by_flag(mesh.e2n_con, mt_int(-1));
  mesh.e2n_con.reallocate();
}

void readElementsBinary(mt_meshdata & mesh, std::string file)
{
  FILE* ele_file = fopen(file.c_str(), MT_FOPEN_READ);
  if(ele_file == NULL) treat_file_open_error(file, errno);

  int max_elem_size = 8;
  unsigned long int numele, numentries = 0;
  int n[max_elem_size]; memset(n, 0, max_elem_size*sizeof(int));

  char header[BIN_HDR_SIZE]; memset(header, 0, sizeof(header));
  fread(header, sizeof(char), BIN_HDR_SIZE, ele_file);
  sscanf(header, "%lu", &numele);

  mesh.e2n_cnt.resize(numele);
  mesh.etags.resize(numele);
  mesh.etype.resize(numele);
  mesh.e2n_con.assign(numele*max_elem_size, -1);
  int etbuff;
  mt_int* elem = mesh.e2n_con.data();
  PROGRESS<size_t> prg(numele, "Reading elements in binary CARP format: ");

  for(size_t i=0; i<numele; i++)
  {
    prg.next();

    fread(&etbuff, sizeof(int), 1, ele_file);
    mesh.etype[i] = (elem_t)etbuff;
    int nodes = 0;
    switch(mesh.etype[i])
    {
      case Line:
        nodes = 2;
        break;

      case Tri:
        nodes = 3;
        break;

      case Tetra:
        nodes = 4;
        break;

      case Quad:
        nodes = 4;
        break;

      case Pyramid:
        nodes = 5;
        break;

      case Prism:
        nodes = 6;
        break;

      case Hexa:
        nodes = 8;
        break;
    }
    mesh.e2n_cnt[i] = nodes;

    int nread = fread(n, sizeof(int), nodes, ele_file);
    if(nread != nodes)
      fprintf(stderr, "%s error: read only %d nodes instead of %d.\n", __func__, nread, nodes);

    // copy the element connectivity
    for(int j=0; j<nodes; j++) elem[j] = n[j];
    elem += nodes;
    numentries += nodes;

    // read tag
    nread = fread(n, sizeof(int), 1, ele_file);
    mesh.etags[i] = n[0];

    if(nread != 1)
      fprintf(stderr, "%s error: could not read fiber.\n", __func__);
  }
  prg.finish();

  mesh.e2n_con.resize(numentries);
  mesh.e2n_con.reallocate();

  fclose(ele_file);
}

void readElements_general(mt_meshdata & mesh,
                          std::string basename)
{
  std::string binname = basename + CARPBIN_ELEM_EXT;
  std::string txtname = basename + CARPTXT_ELEM_EXT;

  if(file_exists(binname))
    readElementsBinary(mesh, binname);
  else
    readElements(mesh, txtname);
}

void readElementTags_general(mt_vector<mt_int> & etags, std::string file)
{
  if(!file_exists(file)) {
    printf("Element-tags file %s does not exist!", file.c_str());
    return;
  }

  if (endswith(file, BIN_TAGS_EXT))
    binary_read(etags, file);
  else if (endswith(file, TXT_TAGS_EXT))
    read_vector_ascii(etags, file, false);
  else
    printf("Element-tags file %s has unknown extension!", file.c_str());
}

size_t readNumPoints(std::string file)
{
  FILE* pts_file = fopen(file.c_str(), MT_FOPEN_READ);
  if(pts_file == NULL) treat_file_open_error(file, errno);

  size_t numpts = 0;
  int bufsize = 2056;
  char buffer[bufsize];
  char *ptr = fgets( buffer, bufsize, pts_file);
  if(ptr != NULL) sscanf(buffer, "%lu", &numpts);

  fclose(pts_file);
  return numpts;
}

void readPoints(mt_vector<mt_real> & xyz, std::string file)
{
  int bufsize = 2048;
  char buffer[bufsize];
  unsigned long int numpts = 0;
  asciireader reader;
  reader.read_file(file);

  if(reader.get_line(0, buffer, bufsize))
    sscanf(buffer, "%lu", &numpts);

  xyz.resize(numpts*3);

  PROGRESS<size_t> prg(numpts, "Reading points in text CARP format: ");

  #ifdef OPENMP
  #pragma omp parallel private(buffer)
  #endif
  {
    float pts[3];
    int num_read = 0;
    const int check_intervall = 100;

    #ifdef OPENMP
    #pragma omp for schedule(dynamic, 100)
    #endif
    for(size_t i=0; i<numpts; i++)
    {

      if(reader.get_line(i+1, buffer, bufsize))
      {
        sscanf(buffer, "%f %f %f", pts, pts+1, pts+2);
        xyz[i*3+0] = pts[0]; xyz[i*3+1] = pts[1]; xyz[i*3+2] = pts[2];
        num_read++;
      }

      if(num_read == check_intervall) {
        #ifdef OPENMP
        #pragma omp critical
        #endif
        {
          prg.increment(num_read);
        }
        num_read = 0;
      }
    }
  }

  prg.finish();
}

void readPointsBinary(mt_vector<mt_real> & xyz, std::string file)
{
  FILE* pts_file = fopen(file.c_str(), MT_FOPEN_READ);
  if(pts_file == NULL) treat_file_open_error(file, errno);

  unsigned long int numpts;
  char header[BIN_HDR_SIZE]; memset(header, 0, sizeof(header));
  fread(header, sizeof(char), BIN_HDR_SIZE, pts_file);
  sscanf(header, "%lu", &numpts);

  xyz.resize(numpts*3);
  mt_real* wp = xyz.data();
  float pts[3];
  PROGRESS<size_t> prg(size_t(numpts), "Reading points in binary CARP format: ");

  for(unsigned long int i=0; i<numpts; i++)
  {
    prg.next();

    fread(pts, sizeof(float), 3, pts_file);
    wp[0] = pts[0]; wp[1] = pts[1]; wp[2] = pts[2];
    wp += 3;
  }
  prg.finish();

  fclose(pts_file);
}

void readPoints_general(mt_vector<mt_real> & xyz, std::string basename)
{
  std::string binname = basename + CARPBIN_PTS_EXT;
  std::string txtname = basename + CARPTXT_PTS_EXT;

  if(file_exists(binname))
    readPointsBinary(xyz, binname);
  else
    readPoints(xyz, txtname);
}

void readUVCPoints(mt_vector<mt_real> & xyz_lv,
                   mt_vector<mt_real> & xyz_rv,
                   mt_vector<mt_int> & pos_lv,
                   mt_vector<mt_int> & pos_rv,
                   mt_vector<mixed_tuple<short, size_t> > & idx_to_uvc_point,
                   std::string filename)
{
  int bufsize = 2048;
  char buffer[bufsize];
  unsigned long int numpts = 0;
    asciireader reader;
  reader.read_file(filename);

  if(reader.get_line(0, buffer, bufsize))
    sscanf(buffer, "%lu", &numpts);
  pos_lv.resize(0);
  pos_lv.reserve(2*numpts);
  pos_rv.resize(0);
  pos_rv.reserve(2*numpts);
  xyz_lv.resize(0);
  xyz_lv.reserve(6*numpts);
  xyz_rv.resize(0);
  xyz_rv.reserve(6*numpts);
  idx_to_uvc_point.resize(numpts);

  PROGRESS<size_t> prg(numpts, "Reading UVC points: ");
  float pts[3];
  float v;
  int num_read = 0;
  const int check_intervall = 100;
  const mt_real brsz = 2.0 * MT_PI;
  const mt_real brlim0 = -MT_PI + brsz * 0.05;
  const mt_real brlim1 = MT_PI - brsz * 0.05;
  for(size_t i=0; i<numpts; i++)
  {
    if(reader.get_line(i+1, buffer, bufsize))
    {
      sscanf(buffer, "%f %f %f %f", pts, pts+1, pts+2, &v);
      if(v > 0) {
        xyz_lv.push_back(pts[1]);
        xyz_lv.push_back(pts[2]);
        xyz_lv.push_back(pts[0]);
        pos_lv.push_back(i);
        idx_to_uvc_point[i].v1 = 1;
        idx_to_uvc_point[i].v2 = pos_lv.size()-1;
        //check for proximity to branch cut
        if(pts[2] < brlim0 || pts[2] > brlim1)
        {
          xyz_lv.push_back(pts[1]);
          const mt_real fac = pts[2] < brlim0 ? 1. : -1.;
          xyz_lv.push_back(pts[2] + fac * brsz);
          xyz_lv.push_back(pts[0]);
          pos_lv.push_back(i);
        } 
      }
      else
      {
        xyz_rv.push_back(pts[1]);
        xyz_rv.push_back(pts[2]);
        xyz_rv.push_back(pts[0]);
        pos_rv.push_back(i);
        idx_to_uvc_point[i].v1 = -1;
        idx_to_uvc_point[i].v2 = pos_rv.size()-1;
        if(pts[2] < brlim0 || pts[2] > brlim1)
        {
          xyz_rv.push_back(pts[1]);
          const mt_real fac = pts[2] < brlim0 ? 1. : -1.;
          xyz_rv.push_back(pts[2] + fac * brsz);
          xyz_rv.push_back(pts[0]);
          pos_rv.push_back(i);
        } 
      }
      num_read++;
    }
    if(num_read == check_intervall) {
      prg.increment(num_read);
      num_read = 0;
    }
  }
  prg.finish();
}

void writePoints(mt_vector<mt_real> & xyz, std::string file)
{
  FILE* pts_file = fopen(file.c_str(), MT_FOPEN_WRITE);
  if(pts_file == NULL) treat_file_open_error(file, errno);

  unsigned long int numpts = xyz.size() / 3;
  fprintf(pts_file, "%lu\n", numpts);

  mt_real* wp = xyz.data();
  float pts[3];
  PROGRESS<size_t> prg(numpts, "Writing points in text CARP format: ");

  for(size_t i=0; i<numpts; i++)
  {
    prg.next();

    pts[0] = wp[0]; pts[1] = wp[1]; pts[2] = wp[2];
    fprintf(pts_file, "%f %f %f\n", pts[0], pts[1], pts[2]);
    wp += 3;
  }
  prg.finish();

  fclose(pts_file);
}

void writePointsBinary(mt_vector<mt_real> & xyz, std::string file)
{
  FILE* pts_file = fopen(file.c_str(), MT_FOPEN_WRITE);
  if(pts_file == NULL) treat_file_open_error(file, errno);

  // we need to write the CARP header
  unsigned long int numpts = xyz.size() / 3;
  char header[BIN_HDR_SIZE]; memset(header, 0, sizeof(header));
  int checksum = 666;

  sprintf(header, "%lu %d %d", numpts, MT_ENDIANNESS, checksum);
  fwrite(header, sizeof(char), BIN_HDR_SIZE, pts_file);

  mt_real* wp = xyz.data();
  float pts[3];
  PROGRESS<size_t> prg(numpts, "Writing points in binary CARP format: ");

  for(unsigned long int i=0; i<numpts; i++)
  {
    prg.next();

    pts[0] = wp[0]; pts[1] = wp[1]; pts[2] = wp[2];
    fwrite(pts, sizeof(float), 3, pts_file);
    wp += 3;
  }
  prg.finish();

  fclose(pts_file);
}

int readFibers(mt_vector<mt_real> & lon, size_t numelem, std::string file)
{
  int bufsize = 2048;
  char buffer[bufsize];
  int numfibres = 0;

  asciireader reader;
  bool exit_on_error = false;
  reader.read_file(file, exit_on_error);

  // we may not have fibers
  if(reader.num_lines() == 0) return 0;

  // the increment we need from fiber number to fibre file line
  int line_inc = 1;

  if(reader.get_line(0, buffer, bufsize)) {
    float pts[6];
    size_t nread = sscanf(buffer, "%f %f %f %f %f %f",
                          pts, pts+1, pts+2, pts+3, pts+4, pts+5);
    switch(nread) {
      case 1: numfibres = pts[0]; break;
      case 3: numfibres = 1; line_inc = 0; break;
      case 6: numfibres = 2; line_inc = 0; break;
      default:
        fprintf(stderr, "%s error parsing first line of %s ! Aborting!\n", __func__,
                file.c_str());
        exit(EXIT_FAILURE);
    }
  }

  lon.resize(numelem*numfibres*3);

  char msgstring[1024];
  sprintf(msgstring, "Reading fibers (%d directions) in text CARP format: ", numfibres);
  PROGRESS<size_t> prg(numelem, msgstring);
  std::string errdeco = std::string(80, '~');
  int global_err = 0;

  #ifdef OPENMP
  #pragma omp parallel private(buffer)
  #endif
  {
    float pts[6];
    bool readerror = false;
    int num_read = 0;
    const int check_intervall = 100;

    #ifdef OPENMP
    #pragma omp for schedule(dynamic, 100)
    #endif
    for(size_t i=0; i<numelem; i++)
    {
      readerror = !reader.get_line(i+line_inc, buffer, bufsize);
      if(!readerror) {
        if(numfibres == 2) {
          sscanf(buffer, "%f %f %f %f %f %f", pts, pts+1, pts+2, pts+3, pts+4, pts+5);
          lon[i*6+0] = pts[0]; lon[i*6+1] = pts[1]; lon[i*6+2] = pts[2];
          lon[i*6+3] = pts[3]; lon[i*6+4] = pts[4]; lon[i*6+5] = pts[5];
        }
        else {
          sscanf(buffer, "%f %f %f", pts, pts+1, pts+2);
          lon[i*3+0] = pts[0]; lon[i*3+1] = pts[1]; lon[i*3+2] = pts[2];
        }
        num_read++;
      }

      if(num_read == check_intervall) {
        #ifdef OPENMP
        #pragma omp critical
        #endif
        {
          global_err += int(readerror);
          prg.increment(num_read);
        }
        num_read = 0;
      }
    }
  }
  prg.finish();

  if (global_err)
  {
    std::cerr << std::endl << errdeco << std::endl;
    std::cerr << "ATTENTION: fibre number does not match element number!" << std::endl;
    std::cerr << errdeco << std::endl;
  }

  return numfibres;
}

int readFibersBinary(mt_vector<mt_real> & lon, size_t numelem, std::string file)
{
  FILE* _file = fopen(file.c_str(), MT_FOPEN_READ);
  if(_file == NULL) treat_file_open_error(file, errno, false);

  int numfibres;

  char header[BIN_HDR_SIZE]; memset(header, 0, sizeof(header));
  fread(header, sizeof(char), BIN_HDR_SIZE, _file);
  sscanf(header, "%d", &numfibres);

  bool errmsg = false;
  std::string errdeco = std::string(80, '~');
  if(numfibres == 2)
  {
    lon.resize(numelem*6);
    mt_real* wp = lon.data();
    float pts[6];
    PROGRESS<size_t> prg(numelem, "Reading fibers (2 directions) in binary CARP format: ");

    for(size_t i=0; i<numelem; i++)
    {
      prg.next();
      if (fread(pts, sizeof(float), 6, _file) == 6)
      {
        wp[0] = pts[0]; wp[1] = pts[1]; wp[2] = pts[2];
        wp[3] = pts[3]; wp[4] = pts[4]; wp[5] = pts[5];
        wp += 6;
      }
      else if (errmsg == false)
      {
        std::cerr << std::endl << errdeco << std::endl;
        std::cerr << std::string() << "WARNING: only " << i << " / " << numelem << " fibre directions read!" << std::endl;
        std::cerr << errdeco << std::endl;
        errmsg = true;
      }
    }
    prg.finish();
  }
  else {
    lon.resize(numelem*3);
    mt_real* wp = lon.data();
    float pts[3];
    PROGRESS<size_t> prg(numelem, "Reading fibers (1 directions) in binary CARP format: ");

    for(size_t i=0; i<numelem; i++)
    {
      prg.next();

      if (fread(pts, sizeof(float), 3, _file) == 3)
      {
        wp[0] = pts[0]; wp[1] = pts[1]; wp[2] = pts[2];
        wp += 3;
      }
      else if (errmsg == false)
      {
        std::cerr << std::endl << errdeco << std::endl;
        std::cerr << std::string() << "WARNING: only " << i << " / " << numelem << " fibre directions read!" << std::endl;
        std::cerr << errdeco << std::endl;
        errmsg = true;
      }
    }
    prg.finish();
  }

  fclose(_file);
  return numfibres;
}

void readFibers_general(mt_vector<mt_real> & lon, size_t numelem, std::string basename)
{
  std::string binname = basename + CARPBIN_LON_EXT;
  std::string txtname = basename + CARPTXT_LON_EXT;

  if(file_exists(binname))
    readFibersBinary(lon, numelem, binname);
  else
    readFibers(lon, numelem, txtname);
}

void writeFibers(mt_vector<mt_real> & lon, size_t numelem, std::string file)
{
  FILE* _file = fopen(file.c_str(), MT_FOPEN_WRITE);
  if(_file == NULL) treat_file_open_error(file, errno);

  int numfibres = (lon.size() / 6) == numelem ? 2 : 1;
  fprintf(_file, "%d\n", numfibres);

  if(numfibres == 2)
  {
    mt_real* wp = lon.data();
    float pts[6];
    PROGRESS<size_t> prg(numelem, "Writing fibers (2 directions) in text CARP format: ");

    for(size_t i=0; i<numelem; i++)
    {
      prg.next();

      pts[0] = wp[0]; pts[1] = wp[1]; pts[2] = wp[2];
      pts[3] = wp[3]; pts[4] = wp[4]; pts[5] = wp[5];
      fprintf(_file, "%f %f %f %f %f %f\n", pts[0], pts[1], pts[2], pts[3], pts[4], pts[5]);
      wp += 6;
    }
    prg.finish();
  }
  else {
    mt_real* wp = lon.data();
    float pts[3];
    PROGRESS<size_t> prg(numelem, "Writing fibers (1 direction) in text CARP format: ");

    for(size_t i=0; i<numelem; i++)
    {
      prg.next();

      pts[0] = wp[0]; pts[1] = wp[1]; pts[2] = wp[2];
      fprintf(_file, "%f %f %f\n", pts[0], pts[1], pts[2]);
      wp += 3;
    }
    prg.finish();
  }

  fclose(_file);
}

void writeFibersBinary(mt_vector<mt_real> & lon, size_t numelem, std::string file)
{
  FILE* _file = fopen(file.c_str(), MT_FOPEN_WRITE);
  if(_file == NULL) treat_file_open_error(file, errno);

  int numfibres = (lon.size() / 6) == numelem ? 2 : 1;

  // we need to write the CARP header
  char header[BIN_HDR_SIZE]; memset(header, 0, sizeof(header));
  int checksum = 666;

  // first numaxes, second numelems, third endiannes, last checksum
  sprintf(header, "%d %lu %d %d", numfibres, (unsigned long int)numelem, MT_ENDIANNESS, checksum);
  fwrite(header, sizeof(char), BIN_HDR_SIZE, _file);

  if(numfibres == 2)
  {
    mt_real* wp = lon.data();
    float pts[6];
    PROGRESS<size_t> prg(numelem, "Writing fibers (2 directions) in binary CARP format: ");

    for(size_t i=0; i<numelem; i++)
    {
      prg.next();

      pts[0] = wp[0]; pts[1] = wp[1]; pts[2] = wp[2];
      pts[3] = wp[3]; pts[4] = wp[4]; pts[5] = wp[5];
      fwrite(&pts, sizeof(float), 6, _file);
      wp += 6;
    }
    prg.finish();
  }
  else {
    mt_real* wp = lon.data();
    float pts[3];
    PROGRESS<size_t> prg(numelem, "Writing fibers (1 direction) in binary CARP format: ");

    for(size_t i=0; i<numelem; i++)
    {
      prg.next();

      pts[0] = wp[0]; pts[1] = wp[1]; pts[2] = wp[2];
      fwrite(&pts, sizeof(float), 3, _file);
      wp += 3;
    }
    prg.finish();
  }

  fclose(_file);
}

void write_surf(const mt_vector<mt_int> & cnt, const mt_vector<mt_int> & con,
               const std::string filename)
{
  FILE* file = fopen(filename.c_str(), MT_FOPEN_WRITE);
  if(file == NULL) treat_file_open_error(filename, errno);

  size_t numsurf = cnt.size();
  fprintf(file, "%lu\n", numsurf);

  for(size_t i=0,k=0; i<numsurf; i++) {
    switch(cnt[i]) {
      case 3:
        fprintf(file, "Tr %ld %ld %ld\n", long(con[k]), long(con[k+1]), long(con[k+2]));
        k += 3;
        break;
      case 4:
        fprintf(file, "Qd %ld %ld %ld %ld\n",
                long(con[k]), long(con[k+1]), long(con[k+2]), long(con[k+3]));
        k += 4;
        break;
      default: break;
    }
  }

  fclose(file);
}

void read_vtx(mt_vector<mt_int> & nod, std::string filename)
{
  FILE* file = fopen(filename.c_str(), MT_FOPEN_READ);
  if(file == NULL) treat_file_open_error(filename, errno);

  int bufsize = 2048;
  char buffer[bufsize];
  int numpts = 0;

  // Read just the number of points from the header
  char *fp = fgets(buffer, bufsize, file);
  if(fp != NULL) sscanf(buffer, "%d", &numpts);

  // Discard next line
  fgets(buffer, bufsize, file);

  // Allocate
  nod.resize(numpts);
  long nb;

  for(int i=0; i<numpts; i++) {
    fp = fgets(buffer, bufsize, file);
    if (fp != NULL) {
      sscanf(buffer, "%ld", &nb);
      nod[i] = nb;
    }
  }

  fclose(file);
}

void write_vtx(const mt_vector<mt_int> & nod, const std::string filename,
              bool write_header)
{
  FILE* file = fopen(filename.c_str(), MT_FOPEN_WRITE);
  if(file == NULL) treat_file_open_error(filename, errno);

  // Header
  if(write_header) {
    fprintf(file, "%lu\n", (long unsigned int) nod.size());
    fprintf(file, "extra\n");
  }

  // Content
  for(size_t i=0; i<nod.size(); i++)
    fprintf(file, "%ld\n", long(nod[i]));

  fclose(file);
}

void read_nbc(nbc_data & nbc, mt_vector<mt_int> & con, std::string filename)
{
  FILE* file = fopen(filename.c_str(), MT_FOPEN_READ);
  if(file == NULL) treat_file_open_error(filename, errno);

  int bufsize = 2048;
  char buffer[bufsize];
  int numface = 0;

  // Read just the number of faces from the header
  char *ptr = fgets(buffer, bufsize, file);
  if(ptr != NULL) sscanf(buffer, "%d", &numface);

  // Allocate structures
  nbc.eidx.resize(numface); nbc.sp_vtx.resize(numface); nbc.tag.resize(numface);
  con.resize(numface * 3);

  // Get data pointers
  long rb[6];

  for(int i=0, k=0; i<numface; i++) {
    ptr = fgets(buffer, bufsize, file);

    if (ptr != NULL) {
      sscanf(buffer, "%ld %ld %ld %ld %ld %ld", rb, rb+1, rb+2, rb+3, rb+4, rb+5);
      con[k] = rb[0], con[k+1] = rb[1], con[k+2] = rb[2];
      nbc.tag[i] = rb[3], nbc.sp_vtx[i] = rb[4], nbc.eidx[i] = rb[5];
    }

    k += 3;
  }

  fclose(file);
}

void write_nbc(const nbc_data & nbc, const mt_vector<mt_int> & con, int npts, int nelem,
               const std::string filename)
{
  FILE* file = fopen(filename.c_str(), MT_FOPEN_WRITE);
  if(file == NULL) treat_file_open_error(filename, errno);

  int numsurf = nbc.tag.size();
  fprintf(file, "%d # mesh elast %d %d :\n", numsurf, nelem, npts);

  for(int i=0, k=0; i<numsurf; i++) {
    fprintf(file, "%ld %ld %ld %ld %ld %ld\n", long(con[k]), long(con[k+1]),
            long(con[k+2]), long(nbc.tag[i]), long(nbc.sp_vtx[i]), long(nbc.eidx[i]));
    k += 3;
  }

  fclose(file);
}

void write_mesh_selected(mt_meshdata & mesh,
                         std::string format,
                         std::string basename)
{
  if(basename.size() == 0) {
    fprintf(stderr, "Attempted to write mesh with empty filename! Aborting!\n");
    return;
  }

  // if there are no fibers present we create default (1,0,0) fibers
  if(mesh.lon.size() == 0) {
    fprintf(stderr, "%s warning: No fibers detected. Generating default fibers.\n", __func__);
    mesh.lon.resize(3 * mesh.e2n_cnt.size());
    for(size_t eidx=0; eidx < mesh.e2n_cnt.size(); eidx++) {
      mesh.lon[eidx*3+0] = 1;
      mesh.lon[eidx*3+1] = 0;
      mesh.lon[eidx*3+2] = 0;
    }
  }

  if( (format.size() == 0) || (format.compare(carp_txt_fmt) == 0) )
  {
    std::string filename = basename + CARPTXT_ELEM_EXT;
    writeElements(mesh, filename);
    filename = basename + CARPTXT_PTS_EXT;
    writePoints(mesh.xyz, filename);
    filename = basename + CARPTXT_LON_EXT;
    writeFibers(mesh.lon, mesh.e2n_cnt.size(), filename);
  }
  else if(format.compare(carp_bin_fmt) == 0) {
    std::string filename = basename + CARPBIN_ELEM_EXT;
    writeElementsBinary(mesh, filename);
    filename = basename + CARPBIN_PTS_EXT;
    writePointsBinary(mesh.xyz, filename);
    filename = basename + CARPBIN_LON_EXT;
    writeFibersBinary(mesh.lon, mesh.e2n_cnt.size(), filename);
  }
  else if(format.compare(vtk_fmt) == 0) {
    std::string filename = basename + VTK_EXT;
    writeVTKmesh(mesh, filename);
  }
  else if(format.compare(vtk_bin_fmt) == 0) {
    std::string filename = basename + VTK_EXT;
    writeVTKmesh_binary(mesh, filename);
  }
  else if(format.compare(vtu_fmt) == 0) {
    std::string filename = basename + VTU_EXT;
    writeXMLVTKmesh(mesh, filename);
  }
  else if(format.compare(vtk_polydata_fmt) == 0) {
    std::string filename = basename + VTK_EXT;
    writeVTKPolydataMesh(mesh, filename);
  }
  else if(format.compare(mmg_fmt) == 0) {
    std::string filename = basename + MMG_EXT;
    write_mmg_mesh(mesh, filename);
  }
  else if(format.compare(netgen_fmt) == 0) {
    std::string filename = basename + NETGEN_EXT;
    write_netgen_mesh(mesh, filename);
  }
  else if(format.compare(obj_fmt) == 0) {
    std::string filename = basename + OBJ_EXT;
    write_obj_file(mesh, filename);
  }
  else if(format.compare(stellar_fmt) == 0) {
    std::set<mt_int> surf_vtx;
    mt_vector<mt_int> surf_con;
    mt_vector<mt_int> elem_orig;

    MT_MAP<triple<mt_int>, tri_sele> tri_surf;
    MT_MAP<quadruple<mt_int>, quad_sele> quad_surf;

    std::cout << "Computing mesh surface for stellar format .." << std::endl;
    compute_surface(mesh.etype, mesh.e2n_cnt, mesh.e2n_con, tri_surf, quad_surf);
    surfmap_to_vector(tri_surf, quad_surf, surf_con, elem_orig);
    surf_vtx.insert(surf_con.begin(), surf_con.end());

    write_stellar_mesh(mesh, surf_vtx, basename);
  }
  else if(format.compare(vcflow_fmt) == 0) {
    std::string filename = basename + VCFLOW_ELEM_EXT;
    write_vc_flow_elements(mesh, filename);
    filename = basename + VCFLOW_PTS_EXT;
    write_vc_flow_points(mesh.xyz, filename);
    filename = basename + VCFLOW_ADJ_EXT;
    write_vc_flow_adjacency(mesh, filename);
  }
  else if(format.compare(ens_txt_fmt) == 0 || format.compare(ens_bin_fmt) == 0) {
    std::string filename = basename + ENSIGHT_MESH_EXT;
    std::string casefile = basename + ENSIGHT_CASE_EXT;

    bool write_binary = format.compare(ens_bin_fmt) == 0;
    ens_typeinfo ti;
    get_ensight_type_info(mesh, ti);
    write_ensight_mesh(mesh, ti, write_binary, filename);
    write_ensight_case_header(casefile, filename,
                              mesh.lon.size() == mesh.e2n_cnt.size()*6);
  }
  else {
    std::cerr << "Mesh output error: Unrecognized output format. Aborting!" << std::endl;
    return;
  }
}

void read_mesh_selected(mt_meshdata & mesh,
                        std::string format,
                        std::string basename,
                        short readmask)
{
  if(format.size() == 0)
  {
    if( readmask & CRP_READ_ELEM ) readElements_general(mesh, basename);
    if( readmask & CRP_READ_PTS  ) readPoints_general(mesh.xyz, basename);
    if( readmask & CRP_READ_LON  ) readFibers_general(mesh.lon, mesh.e2n_cnt.size(), basename);
  }
  else if(format.compare(carp_txt_fmt) == 0)
  {
    if( readmask & CRP_READ_ELEM ) readElements(mesh, basename + CARPTXT_ELEM_EXT);
    if( readmask & CRP_READ_PTS  ) readPoints(mesh.xyz, basename + CARPTXT_PTS_EXT);
    if( readmask & CRP_READ_LON  ) readFibers(mesh.lon, mesh.e2n_cnt.size(), basename + CARPTXT_LON_EXT);
  }
  else if(format.compare(carp_bin_fmt) == 0) {
    if( readmask & CRP_READ_ELEM ) readElementsBinary(mesh, basename + CARPBIN_ELEM_EXT);
    if( readmask & CRP_READ_PTS  ) readPointsBinary(mesh.xyz, basename + CARPBIN_PTS_EXT);
    if( readmask & CRP_READ_LON  ) readFibersBinary(mesh.lon, mesh.e2n_cnt.size(), basename + CARPBIN_LON_EXT);
  }
  else if( (format.compare(vtk_fmt) == 0) || (format.compare(vtk_bin_fmt) == 0) ) {
    readVTKmesh(mesh, basename + VTK_EXT);
  }
  else if( (format.compare(vtu_fmt) == 0) ) {
    readXMLVTKmesh(mesh, basename + VTU_EXT);
  }
  else if(format.compare(mmg_fmt) == 0) {
    read_mmg_mesh(mesh, basename + MMG_EXT);
  }
  else if(format.compare(netgen_fmt) == 0) {
    read_netgen_mesh(mesh, basename + NETGEN_EXT);
  }
  else if(format.compare(stellar_fmt) == 0) {
    std::set<mt_int> surf_vtx;
    read_stellar_mesh(mesh, surf_vtx, basename);
  }
  else if(format.compare(obj_fmt) == 0) {
    mt_vector<mt_int> nrml_idx;
    mt_vector<mt_real> nrml_xyz;
    read_obj_file(mesh, nrml_idx, nrml_xyz, basename + OBJ_EXT);
  }
  else if(format.compare(purk_fmt) == 0) {
    mt_psdata ps;
    read_purkinje(ps, basename + PURK_EXT);
    psdata_to_mesh(ps, mesh);
  }
  else if(format.compare(vcflow_fmt) == 0) {
    read_vc_flow_elements(mesh, basename + VCFLOW_ELEM_EXT);
    read_vc_flow_points(mesh.xyz, basename + VCFLOW_PTS_EXT);
  }
  else {
    std::cerr << "Mesh input error: Unrecognized input format. Aborting!" << std::endl;
    exit(1);
  }

  // if there are no fibers present we create default (1,0,0) fibers
  if(mesh.lon.size() == 0 && (readmask & CRP_READ_LON) ) {
    fprintf(stderr, "%s warning: No fibers detected. Generating default fibers.\n", __func__);
    mesh.lon.resize(3 * mesh.e2n_cnt.size());
    for(size_t eidx=0; eidx < mesh.e2n_cnt.size(); eidx++) {
      mesh.lon[eidx*3+0] = 1;
      mesh.lon[eidx*3+1] = 0;
      mesh.lon[eidx*3+2] = 0;
    }
  }
}

void read_purkinje(mt_psdata & ps, std::string file)
{
  FILE* psfile = fopen(file.c_str(), MT_FOPEN_READ);
  if(psfile == NULL) treat_file_open_error(file, errno);

  const int bufsize = 2048;
  char buffer[bufsize];
  long unsigned int rcidx, npt;
  double pt[3];
  int    n[2];

  unsigned long int numcab = 0;
  char *ptr = fgets( buffer, bufsize, psfile);
  if(ptr != NULL) sscanf(buffer, "%lu", &numcab);

  ps.cables.resize(numcab);

  for(size_t cidx = 0; cidx < numcab; cidx++)
  {
    mt_pscable & ccab = ps.cables[cidx];

    // read cable index
    ptr = fgets( buffer, bufsize, psfile);
    if(ptr != NULL) sscanf(buffer, "Cable %lu", &rcidx);
    if(rcidx != cidx) {
      fprintf(stderr, "Error: Wrong PS cable order!\n");
      return;
    }
    // read parents
    ptr = fgets( buffer, bufsize, psfile);
    if(ptr != NULL) sscanf(buffer, "%d %d", n, n+1);
    ccab.par1 = n[0], ccab.par2 = n[1];
    // read branches
    ptr = fgets( buffer, bufsize, psfile);
    if(ptr != NULL) sscanf(buffer, "%d %d", n, n+1);
    ccab.br1 = n[0], ccab.br2 = n[1];

    // read number of points
    ptr = fgets( buffer, bufsize, psfile);
    if(ptr != NULL) sscanf(buffer, "%lu", &npt);
    ccab.pts.resize(npt*3);
    // read fiber data
    if(fgets( buffer, bufsize, psfile) != NULL) sscanf(buffer, "%lf", pt+0);
    if(fgets( buffer, bufsize, psfile) != NULL) sscanf(buffer, "%lf", pt+1);
    if(fgets( buffer, bufsize, psfile) != NULL) sscanf(buffer, "%lf", pt+2);
    ccab.size = pt[0], ccab.rgj = pt[1], ccab.cond = pt[2];
    // read points
    for(size_t i=0; i<npt; i++) {
      ptr = fgets( buffer, bufsize, psfile);
      if(ptr != NULL) sscanf(buffer, "%lf %lf %lf", pt, pt+1, pt+2);
      ccab.pts[i*3+0] = pt[0];
      ccab.pts[i*3+1] = pt[1];
      ccab.pts[i*3+2] = pt[2];
    }
  }
  fclose(psfile);
}

void write_purkinje(mt_psdata & ps, std::string file)
{
  FILE* psfile = fopen(file.c_str(), MT_FOPEN_WRITE);
  if(psfile == NULL) treat_file_open_error(file, errno);

  float pt[3];

  unsigned long int numcab = ps.cables.size();
  fprintf(psfile, "%lu\n", numcab);

  for(size_t cidx = 0; cidx < numcab; cidx++)
  {
    mt_pscable & ccab = ps.cables[cidx];
    long unsigned int npt = ccab.pts.size()/3;

    // write cable index
    fprintf(psfile, "Cable %lu\n", (long unsigned int)cidx);
    // write parents
    fprintf(psfile, "%ld %ld\n", long(ccab.par1), long(ccab.par2));
    // write branches
    fprintf(psfile, "%ld %ld\n", long(ccab.br1), long(ccab.br2));
    // write number of points
    fprintf(psfile, "%lu\n", long(npt));

    // write fiber data
    fprintf(psfile, "%f\n", ccab.size);
    fprintf(psfile, "%f\n", ccab.rgj);
    fprintf(psfile, "%f\n", ccab.cond);
    // write points
    for(size_t i=0; i<npt; i++) {
      pt[0] = ccab.pts[i*3+0];
      pt[1] = ccab.pts[i*3+1];
      pt[2] = ccab.pts[i*3+2];
      fprintf(psfile, "%f %f %f\n", pt[0], pt[1], pt[2]);
    }
  }
  fclose(psfile);
}

mt_filename::mt_filename(std::string fname, std::string fmt)
{
  this->assign(fname, fmt);
}

void mt_filename::assign(std::string fname, std::string fmt)
{
  std::string ext_format = "";
  size_t len = 0;

  if(find_extension(fname, VTK_EXT, len))
    ext_format = vtk_bin_fmt;
  else if(find_extension(fname, VTU_EXT, len))
    ext_format = vtu_fmt;
  else if(find_extension(fname, CARPTXT_ELEM_EXT, len))
    ext_format = carp_txt_fmt;
  else if(find_extension(fname, CARPBIN_ELEM_EXT, len))
    ext_format = carp_bin_fmt;
  else if(find_extension(fname, MMG_EXT, len))
    ext_format = mmg_fmt;
  else if(find_extension(fname, OBJ_EXT, len))
    ext_format = obj_fmt;
  else if(find_extension(fname, PURK_EXT, len))
    ext_format = purk_fmt;
  else if(find_extension(fname, ENSIGHT_CASE_EXT, len) ||
          find_extension(fname, ENSIGHT_MESH_EXT, len))
    ext_format = ens_bin_fmt;

  if(ext_format.size())
    fname.resize(len);

  if(fmt.size() > 0)
    format = fmt;
  else
    format = ext_format;

  base   = fname;
  fixBasename(base);
}

bool mt_filename::isSet() const
{
  return base.size() > 0;
}


void remove_in_place(std::string & str, const char c)
{
  size_t widx=0, ridx=0;
  for(ridx=0; ridx < str.size(); ridx++)
    if(str[ridx] != c) str[widx++] = str[ridx];

  str.resize(widx);
}

void remove_indent(std::string & str, const char c)
{
  size_t widx=0, ridx=0;
  while(str[ridx] == c) ridx++;

  for(; ridx < str.size(); ridx++)
    str[widx++] = str[ridx];

  str.resize(widx);
}

void parse_xml_line(std::string line, MT_MAP<std::string, std::string> & fields)
{
  fields.clear();
  remove_indent(line);
  remove_in_place(line, '"');
  remove_in_place(line, '\n');

  mt_vector<std::string> words;
  split_string(line, ' ', words);
  int numwords = words.size();

  if(numwords == 1 && words[0][0] == '<') {
    remove_in_place(words[0], '<');
    remove_in_place(words[0], '>');
    fields["tag"] = words[0];
  }
  else if(numwords > 1) {
    remove_in_place(words[0], '<');
    fields["tag"] = words[0];

    for(int i=1; i<numwords-1; i++) {
      mt_vector<std::string> kvpair;
      split_string(words[i], '=', kvpair);
      if(kvpair.size() == 2)
        fields[kvpair[0]] = kvpair[1];
    }

    if(words[numwords-1].size() != 2)
    {
      remove_in_place(words[numwords-1], '>');
      mt_vector<std::string> kvpair;
      split_string(words[numwords-1], '=', kvpair);
      if(kvpair.size() == 2)
        fields[kvpair[0]] = kvpair[1];
    }
  }
  else
    fields["tag"] = "NONE";
}

void seek_over_char(FILE* fd, const char c)
{
  char buff;
  do {
    fread(&buff, 1, 1, fd);
  } while(buff == c);
  fseek(fd, -1, SEEK_CUR);
}

void setup_data_format(std::string inp_dat, std::string out_dat, short & data_idx,
                       igb_header & inp_header, igb_header & out_header)
{
  // data_idx may be:
  // 0 : scalar dat file
  // 1 : scalar igb file
  // 2 : vec file
  // 3 : vec3 igb file
  // 4 : vec9 igb file

  data_idx = -1;

  if(endswith(inp_dat, DAT_EXT))
    data_idx = 0;
  else if(endswith(inp_dat, VEC_EXT))
    data_idx = 2;
  else if( endswith(inp_dat, IGB_EXT) || endswith(inp_dat, DYNPT_EXT) )
  {
    // setup input igb header
    init_igb_header(inp_dat, inp_header);
    read_igb_header(inp_header);

    if(file_exists(out_dat)) {
      init_igb_header(out_dat, out_header);
      read_igb_header(out_header);
    }
    else {
      // copy most values from input igb header, only few will be set correctly later
      out_header = inp_header;
      out_header.filename = out_dat;
    }

    switch(Num_Components[inp_header.v_type]) {
      case 1: data_idx = 1; break;
      case 3: data_idx = 3; break;
      case 9: data_idx = 4; break;
     default: break;
    }
  }

  if(data_idx == -1) {
    fprintf(stderr, "%s error: Input data type could not be identified! Aborting!\n", __func__);
    exit(1);
  }
}

// defines used in the get_face_comp function
#define COMP_HAVE_NO_SLASHES 0
#define COMP_HAVE_VERT 1
#define COMP_HAVE_TEXT 2
#define COMP_HAVE_NORM 4

// helper function to determine which components of a obj-format face definition are set
int get_face_comp(const char* line) {

  int first_space, second_space;
  int first_slash, second_slash;
  int ret = 0;
  int i=0;

  int line_len = 0;
  while(line[line_len] != 0 && line[line_len] != '\n') line_len++;

  while(line[i] != ' ') i++;
  first_space = i++;

  // we are looking for first slash
  while(line[i] != 0 && line[i] != '\n' && line[i] != '/') i++;

  // it may be that there are no slashes, then we return 0
  if(i < line_len) {
    first_slash = i++;
    while(line[i] != '/') i++;
    second_slash = i++;
    while(line[i] != ' ') i++;
    second_space = i;

    if(first_slash != (first_space + 1))  ret |= COMP_HAVE_VERT;
    if(first_slash != (second_slash - 1)) ret |= COMP_HAVE_TEXT;
    if(second_slash != (second_space - 1)) ret |= COMP_HAVE_NORM;
  }

  return ret;
}

void read_obj_file(mt_meshdata & mesh, mt_vector<mt_int> & nrml_idx,
                   mt_vector<mt_real> & nrml_xyz, std::string file)
{
  int current_material = 0;
  const int linelength = 2048;
  char current_line[linelength];
  int line_number = 0;

  // open scene
  FILE* fd = fopen( file.c_str(), MT_FOPEN_READ);
  if(fd == NULL)
    treat_file_open_error(file, errno);

  //parser loop
  while( fgets(current_line, linelength, fd) )
  {
    line_number++;

    //skip comments
    if(current_line[0] == '#') continue;

    // increment material index for every usemtl
    else if(current_line[0] == 'u' && current_line[1] == 's' && current_line[2] == 'e')
      current_material++;

    //parse objects
    else if(current_line[0] == 'v' && current_line[1] == ' ') //process vertex
    {
      float e1, e2, e3;
      sscanf(current_line, "v %f %f %f", &e1, &e2, &e3);
      mesh.xyz.push_back(e1);
      mesh.xyz.push_back(e2);
      mesh.xyz.push_back(e3);
    }

    else if(current_line[0] == 'v' && current_line[1] == 'n') //process vertex normal
    {
      float e1, e2, e3;
      sscanf(current_line, "vn %f %f %f", &e1, &e2, &e3);
      nrml_xyz.push_back(e1);
      nrml_xyz.push_back(e2);
      nrml_xyz.push_back(e3);
    }

    else if(current_line[0] == 'f') //process face
    {
      int v0 = 0, t0 = 0, n0 = 0, v1 = 0, t1 = 0, n1 = 0, v2 = 0, t2 = 0, n2 = 0;
      int comp_flags = get_face_comp(current_line);

      switch(comp_flags) {
        case COMP_HAVE_NO_SLASHES:
          sscanf(current_line, "f %d %d %d", &v0, &v1, &v2);
          break;
        case COMP_HAVE_VERT:
          sscanf(current_line, "f %d// %d// %d//", &v0, &v1, &v2);
          break;
        case COMP_HAVE_VERT | COMP_HAVE_NORM:
          sscanf(current_line, "f %d//%d %d//%d %d//%d", &v0, &n0, &v1, &n1, &v2, &n2);
          break;
        case COMP_HAVE_VERT | COMP_HAVE_TEXT:
          sscanf(current_line, "f %d/%d/ %d/%d/ %d/%d/", &v0, &t0, &v1, &t1, &v2, &t2);
          break;
        case COMP_HAVE_VERT | COMP_HAVE_NORM | COMP_HAVE_TEXT:
          sscanf(current_line, "f %d/%d/%d %d/%d/%d %d/%d/%d", &v0, &t0, &n0, &v1, &t1, &n1, &v2, &t2, &n2);
          break;
        default:
          fprintf(stderr, "Error: Cannot parse face definition in line %d! Aborting!\n", line_number);
          exit(1);
      }

      mesh.e2n_cnt.push_back(3);
      mesh.etags.push_back(current_material);
      mesh.etype.push_back(Tri);

      mesh.e2n_con.push_back(v0-1);
      mesh.e2n_con.push_back(v1-1);
      mesh.e2n_con.push_back(v2-1);

      nrml_idx.push_back(n0-1);
      nrml_idx.push_back(n1-1);
      nrml_idx.push_back(n2-1);
    }
  }

  fclose(fd);

  generate_default_fib(mesh.e2n_cnt.size(), mesh.lon);
}

void write_obj_file(const mt_meshdata & mesh, std::string filename)
{
  if(mesh.etype[0] != Tri) {
    fprintf(stderr, "Error: Wavefront obj files support only triangle elements! Aborting!\n");
    exit(1);
  }

  FILE* fd = fopen(filename.c_str(), MT_FOPEN_WRITE);
  if(!fd) treat_file_open_error(filename, errno);

  size_t nele = mesh.e2n_cnt.size();
  mt_vector<mt_real> nrml(nele * 3);
  for(size_t eidx = 0; eidx < nele; eidx++)
  {
    mt_int v1 = mesh.e2n_con[eidx*3+0];
    mt_int v2 = mesh.e2n_con[eidx*3+1];
    mt_int v3 = mesh.e2n_con[eidx*3+2];

    vec3r e1 = vec3r(mesh.xyz.data() + v2*3) - vec3r(mesh.xyz.data() + v1*3);
    vec3r e2 = vec3r(mesh.xyz.data() + v3*3) - vec3r(mesh.xyz.data() + v1*3);
    vec3r n  = unit_vector(e1.crossProd(e2));

    nrml[eidx*3+0] = n.x;
    nrml[eidx*3+1] = n.y;
    nrml[eidx*3+2] = n.z;
  }

  fprintf(fd, "# Wavefront OBJ file\n");
  fprintf(fd, "# Generated by meshtool\n");
  fprintf(fd, "\n");

  size_t nnod = mesh.xyz.size() / 3;
  for(size_t i=0; i<nnod; i++)
    fprintf(fd, "v %f %f %f\n",
            float(mesh.xyz[i*3+0]), float(mesh.xyz[i*3+1]), float(mesh.xyz[i*3+2]));
  fprintf(fd, "\n");

  for(size_t i=0; i<nele; i++)
    fprintf(fd, "vn %f %f %f\n",
            float(nrml[i*3+0]), float(nrml[i*3+1]), float(nrml[i*3+2]));
  fprintf(fd, "\n");

  for(size_t i=0; i<nele; i++)
    fprintf(fd, "f %d//%d %d//%d %d//%d\n",
            int(mesh.e2n_con[i*3+0]+1), int(i+1), int(mesh.e2n_con[i*3+1]+1), int(i+1),
            int(mesh.e2n_con[i*3+2]+1), int(i+1));
  fprintf(fd, "\n");
}

bool is_integer(const std::string & s)
{
  int buff;
  char tailbuff[1024];

  int nread = sscanf(s.c_str(), "%d%s", &buff, tailbuff);

  if(nread == 1) return true;
  else           return false;
}

bool is_float(const std::string & s)
{
  float buff;
  char tailbuff[1024];

  int nread = sscanf(s.c_str(), "%f%s", &buff, tailbuff);

  if(nread == 1) return true;
  else           return false;
}

double timediff_sec(struct timeval & t1, struct timeval & t2)
{
  return ((t2.tv_sec - t1.tv_sec)*1e6 + (t2.tv_usec - t1.tv_usec))/1e6;
}

void write_full_mesh_connectivity(const mt_meshdata & mesh, const std::string filename)
{
  std::cout << "Writing full mesh connectivity to: " << filename << std::endl;

  FILE* fd = fopen(filename.c_str(), MT_FOPEN_WRITE);
  if(fd == NULL) treat_file_open_error(filename, errno);

  mesh.n2e_cnt.write(fd);
  mesh.n2e_con.write(fd);
  mesh.n2n_cnt.write(fd);
  mesh.n2n_con.write(fd);

  fclose(fd);
}

void read_full_mesh_connectivity(mt_meshdata & mesh, const std::string filename)
{
  std::cout << "Reading full mesh connectivity from: " << filename << std::endl;

  FILE* fd = fopen(filename.c_str(), MT_FOPEN_READ);
  if(fd == NULL) treat_file_open_error(filename, errno);

  mesh.n2e_cnt.read(fd);
  mesh.n2e_con.read(fd);
  mesh.n2n_cnt.read(fd);
  mesh.n2n_con.read(fd);

  fclose(fd);

  bucket_sort_offset(mesh.e2n_cnt, mesh.e2n_dsp);
  bucket_sort_offset(mesh.n2e_cnt, mesh.n2e_dsp);
  bucket_sort_offset(mesh.n2n_cnt, mesh.n2n_dsp);
}

void write_split_file(mt_vector<split_item> & splitlist, std::string filename)
{
  FILE* fd = fopen(filename.c_str(), MT_FOPEN_WRITE);
  if(!fd) treat_file_open_error(filename, errno);

  fprintf(fd, "%lu\n", (unsigned long int)splitlist.size());

  for(auto it = splitlist.begin(); it != splitlist.end(); ++it)
  {
    fprintf(fd, "%d %d %d\n", it->oldIdx, it->elemIdx, it->newIdx);
  }

  fclose(fd);
}

void read_split_file(mt_vector<split_item> & splitlist, std::string filename)
{
  FILE* fd = fopen(filename.c_str(), MT_FOPEN_READ);
  if(!fd) treat_file_open_error(filename, errno);

  long int size = 0;
  fscanf(fd, "%lu\n", &size);

  splitlist.resize(size);

  for(size_t i = 0; i < splitlist.size(); ++i)
    fscanf(fd, "%d %d %d\n",
            &splitlist[i].oldIdx, &splitlist[i].elemIdx, &splitlist[i].newIdx);

  fclose(fd);
}

