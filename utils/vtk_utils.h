/**
* @file vtk_utils.h
* @brief VTK format specific IO utils.
* @author Aurel Neic
* @version
* @date 2016-12-13
*/


#ifndef _VTK_UTILS
#define _VTK_UTILS

/// often needed parameters
struct vtk_parameters {
  bool binary_data                = false;
  unsigned long int numelems      = 0;
  unsigned long int numpoints     = 0;
  unsigned long int numelemdata   = 0;
  bool points_read                = false;
  bool elems_read                 = false;
  bool elem_types_read            = false;
  bool is_polydata                = false;
};

/// read the next non-empty line
char* next_line(char* buff, const int buffsize, FILE* fin);

/**
* @brief Read the point coordinates
*
* @param mesh     the mesh data struct.
* @param par      the vtk parameters struct.
* @param buff     the string buffer used for keyword parsing.
* @param bufsize  the buffer size.
* @param fin      the input file-stream.
*/
void vtk_process_points(mt_meshdata & mesh,
                        vtk_parameters & par,
                        char* buff, const int bufsize, FILE* fin);

/**
* @brief read the element connectivities
*
* @param mesh     the mesh data struct.
* @param par      the vtk parameters struct.
* @param buff     the string buffer used for keyword parsing.
* @param bufsize  the buffer size.
* @param fin      the input file-stream.
*/
void vtk_process_cells(mt_meshdata & mesh,
                       vtk_parameters & par,
                       char* buff, const int bufsize, FILE* fin);

/**
* @brief read the element types
*
* @param mesh     the mesh data struct.
* @param par      the vtk parameters struct.
* @param buff     the string buffer used for keyword parsing.
* @param bufsize  the buffer size.
* @param fin      the input file-stream.
*/
void vtk_process_cellTypes(mt_meshdata & mesh,
                           vtk_parameters & par,
                           char* buff, const int bufsize, FILE* fin);

/**
* @brief read in element-based data (tags, fibers)
*
* @param mesh     the mesh data struct.
* @param par      the vtk parameters struct.
* @param buff     the string buffer used for keyword parsing.
* @param bufsize  the buffer size.
* @param fin      the input file-stream.
*/
void vtk_process_cellData(mt_meshdata & mesh,
                          vtk_parameters & par,
                          char* buff, const int bufsize, FILE* fin);

/**
* @brief Process the next keyword in the vtk file
*
* @param mesh     the mesh data struct.
* @param par      the vtk parameters struct.
* @param buff     the string buffer used for keyword parsing.
* @param bufsize  the buffer size.
* @param fin      the input file-stream.
* @param progress a progress class.
*/
void vtk_process(mt_meshdata & mesh,
                 vtk_parameters & par,
                 char* buff, const int bufsize, FILE* fin,
                 progress<int> & progress);

/**
* @brief Read a txt or binary vtk file.
*
* @param mesh  the mesh data struct.
* @param file  the file name.
*/
void readVTKmesh(mt_meshdata & mesh, std::string file);

/**
* @brief Read a txt vtk file.
*
* @param mesh  the mesh data struct.
* @param file  the file name.
*/
void writeVTKmesh(const mt_meshdata & mesh, std::string file);

/**
* @brief Read a binary vtk file.
*
* @param mesh  the mesh data struct.
* @param file  the file name.
*/
void writeVTKmesh_binary(const mt_meshdata & mesh, std::string file);

void writeVTKPolydataMesh(const mt_meshdata & mesh, std::string file);


/**
* @brief Write an array in base64-encoding
*
* @param vtk_file    The file descriptor we write to.
* @param data_buf    A buffer for the un-encrypted data.
* @param encode_buf  A buffer for teh encrypted data.
* @param data        The data we want to write.
*/
template<class T>
void write_vtu_data(FILE* vtk_file,
                    mt_vector<unsigned char> & data_buf,
                    mt_vector<unsigned char> & encode_buf,
                    mt_vector<T> & data)
{
  long int data_len = data.size()*sizeof(T);
  //data_buf.resize((data.size()*sizeof(T) + sizeof(long int) + 2) / 3 * 4 + 2);
  data_buf.resize(data.size()*sizeof(T) + sizeof(long int));
  char* start = (char*)data_buf.data();
  memcpy(start, &data_len, sizeof(long int));
  start += sizeof(long int);
  memcpy(start, data.data(), data.size()*sizeof(T));
  encode_base64(data_buf, encode_buf);
  fprintf(vtk_file, "%s\n", encode_buf.data());
}

/**
* @brief 
*
* @tparam V data type
* @param data       pointer to data.
* @param size       number of components in dataset.
* @param dpn        d.o.f. per components
* @param name       name of dataset
* @param is_float   whether data is of floating point type
* @param binary     whether to write binary
* @param out        output fd
*/
template<typename V>
void write_VTU_dataset(const V* data, const size_t size, const int dpn,
                       const char *name, bool is_float, bool binary, FILE* out)
{
  mt_vector<float> fbuff;
  mt_vector<int>   ibuff;

  const size_t datasize = size*dpn;

  if(is_float) {
    fbuff.assign(data, data + datasize);

    if(binary == false) {
      fprintf(out, "        <DataArray type=\"Float32\"");
      fprintf(out, " format=\"ascii\" NumberOfComponents=\"%d\"", dpn);
      fprintf(out, " Name=\"%s\">", name);
      for(size_t i=0; i < datasize; i++) {
        if(!(i%dpn))
          fprintf(out, "\n         ");

        fprintf(out, " %g", fbuff[i]);
      }
    }
    else {
      fprintf(out, "        <DataArray type=\"Float32\"");
      fprintf(out, " format=\"binary\" NumberOfComponents=\"%d\"", dpn);
      fprintf(out, " Name=\"%s\">\n", name);

      mt_vector<unsigned char> data_buf, encode_buf;
      write_vtu_data(out, data_buf, encode_buf, fbuff);
    }
  }
  else {
    ibuff.assign(data, data + datasize);

    if(binary == false) {
      fprintf(out, "        <DataArray type=\"Int32\"");
      fprintf(out, " format=\"ascii\" NumberOfComponents=\"%d\"", dpn);
      fprintf(out, " Name=\"%s\">", name);
      for(size_t i=0; i < datasize; i++) {
        if(!(i%dpn))
          fprintf(out, "\n         ");

        fprintf(out, " %d", ibuff[i]);
      }
    }
    else {
      fprintf(out, "        <DataArray type=\"Int32\"");
      fprintf(out, " format=\"binary\" NumberOfComponents=\"%d\"", dpn);
      fprintf(out, " Name=\"%s\">\n", name);

      mt_vector<unsigned char> data_buf, encode_buf;
      write_vtu_data(out, data_buf, encode_buf, ibuff);
    }
  }
  fprintf(out, "\n        </DataArray>\n");
}

void write_VTU_londata(const mt_vector<mt_real> & lon, const size_t nelem,
                       const int num_axes, const char* prefix, bool binary,
                       FILE* output_file);

void write_VTU_elemtags (const mt_vector<mt_int> & tags, const char *name, bool binary,
                        FILE* output_file);

/**
* @brief Write a mesh in the binary (base64-encoded) xml-vtk format.
*
* @param mesh  The mesh is stored.
* @param file  The file name.
*/
void writeXMLVTKmesh(mt_meshdata & mesh, std::string file);  // deprecated
void write_VTU_mesh(const mt_meshdata & mesh, std::string file, bool binary);
void write_VTU_mesh(const mt_meshdata & mesh, bool binary, bool finish_file, FILE* vtk_file);
void finish_VTU(FILE* output_file);

void parse_xml_tag(const std::string & tag, MT_MAP<std::string, std::string> & fields);
void readXMLVTKmesh(mt_meshdata & mesh, std::string file);

/**
* @brief Convert connectivity ordering from CARP to VTK ordering
*
* @param [in]  mesh     The mesh. Connectivity has CARP ordering.
* @param [out] e2n_con  The connectivity in VTK ordering.
*/
void perm_con_carp_to_vtk(const mt_meshdata & mesh, mt_vector<mt_int> & e2n_con);
/**
* @brief Convert connectivity ordering from VTK to CARP ordering
*
* @param [out]  mesh     The mesh. Connectivity has CARP ordering.
* @param [in]   e2n_con  The connectivity in VTK ordering.
*/
void perm_con_vtk_to_carp(mt_meshdata & mesh, const mt_vector<mt_int> & e2n_con);
#endif

