/**
* @file io_utils.h
* @brief General IO utils.
* @author Aurel Neic
* @version
* @date 2016-12-13
*/



#ifndef _IO_UTILS
#define _IO_UTILS

#define BIN_HDR_SIZE 1024

#if __BYTE_ORDER == __LITTLE_ENDIAN
#define htobe(x) byte_swap(x)
#define betoh(x) byte_swap(x)
#define htole(x) (x)
#define letoh(x) (x)
#define MT_ENDIANNESS 0
#else
#define htobe(x) (x)
#define betoh(x) (x)
#define htole(x) byte_swap(x)
#define letoh(x) byte_swap(x)
#define MT_ENDIANNESS 1
#endif

#define CRP_READ_ELEM 1
#define CRP_READ_PTS  2
#define CRP_READ_LON  4

#define MT_FOPEN_READ   "rb"
#define MT_FOPEN_WRITE  "wb"
#define MT_FOPEN_APPEND "ab"

#include "mt_utils_base.h"

static const std::string carp_txt_fmt     = "carp_txt";
static const std::string carp_bin_fmt     = "carp_bin";
static const std::string vtk_fmt          = "vtk";
static const std::string vtk_bin_fmt      = "vtk_bin";
static const std::string vtu_fmt          = "vtu";
static const std::string vtk_polydata_fmt = "vtk_polydata";
static const std::string stellar_fmt      = "stellar";
static const std::string purk_fmt         = "purk";
static const std::string mmg_fmt          = "mmg";
static const std::string netgen_fmt       = "neu";
static const std::string obj_fmt          = "obj";
static const std::string vcflow_fmt       = "vcflow";
static const std::string ens_txt_fmt      = "ensight_txt";
static const std::string ens_bin_fmt      = "ensight_bin";

#define  CARPTXT_ELEM_EXT  ".elem"
#define  CARPTXT_LON_EXT   ".lon"
#define  CARPTXT_PTS_EXT   ".pts"
#define  CARPBIN_ELEM_EXT  ".belem"
#define  CARPBIN_LON_EXT   ".blon"
#define  CARPBIN_PTS_EXT   ".bpts"
#define  VTK_EXT           ".vtk"
#define  VTU_EXT           ".vtu"
#define  MMG_EXT           ".mesh"
#define  NETGEN_EXT        ".neu"
#define  PURK_EXT          ".pkje"
#define  EIDX_EXT          ".eidx"
#define  NOD_EXT           ".nod"
#define  SURF_EXT          ".surf"
#define  NBC_EXT           ".neubc"
#define  VTX_EXT           ".vtx"
#define  DAT_EXT           ".dat"
#define  VEC_EXT           ".vec"
#define  IGB_EXT           ".igb"
#define  DYNPT_EXT         ".dynpt"
#define  OBJ_EXT           ".obj"
#define  VCFLOW_ELEM_EXT   "_connectivity.bin"
#define  VCFLOW_PTS_EXT    "_coordinates.bin"
#define  VCFLOW_ADJ_EXT    "_adjacency.bin"
#define  UVC_PTS_EXT       ".uvc_pts"
#define  ENSIGHT_MESH_EXT  ".geo"
#define  ENSIGHT_CASE_EXT  ".case"
#define  ENSIGHT_DATA_EXT  ".ens"
#define  FULL_CON_EXT      ".fcon"
#define  SPLIT_EXT         ".splt"
#define  TXT_TAGS_EXT      ".tags"
#define  BIN_TAGS_EXT      ".btags"

// input / output formats
static const std::string input_formats = carp_txt_fmt+", "+carp_bin_fmt+", "+vtk_fmt+", "+
        vtk_bin_fmt +", "+mmg_fmt+", "+netgen_fmt+", "+purk_fmt+", "+obj_fmt+", "+stellar_fmt+", "+vcflow_fmt;

static const std::string output_formats =carp_txt_fmt+", "+carp_bin_fmt+", "+vtk_fmt+", "+
        vtk_bin_fmt+", "+vtk_polydata_fmt + ", "+mmg_fmt+", "+netgen_fmt+", "+obj_fmt+", "+
        stellar_fmt+", "+vcflow_fmt+", "+ens_txt_fmt;

// include igb function declarations
#define IGB_EXCLUDE_DEFS
#include "igb_utils.hpp"

/**
* @brief Convert an element type string into an enum
*
* @param eletype element type string
*
* @return element type enum
*/
elem_t getElemTypeID(char *eletype);

template<typename T>
T byte_swap(T in)
{
  T out; // output buffer

  char* inp   = ( char* ) (& in );
  char* outp  = ( char* ) (& out);
  size_t size = sizeof(T);

  switch(size)
  {
    case 4:
      outp[0] = inp[3];
      outp[1] = inp[2];
      outp[2] = inp[1];
      outp[3] = inp[0];
      break;

    case 2:
      outp[0] = inp[1];
      outp[1] = inp[0];
      break;

    case 8:
      outp[0] = inp[7], outp[1] = inp[6];
      outp[2] = inp[5], outp[3] = inp[4];
      outp[4] = inp[3], outp[5] = inp[2];
      outp[6] = inp[1], outp[7] = inp[0];
      break;

    default:
      for(size_t i=0; i<size; i++)
        outp[i] = inp[size-1-i];
  }

  return out;
}

/// give explanation why a file open failed
void treat_file_open_error(std::string file, int errnum, bool do_exit = true);
/// give explanation why a file read failed
void treat_file_read_error(std::string file, FILE* fd);
/// some general error happend that we describe in msg
void treat_file_error(const char* msg, FILE* fd);

/**
* @brief Check if file exists.
*
* @param filename [in]   The file to check
*
* @return True if file exists.
*/
bool file_exists(std::string filename);

/// get the size of a file in bytes
size_t file_size(FILE* fd);
size_t file_size(const char* file);

void generate_default_fib(size_t numelem, mt_vector<mt_real> & lon);


/**
* @brief Write element data in text CARP format.
*
* @param mesh [in]  The mesh.
* @param file [in]  Element file name.
*/
void writeElements(mt_meshdata & mesh, std::string file);


/**
* @brief Write element data in binary CARP format.
*
* @param mesh [in]  The mesh.
* @param file [in]  Element file name.
*/
void writeElementsBinary(mt_meshdata & mesh, const std::string file);

/**
* @brief Read element data from text CARP format.
*
* @param mesh [out] The mesh.
* @param file [in]  Element file name.
*/
void readElements(mt_meshdata & mesh, std::string file);

/**
* @brief Read element data from binary CARP format.
*
* @param mesh [out] The mesh.
* @param file [in]  Element file name.
*/
void readElementsBinary(mt_meshdata & mesh, std::string file);

void readElements_general(mt_meshdata & mesh, std::string basename);

/**
* @brief Read element-tags.
*
* @param etags [out] preallocated vector to read the tag-data into
* @param file  [in]  Element-tags file name.
*/
void readElementTags_general(mt_vector<mt_int> & etags, std::string file);

size_t readNumPoints(std::string file);

void readPoints(mt_vector<mt_real> & xyz, std::string file);
void readPointsBinary(mt_vector<mt_real> & xyz, std::string file);
void readPoints_general(mt_vector<mt_real> & xyz, std::string basename);

/**
* @brief Read the UVC coordinates of a mesh and store them in two vectors (one for LV, one for RV).
*
* @param xyz_lv   uvc coordinates of LV
* @param xyz_rv   uvc coordinates of RV
* @param pos_lv   index of LV uvc coordinates w.r.t. global UVC list
* @param pos_rv   index of RV uvc coordinates w.r.t. global UVC list
* @param idx_to_uvc_points  map between global index and local LV/RV list index
* @param file     The file we read from.
*/
void readUVCPoints(mt_vector<mt_real> & xyz_lv,
                   mt_vector<mt_real> & xyz_rv,
                   mt_vector<mt_int> & pos_lv,
                   mt_vector<mt_int> & pos_rv,
                   mt_vector<mixed_tuple<short, size_t> > & idx_to_uvc_point,
                   std::string filename);

void writePoints(mt_vector<mt_real> & xyz, std::string file);
void writePointsBinary(mt_vector<mt_real> & xyz, std::string file);

int readFibers(mt_vector<mt_real> & lon, size_t numelem, std::string file);
int readFibersBinary(mt_vector<mt_real> & lon, size_t numelem, std::string file);
void readFibers_general(mt_vector<mt_real> & lon, size_t numelem, std::string basename);

void writeFibers(mt_vector<mt_real> & lon, size_t numelem, std::string file);
void writeFibersBinary(mt_vector<mt_real> & lon, size_t numelem, std::string file);

void write_surf(const mt_vector<mt_int> & cnt, const mt_vector<mt_int> & con, const std::string filename);

void read_vtx(mt_vector<mt_int> & nod, std::string filename);
void write_vtx(const mt_vector<mt_int> & nod, const std::string filename,
              bool write_header = true);

void read_nbc(nbc_data & nbc, mt_vector<mt_int> & con, std::string filename);
void write_nbc(const nbc_data & nbc, const mt_vector<mt_int> & con, int npts, int nelem, const std::string filename);


/**
* @brief Write a vector to disk as ascii data.
*
* @tparam T   Data type.
* @param [in] vec  Data vector.
* @param [in] file File name.
*/
template<class T>
void write_vector_ascii(const mt_vector<T> & vec, const std::string file, short dpn = 1)
{
  FILE* out_file = fopen(file.c_str(), MT_FOPEN_WRITE);
  if(out_file == NULL) treat_file_open_error(file, errno);

  char msg[1024];
  sprintf(msg, "Writing %s: ", file.c_str());
  PROGRESS<size_t> prg(vec.size(), msg);

  for (size_t i=0; i<vec.size() / dpn; i++ ) {
    for(short j=0; j<dpn-1; j++)
      fprintf(out_file, "%g ", double(vec[i*dpn + j]) );
    fprintf(out_file, "%g\n", double(vec[i*dpn + (dpn-1)]) );
    prg.next();
  }

  prg.finish();
  fclose(out_file);
}
template<class T>
void write_vector_ascii(const mt_vector<mt_point<T> > & vec, const std::string file)
{
  FILE* out_file = fopen(file.c_str(), MT_FOPEN_WRITE);
  if(out_file == NULL) treat_file_open_error(file, errno);

  char msg[1024];
  sprintf(msg, "Writing %s: ", file.c_str());
  PROGRESS<size_t> prg(vec.size(), msg);

  for (size_t i=0; i<vec.size(); i++ ) {
    fprintf(out_file, "%g %g %g\n", double(vec[i].x), double(vec[i].y), double(vec[i].z));
    prg.next();
  }

  prg.finish();
  fclose(out_file);
}

/**
* @brief Read an ascii data vector
*
* @tparam T     Integer or floating point type
*
* @param [out] vec       Data vector.
* @param [in]  file      File to read from
* @param [in]  readFloat Whether to read floating point or integer values.
*/
template<class T>
void read_vector_ascii(mt_vector<T> & vec, const std::string file, bool readFloat = true)
{
  FILE* in_file = fopen(file.c_str(), MT_FOPEN_READ);
  if(in_file == NULL) treat_file_open_error(file, errno);

  const int bufsize = 2056;
  char buffer[bufsize+1];
  char* ptr;
  double fbuf[3];
  long int ibuf[3];
  short entr_per_line = 0, byte_per_line = 0;
  size_t num_lines = 0;

  size_t size_of_file = file_size(in_file);

  // first get the number of lines and entries per line
  ptr = fgets( buffer, bufsize, in_file);
  if(ptr != NULL) {
    byte_per_line = strlen(ptr);
    entr_per_line = sscanf(buffer, "%lf %lf %lf", fbuf, fbuf+1, fbuf+2);
  }

  // we estimate file size to be able to have a progress bar
  size_t num_lines_estimate = (size_of_file / byte_per_line) * 1.25f;
  char msg[1024];
  sprintf(msg, "Reading %s: ", file.c_str());
  PROGRESS<size_t> prg(num_lines_estimate*2, msg);

  if(entr_per_line == 0)
    treat_file_read_error(file, in_file);

  while(ptr) {
    num_lines++;
    ptr = fgets( buffer, bufsize, in_file);
    prg.next();
  }

  vec.resize(num_lines * entr_per_line);
  fseek(in_file, 0L, SEEK_SET);

  // read file into a list so that we dont need to know its size
  if(readFloat) {
    for(size_t l=0, widx=0; l<num_lines; l++)
    {
      fgets( buffer, bufsize, in_file);
      short nread = sscanf(buffer, "%lf %lf %lf", fbuf, fbuf+1, fbuf+2);
      for(short i=0; i < nread; i++)
        vec[widx++] = T(fbuf[i]);

      prg.next();
    }
  }
  else {
    for(size_t l=0, widx=0; l<num_lines; l++)
    {
      fgets( buffer, bufsize, in_file);
      short nread = sscanf(buffer, "%ld %ld %ld", ibuf, ibuf+1, ibuf+2);
      for(short i=0; i < nread; i++)
        vec[widx++] = T(ibuf[i]);

      prg.next();
    }
  }

  prg.finish();
  fclose(in_file);
}

/**
 * The binary_read procedure reads binary data from a file into a memory block.
 * \param vec Output: Vector we read into.
 * \param filename Input: Name of the input file.
 * \param seek Input: Offset into the data file
 */
template<class T> inline
void binary_read(mt_vector<T> & vec, const std::string filename, const int seek = 0)
{
  FILE* fd = fopen(filename.c_str(), MT_FOPEN_READ);
  if(fd == NULL){
    fprintf(stderr, "%s error: Cannot read file %s.\n", __func__, filename.c_str());
    exit(1);
  }

  if(vec.size() == 0) {
    size_t fsize = file_size(fd);
    size_t vsize = fsize / sizeof(T);
    if(fsize % sizeof(T)) {
      fprintf(stderr, "%s error: File size and requested datatype dont match.\n", __func__);
      exit(1);
    }
    vec.resize(vsize);
  }

  if(seek)
    fseek(fd, seek, SEEK_SET);

  fread(vec.data(), sizeof(T), vec.size(), fd);
  fclose(fd);
}

/**
 * The binary_write procedure writes binary data to a file from a memory block.
 * \param _First Input: Iterator pointing to the first element.
 * \param _Last Input: Iterator pointing past the last element.
 * \param _Filename Input: Name of the output file.
 */
template<class _InIt>
inline void binary_write(_InIt _First, _InIt _Last, const std::string _Filename)
{
  std::ofstream os(_Filename.c_str(), std::ios::binary);
  if(! os.good()){
    std::cerr << "ERROR: Cannot write to file " << _Filename << std::endl;
    return;
  }
  os.write((char*)&*_First, std::streamsize(_Last - _First) * sizeof(*_First));
  os.close();
}


/**
* @brief Write a mesh in a given format.
*
* If the given format is an empty string, the function defaults to carp_txt.
* The output format is compared against the *_fmt string constants from
* mesh_mode_funcs.h.
*
* @param mesh     [in]  The mesh
* @param format   [in]  The output format.
* @param basename [in]  Basename of the mesh to output.
*/
void write_mesh_selected(mt_meshdata & mesh,
                         std::string format,
                         std::string basename);
/**
* @brief Read a mesh in a given format.
*
* If the given format is an empty string, the function defaults to carp_txt.
* The output format is compared against the *_fmt string constants from
* mesh_mode_funcs.h.
*
* @param mesh     [out] The mesh
* @param format   [in]  The output format.
* @param basename [in]  Basename of the mesh to output.
* @param readmask [in]  Bitmask for what to read (only affect CARP formats).
*/
void read_mesh_selected(mt_meshdata & mesh,
                        std::string format,
                        std::string basename,
                        short readmask = (CRP_READ_ELEM | CRP_READ_PTS | CRP_READ_LON) );

/**
* @brief Read the data of a purkinje file.
*
* @param [out] ps   PS data struct
* @param [in]  file PS data file name.
*/
void read_purkinje(mt_psdata & ps, std::string file);
/**
* @brief Write the data of a purkinje file.
*
* @param [in]  ps   PS data struct
* @param [in]  file PS data file name.
*
* @post  PS has been writen to disk.
*/
void write_purkinje(mt_psdata & ps, std::string file);

/**
* @brief Class for filename and format management
*/
class mt_filename {
  public:
  std::string base;   ///< the filename without extension
  std::string format; ///< the file format

  /// construct empty filename
  mt_filename(){}
  /// construct file, calls assign()
  mt_filename(std::string fname, std::string fmt);

  /// construct file, format may be empty
  void assign(std::string fname, std::string fmt);
  /// check whether filename is set
  bool isSet() const;
};

/// chars used in the base64 encoding
static const std::string base64_chars =
             "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
             "abcdefghijklmnopqrstuvwxyz"
             "0123456789+/";


/// test if a char is base64 encoded
static inline bool is_base64(unsigned char c) {
  return (isalnum(c) || (c == '+') || (c == '/'));
}

/**
* @brief Encode a data vector to base64 encoding.
*
* @tparam T   Data type of the input data
* @param in   Input data vector
* @param out  Output data vector
*/
template<class T>
void encode_base64(const mt_vector<T> & in, mt_vector<unsigned char> & out)
{
  const unsigned char* input  = (const unsigned char*) in.data();
  size_t               in_len = in.size() * sizeof(T);

  int i = 0;
  int j = 0;
  unsigned char buff3[3], buff4[4];

  size_t widx = 0;
  out.resize((in.size()*sizeof(T) + 2) / 3 * 4 + 2);

  while (in_len--)
  {
    buff3[i++] = *(input++);
    if (i == 3) {
      buff4[0] = (buff3[0] & 0xfc) >> 2;
      buff4[1] = ((buff3[0] & 0x03) << 4) + ((buff3[1] & 0xf0) >> 4);
      buff4[2] = ((buff3[1] & 0x0f) << 2) + ((buff3[2] & 0xc0) >> 6);
      buff4[3] = buff3[2] & 0x3f;

      for(i = 0; (i <4) ; i++)
        out[widx++] = base64_chars[buff4[i]];
      i = 0;
    }
  }

  if(i) {
    for(j = i; j < 3; j++) buff3[j] = '\0';

    buff4[0] = (buff3[0] & 0xfc) >> 2;
    buff4[1] = ((buff3[0] & 0x03) << 4) + ((buff3[1] & 0xf0) >> 4);
    buff4[2] = ((buff3[1] & 0x0f) << 2) + ((buff3[2] & 0xc0) >> 6);
    buff4[3] = buff3[2] & 0x3f;

    for (j = 0; (j < i + 1); j++)
      out[widx++] = base64_chars[buff4[j]];

    while((i++ < 3)) out[widx++] = '=';
  }

  out[widx++] = '\0';
  out.resize(widx);
}

/**
* @brief Decode a data vector from base64 encoding.
*
* @tparam T   Data type of the output data
* @param in   Input data vector
* @param out  Output data vector
*/
template<class T>
void decode_base64(const mt_vector<unsigned char> & in, mt_vector<T> & out)
{
  size_t in_len = in.size();
  int i = 0;
  int j = 0;
  size_t ridx = 0, widx = 0;
  unsigned char buff3[3], buff4[4];

  out.resize( (in.size() + sizeof(T) - 1) / sizeof(T) );
  unsigned char * op = (unsigned char*) out.data();

  while (in_len-- && ( in[ridx] != '=') && is_base64(in[ridx])) {
    buff4[i++] = in[ridx++];

    if (i ==4) {
      for (i = 0; i <4; i++) buff4[i] = base64_chars.find(buff4[i]);

      buff3[0] = (buff4[0] << 2) + ((buff4[1] & 0x30) >> 4);
      buff3[1] = ((buff4[1] & 0xf) << 4) + ((buff4[2] & 0x3c) >> 2);
      buff3[2] = ((buff4[2] & 0x3) << 6) + buff4[3];

      for (i = 0; (i < 3); i++) op[widx++] = buff3[i];
      i = 0;
    }
  }

  if(i) {
    for (j = i; j <4; j++) buff4[j] = 0;
    for (j = 0; j <4; j++) buff4[j] = base64_chars.find(buff4[j]);

    buff3[0] = (buff4[0] << 2) + ((buff4[1] & 0x30) >> 4);
    buff3[1] = ((buff4[1] & 0xf) << 4) + ((buff4[2] & 0x3c) >> 2);
    buff3[2] = ((buff4[2] & 0x3) << 6) + buff4[3];

    for (j = 0; (j < i - 1); j++) op[widx++] = buff3[j];
  }
  out.resize(widx / sizeof(T));
}

/**
* @brief In-place character removal.
*/
void remove_in_place(std::string & str, const char c);
/**
* @brief In-place indentation removal.
*/
void remove_indent(std::string & str, const char c = ' ');

void seek_over_char(FILE* fd, const char c);
/**
* @brief Parse the key value pairs of an xml expression.
*
* The name of the expression is stored as "tag"
*
* @param line   A line with an xml expression.
* @param fields A dictionary with the key value pairs.
*/
void parse_xml_line(std::string line, MT_MAP<std::string, std::string> & fields);

/**
* @brief Find out the type of input data (.dat, .vec, .igb).
*
* @param [in]  inp_dat      Name of input data file.
* @param [in]  out_dat      Name of output data file.
* @param [out] data_idx     Data type index.
* @param [out] inp_header   Input igb header (if data is igb).
* @param [out] out_header   Output igb header (if data is igb).
*/
void setup_data_format(std::string inp_dat, std::string out_dat, short & data_idx,
                       igb_header & inp_header, igb_header & out_header);


/**
* @brief Reader function for the Wavefront obj format.
*
* This format supports only triangle (=surface) meshes.
*
* @param mesh       The mesh we read into.
* @param nrml_idx   Vector holding the surface normal index for each node of an element
* @param nrml_xyz   Vector holding the surface normals indexed by nrml_idx
* @param file       The file name we read from
*/
void read_obj_file(mt_meshdata & mesh, mt_vector<mt_int> & nrml_idx,
                   mt_vector<mt_real> & nrml_xyz, std::string file);


void write_obj_file(const mt_meshdata & mesh, std::string filename);


/**
* @brief Compute an evend distribution of the parts of a value among the entries of a
*        vector.
*
* @tparam T   Integer type.
* @param vec  The vector we distribute in.
* @param val  The value we distribute.
*/
template<class T> inline
void even_distribution(mt_vector<T> & vec, T val)
{
  for(size_t i=0; i<vec.size(); i++) vec[i] = val / vec.size();
  for(size_t i=0; i<val % vec.size(); i++) vec[i]++;
}

/**
* @brief Check wheter a string can be interpreted as an integer.
*
* @param s  The string.
*
* @return   Whether the string can be interpreted as an integer.
*/
bool is_integer(const std::string & s);

/**
* @brief Check wheter a string can be interpreted as an float.
*
* @param s  The string.
*
* @return   Whether the string can be interpreted as an float.
*/
bool is_float(const std::string & s);


/**
* @brief Return the timedifference t2 - t1 in seconds.
*
* @param t1 Time one.
* @param t2 Time two.
*
* @return The time difference.
*/
double timediff_sec(struct timeval & t1, struct timeval & t2);

/**
* @brief Write the compute_full_mesh_connectivity() information to disk. 
*
* @param mesh       The mesh data.
* @param filename   filename used for outputting.
*/
void write_full_mesh_connectivity(const mt_meshdata & mesh, const std::string filename);

/**
* @brief Read the compute_full_mesh_connectivity() information from disk. 
*
* @param mesh       The mesh data.
* @param filename   filename used for reading.
*/
void read_full_mesh_connectivity(mt_meshdata & mesh, const std::string filename);


struct split_item {
  int elemIdx;
  int oldIdx;
  int newIdx;
};

void write_split_file(mt_vector<split_item> & splitlist, std::string filename);
void read_split_file(mt_vector<split_item> & splitlist, std::string filename);

#endif

