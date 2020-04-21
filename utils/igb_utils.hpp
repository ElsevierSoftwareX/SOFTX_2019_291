/**
* @file igb_utils.hpp
*
* IGB utils can either be included in a header-only syle,
* or be compiled into an object.
*
* For using IGB utils as a header-only lib, just include igb_utils.hpp where needed.
*
* For building object code, create a .cpp file that includes igb_utls.hpp.
* Whenever you only need the header information of igb_utils.hpp,
* define IGB_EXCLUDE_DEFS before including igb_utls.hpp.
*
* @brief IGB format specific IO utils.
* @author Aurel Neic and Elias Karabelas
* @version
* @date 2017-08-04
*/

#ifndef _IGB_UTILS_H
#define _IGB_UTILS_H

#include <stdio.h>    /* printf */
#include <stdlib.h>
#include <string.h>
#include <strings.h>  // needed for strcasecmp()
#include <assert.h>
#include <errno.h>

#include <iostream>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>
#include <map>

// -------------- Bits de statut pour Header_Read et Header_Write ------ */
#define MOT_CLEF_INV 2
#define GRANDEUR_INV 4

/*
 * Types de trames
 */
#define MIN_TRAME 0
#define CC8    0
#define CC4    1
#define HEX   2
#define HEXEDGES  3
#define HEXBRIDGES  4
#define HEXLINES  5
#define HEX2    6
#define MAX_TRAME 6
#define NTRAMES   7

#define LF  0x0A
#define FF  0x0C
// #define CR  0x0D
#define MY_IGB_HEADERSIZE 1024
#define MLF  0x0A
#define MFF  0x0C
#define MCR  0x0D
#define NALLOC 100

#define     NUL   0
#define     INCONNU 0

#define     IGB_BIG_ENDIAN    666666666
#define     IGB_LITTLE_ENDIAN 777777777
#define     N_SYSTEMES  2

// -------------- Constantes diverses ---------------------------------- */

#define     MAXL      80  //  Longueur maximale d'une ligne d'entete */
#define     N_MAX_ITEMS     30  //  Nombre maximal d'items optionnels */
#define     L_MAX_ITEM      49  //  Longueur maximale pour un item
/* ------------------------ TYPES definition ------------------------------ */
#define     IGB_BYTE      1  /* -- byte ----------------------------------- */
#define     IGB_CHAR      2  /* -- Char ----------------------------------- */
#define     IGB_SHORT     3  /* -- short ---------------------------------- */
#define     IGB_LONG      4  /* -- long ----------------------------------- */
#define     IGB_FLOAT     5  /* -- float ---------------------------------- */
#define     IGB_DOUBLE    6  /* -- Double --------------------------------- */
#define     IGB_COMPLEX   7  /* -- 2 x float (real part, imaginary part) -- */
#define     IGB_D_COMPLEX 8  /* -- 2 x Double (real part, imaginary part) - */
#define     IGB_RGBA      9  /* -- 4 x byte (red, green, blue, alpha) ----- */
#define     IGB_STRUCTURE 10 /* -- Structure ------------------------------ */
#define     IGB_POINTER   11 /* -- void * --------------------------------- */
#define     IGB_LIST      12 /* -- List   --------------------------------- */
#define     IGB_INT       13 /* -- integer -------------------------------- */
#define     IGB_UINT      14 /* -- unsigned integer------------------------ */
#define     IGB_USHORT    15 /* -- unsigned short integer------------------ */
#define     IGB_VEC3_f    16 /* -- 3 X float ------------------------------ */
#define     IGB_VEC3_d    17 /* -- 3 X double ----------------------------- */
#define     IGB_VEC4_f    18 /* -- 4 X float ------------------------------ */
#define     IGB_VEC4_d    19 /* -- 4 X double ----------------------------- */
#define     IGB_VEC9_f    20 /* -- 9 X float  ----------------------------- */
#define     IGB_MIN_TYPE  1
#define     IGB_MAX_TYPE  20

/* define for endedness */
#define IGB_ENDIAN_VAL -1.24e5
#define IGB_LITTLE_END_REP 0,48,242,199

// error codes
#define ERR_EOF_IN_HEADER      1
#define ERR_LINE_TOO_LONG      2
#define ERR_UNPRINTABLE_CHAR   3
#define ERR_IGB_SYNTAX         4
#define ERR_UNDEFINED_X_Y_TYPE 5
#define ERR_SIZE_REDEFINED     6
#define ERR_SIZE_NOT_DEFINED   7
#define WARN_DIM_INCONSISTENT  256


/*
   Definition des types List, bytes, Char, Double, complex d_complex
   et rgba
   */
struct List
{
  long    nitems;
  char    *items;
};

typedef  unsigned char byte;

struct S_Complex
{
  float real, imag;
};

struct D_Complex
{
  double real, imag;
};

#ifndef _COMPLEX_DEFINED
struct complex
{
  float   reel ;
  float   imag ;
};

struct d_complex
{
  double  reel ;
  double  imag ;
};

#define _COMPLEX_DEFINED
#endif
union rgba {
  unsigned  long    l;
  byte    b[4];
};

/* Indice de chaque composante dans le vecteur b[] de l'union rgba */
#define RGBA_ROUGE 3
#define RGBA_VERT  2
#define RGBA_BLEU  1
#define RGBA_ALPHA 0

struct igb_header
{
  std::string filename;

  FILE*  fileptr;
  long   v_x, v_y, v_z, v_t ;                  //!< dimensions --------------------
  long   v_type ;                              //!< type arithmetique -------------
  long   v_taille ;                            //!< taille des pixels (type STRUCTURE)
  long   v_systeme;                           //!< big or little endian
  long   v_num ;                              //!< numero de la tranche ----------
  long   v_bin ;                              //!< nombre de couleurs ------------
  long   v_trame ;                            //!< trame(connectivite) -----------
  long   v_lut ;                               //!< nombre de bytes table couleurs
  long   v_comp ;                              //!< nombre de bytes table compres.
  float  v_epais ;                            //!< epaiseur d'une tranche --------
  float  v_org_x, v_org_y, v_org_z, v_org_t ; //!< coin sup gauche --------
  float  v_inc_x, v_inc_y, v_inc_z, v_inc_t ; //!< distance entre pixels -
  float  v_dim_x, v_dim_y, v_dim_z, v_dim_t ; //!< dimension totale -------
  //float *v_vect_z ;                           //!< coord z de chaque tranche -----
  char   v_unites_x[41], v_unites_y[41], v_unites_z[41], v_unites_t[41] ;
  //!< unites de mesure --------------
  char   v_unites[41] ;                       //!< unites de mesure pour les valeurs des pixels -----------
  float  v_facteur, v_zero ;                  //!< facteur d'echelle et valeur du zero -
  char   v_struct_desc[41] ;                  //!< description de la structure ---
  char   v_aut_name[41] ;                     //!< nom de l'auteur ---------------
  void*  v_transparent;                        //!< transparent value for data
  // boolean flags to indicate if a default value has been overridden
  bool   bool_x, bool_y, bool_z, bool_t;
  bool   bool_type;
  bool   bool_taille;
  bool   bool_num;
  bool   bool_bin;
  bool   bool_trame;
  bool   bool_lut;
  bool   bool_comp;
  bool   bool_epais;
  bool   bool_org_x, bool_org_y, bool_org_z, bool_org_t;
  bool   bool_inc_x, bool_inc_y, bool_inc_z, bool_inc_t;
  bool   bool_dim_x, bool_dim_y, bool_dim_z, bool_dim_t;
  bool   bool_unites_x, bool_unites_y, bool_unites_z, bool_unites_t;
  bool   bool_unites;
  bool   bool_facteur, bool_zero;
  bool   bool_struct_desc;
  bool   bool_aut_name;
  bool   bool_transparent;
  char   transstr[257];
  bool   igb_header_initialized;
};


// this are some forward declarations for the case that we do not import
// igb_header_utils.hpp
void make_igb_consistent(igb_header & igb_head);
long igb_system_endian_code();
long init_igb_header(const std::string filename, igb_header & igb_head);
void set_igb_header_datatype(const std::string datatype, igb_header & igb_head);
void set_igb_header_systeme(const std::string datatype, igb_header & igb_head);
long read_igb_header(igb_header & igb_head);
long write_igb_header(igb_header & igb_head);

extern const char* Header_type[IGB_MAX_TYPE+1];
extern size_t Data_Size[IGB_MAX_TYPE+1];
extern const long Num_Components[IGB_MAX_TYPE+1];
extern const long Header_Systeme_No[2];
extern const char* Header_Systeme[2];
extern const char* allowed_keys[];
extern const long allowed_keys_in_numbers[];


#ifndef IGB_EXCLUDE_DEFS
#include "igb_utils_defs.hpp"
#endif

/**
* @brief Byte swapping func.
*
* @tparam T Data type.
* @param in Input val.
*
* @return Output val with swapped bytes.
*/
template<typename T>
T igb_byte_swap(T in)
{
  T out; // output buffer

  char*  inp  = ( char* ) (& in );
  char*  outp = ( char* ) (& out);
  size_t size = sizeof(T);

  switch(size) {
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

/// return whether the system is big endian
bool igb_system_is_big_endian();

/**
* @brief Set the dimensions (number of vals per axis) of the igb data
*
* @param x    x dimension.
* @param y    y dimension.
* @param z    z dimension.
* @param t    time dimension.
* @param igb  The igb header.
*/
void set_igb_header_dim(int x, int y, int z, int t, igb_header & igb);

/**
* @brief Set the units of the igb data.
*
* @param space_unit  Spatial unit.
* @param time_unit   Temporal unit.
* @param val_unit    Function value unit.
* @param igb         The igb header.
*/
void set_igb_header_units(const std::string space_unit,
                   const std::string time_unit,
                   const std::string val_unit,
                   igb_header & igb);

/**
* @brief Read a block of igb time-slices.
*
* Note that the function assumes that init_igb_header and read_igb_header have already
* been called!
*
* @tparam V Data type of the output vector
* @param data_vec    The vector we read into.
* @param igb         The igb header, already set up.
*
* @return number of read values.
*/
template<class V>
size_t read_igb_slice(std::vector<V> & data_vec, igb_header & igb)
{
  if( !(igb.igb_header_initialized && igb.fileptr) )
  {
    fprintf(stderr, "ERROR in %s: IGB header was not initialized!\n", __func__);
    exit(EXIT_FAILURE);
  }

  bool machine_big_endian = igb_system_is_big_endian();
  bool file_big_endian    = igb.v_systeme == IGB_BIG_ENDIAN;
  bool endian_match       = machine_big_endian == file_big_endian;

  int    num_comp       = Num_Components[igb.v_type];
  size_t byte_per_entry = Data_Size[igb.v_type] / num_comp;

  //Load data
  size_t nvals = igb.v_x * igb.v_y * igb.v_z * num_comp;
  char*  rbuff = new char[nvals*byte_per_entry];

  size_t count = fread(rbuff, byte_per_entry, nvals , igb.fileptr);

  switch(igb.v_type) {
    case IGB_BYTE:
    case IGB_CHAR:
      data_vec.assign(rbuff, rbuff+nvals);
      break;
    case IGB_SHORT:
      data_vec.assign((short*)rbuff, ((short*)rbuff)+nvals);
      break;
    case IGB_INT:
      data_vec.assign((int*)rbuff, ((int*)rbuff)+nvals);
      break;
    case IGB_LONG:
      data_vec.assign((long*)rbuff, ((long*)rbuff)+nvals);
      break;
    case IGB_FLOAT:
    case IGB_VEC3_f:
    case IGB_VEC4_f:
    case IGB_VEC9_f:
      data_vec.assign((float*)rbuff, ((float*)rbuff)+nvals);
      break;
    case IGB_DOUBLE:
    case IGB_VEC3_d:
    case IGB_VEC4_d:
      data_vec.assign((double*)rbuff, ((double*)rbuff)+nvals);
      break;
    default:
      fprintf(stderr, "Error: Data type reading not implemented! Aborting!\n");
      exit(EXIT_FAILURE);
  }

  if(!endian_match)
    for(size_t j=0; j<nvals; j++)
      data_vec[j] = igb_byte_swap(data_vec[j]);

  delete [] rbuff;

  if (count != nvals) {
    fprintf(stderr, "Warning in %s: Number of values read does not match given blocksize \n",__func__);
    fprintf(stderr, "               read:     %ld\n", long(count));
    fprintf(stderr, "               expected: %ld\n", long(nvals));
  }

  return count;
}

/**
* @brief Read a block of igb time-slices.
*
* Note that the function assumes that init_igb_header and read_igb_header have already
* been called!
*
* @tparam V Data type of the output vector
* @param data_vec    The vector we read into.
* @param blocksize   The number of time-slices we want to read.
* @param igb         The igb header, already set up.
*
* @return number of read values.
*/
template<class V>
size_t read_igb_block(std::vector<std::vector<V>> & data_vec,
                   const int blocksize,
                   igb_header & igb)
{
  size_t count = 0;

  //allocate vector
  data_vec.resize(blocksize);

  // Load data
  for(int i=0; i < blocksize; i++)
    count += read_igb_slice(data_vec[i], igb);

  return count;
}


/**
* @brief Read a full igb dataset (all time-slices).
**
* @tparam V Data type of the vector we read into.
* @param data_vec The vector we read into.
* @param igb      The igb header struct we read into.
* @param filename The filename of the igb data.
*
* @return 0 on success, -1 else.
*/
template<class V>
int read_igb_data(std::vector<std::vector<V> > & data_vec,
                  igb_header & igb,
                  const std::string & filename)
{
  init_igb_header(filename, igb);
  read_igb_header(igb);

  int ret = read_igb_block(data_vec, igb.v_t, igb);
  fclose(igb.fileptr);
  igb.fileptr = NULL;

  return(ret);
}


/**
* @brief Write a block of igb time-slices.
*
* Note that the function assumes that init_igb_header and write_igb_header have already
* been called!
*
* @tparam V data type.
* @param data_vec  A vector of time slices.
* @param igb_head  The igb header
*
* @return 0 on success, -1 else.
*/
template<class V>
int write_igb_slice(const std::vector<V> & data_vec, igb_header & igb_head)
{
  if( !(igb_head.igb_header_initialized && igb_head.fileptr) )
  {
    fprintf(stderr, "ERROR in %s: IGB header was not initialized!\n", __func__);
    exit(EXIT_FAILURE);
  }

  size_t num_comp       = (size_t) Num_Components[igb_head.v_type];
  size_t byte_per_entry = (size_t) Data_Size[igb_head.v_type] / num_comp;
  size_t nvals          = (size_t) (igb_head.v_x * igb_head.v_y * igb_head.v_z * num_comp);
  size_t checksum       = (size_t) (data_vec.size());
  size_t count = 0;
  size_t vals_written = 0;

  assert(data_vec.size() == size_t(nvals));
  vals_written = 0;
  switch(igb_head.v_type) {
    case IGB_BYTE:
    case IGB_CHAR:
    {
      std::vector<char> wbuff;
      wbuff.assign(data_vec.begin(), data_vec.end());
      const char* wp = (const char*) wbuff.data();
      vals_written =  fwrite(wp, byte_per_entry, nvals, igb_head.fileptr);
      count += vals_written;
      break;
    }
    case IGB_SHORT:
    {
      std::vector<short> wbuff;
      wbuff.assign(data_vec.begin(), data_vec.end());
      const char* wp = (const char*) wbuff.data();
      vals_written =  fwrite(wp, byte_per_entry, nvals, igb_head.fileptr);
      count += vals_written;
      break;
    }
    case IGB_INT:
    {
      std::vector<int> wbuff;
      wbuff.assign(data_vec.begin(), data_vec.end());
      const char* wp = (const char*) wbuff.data();
      vals_written =  fwrite(wp, byte_per_entry, nvals, igb_head.fileptr);
      count += vals_written;
      break;
    }
    case IGB_LONG:
    {
      std::vector<long int> wbuff;
      wbuff.assign(data_vec.begin(), data_vec.end());
      const char* wp = (const char*) wbuff.data();
      vals_written =  fwrite(wp, byte_per_entry, nvals, igb_head.fileptr);
      count += vals_written;
      break;
    }
    case IGB_FLOAT:
    case IGB_VEC3_f:
    case IGB_VEC4_f:
    case IGB_VEC9_f:
    {
      std::vector<float> wbuff;
      wbuff.assign(data_vec.begin(), data_vec.end());
      const char* wp = (const char*) wbuff.data();
      vals_written =  fwrite(wp, byte_per_entry, nvals, igb_head.fileptr);
      count += vals_written;
      break;
    }
    case IGB_DOUBLE:
    case IGB_VEC3_d:
    case IGB_VEC4_d:
    {
      std::vector<double> wbuff;
      wbuff.assign(data_vec.begin(), data_vec.end());
      const char* wp = (const char*) wbuff.data();
      vals_written =  fwrite(wp, byte_per_entry, nvals, igb_head.fileptr);
      count += vals_written;
      break;
    }
    default:
      fprintf(stderr, "Error: Data type reading not implemented! Aborting!\n");
      exit(EXIT_FAILURE);
  }

  if(count != (size_t) (igb_head.v_x *  num_comp))
  {
    fprintf(stderr, "Warning in %s: Number of values written does not match blocksize\n",__func__);
    fprintf(stderr, "               written:  %zd\n", count);
    fprintf(stderr, "               expected: %zd\n", checksum);
  }

  return count;
}

/**
* @brief Write a block of igb time-slices.
*
* Note that the function assumes that init_igb_header and write_igb_header have already
* been called!
*
* @tparam V data type.
* @param data_vec  A vector of time slices.
* @param igb_head  The igb header
*
* @return number of written values.
*/
template<class V>
int write_igb_block(const std::vector<std::vector<V>> & data_vec,
                    igb_header & igb_head)
{
  size_t count = 0;

  for(size_t i=0; i < data_vec.size(); i++)
    count += write_igb_slice(data_vec[i], igb_head);

  return count;
}


/**
* @brief Write an igb dataset, all time-slices at once.
*
* Note that the function assumes that init_igb_header and write_igb_header have already
* been called!
*
* @tparam V data type.
* @param data_vec A vector of time-slice data.
* @param igb      The igb header, already initialized.
*
* @return 0 on success, -1 else.
*/
template<class V>
int write_igb_data(const std::vector<std::vector<V> > & data_vec,
                   igb_header & igb)
{
  assert(data_vec.size() == (size_t)igb.v_t);
  int ret = write_igb_block(data_vec, igb);
  fclose(igb.fileptr);
  return ret;
}

#endif

