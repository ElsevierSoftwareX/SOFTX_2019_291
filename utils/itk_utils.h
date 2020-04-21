/**
* @file itk_utils.h
* @brief A collection of ITK utility functions.
* @authors Matthias Gsell, Aurel Neic
* @version
* @date 2017-02-16
*/


#ifndef _ITK_UTILS_H
#define _ITK_UTILS_H

#include <limits.h>

#define ITKSMOOTH_CHECK_RAD 1
#define ITK_MAX_PIXEL_COMP 4

enum ITK_AXES
{
  ITK_AXIS_0 = 0,
  ITK_AXIS_X = 1,
  ITK_AXIS_Y = 2,
  ITK_AXIS_Z = 4
};

/**
* @brief Convert the ITK type string to a type ID
*
* @param typestr The type string.
*
* @return The type ID.
*/
int itk_get_datatype_id(const char * typestr)
{
  int id = -1;
  if (strncmp(typestr, "bit", 3) == 0)                  id = 0;
  else if (strncmp(typestr, "unsigned_char", 13) == 0)  id = 1;
  else if (strncmp(typestr, "char", 4) == 0)            id = 2;
  else if (strncmp(typestr, "unsigned_short", 14) == 0) id = 3;
  else if (strncmp(typestr, "short", 5) == 0)           id = 4;
  else if (strncmp(typestr, "unsigned_int", 12) == 0)   id = 5;
  else if (strncmp(typestr, "int", 3) == 0)             id = 6;
  else if (strncmp(typestr, "unsigned_long", 13) == 0)  id = 7;
  else if (strncmp(typestr, "long", 4) == 0)            id = 8;
  else if (strncmp(typestr, "float", 5) == 0)           id = 9;
  else if (strncmp(typestr, "double", 6) == 0)          id = 10;

  return id;
}

/**
* @brief Convert the type ID into an ITK type string
*
*/
const char * itk_get_datatype_str(const int type)
{
  switch (type)
  {
    case 0: return "bit";
    case 1: return "unsigned_char";
    case 2: return "char";
    case 3: return "unsigned_short";
    case 4: return "short";
    case 5: return "unsigned_int";
    case 6: return "int";
    case 7: return "unsigned_long";
    case 8: return "long";
    case 9: return "float";
    case 10: return "double";
    default: return NULL;
  }
}

/// Get the size in bytes from an type ID.
unsigned int itk_get_datatype_size(const int type)
{
  switch (type)
  {
    case 0: return UINT_MAX;
    case 1: return sizeof(unsigned char);
    case 2: return sizeof(char);
    case 3: return sizeof(unsigned short);
    case 4: return sizeof(short);
    case 5: return sizeof(unsigned int);
    case 6: return sizeof(int);
    case 7: return sizeof(unsigned long);
    case 8: return sizeof(long);
    case 9: return sizeof(float);
    case 10: return sizeof(double);
    default: return UINT_MAX;
  }
}

/**
* @brief Read one ascii pixel.
*
* @param fp          The file descriptor to read from.
* @param ncomp       The number of components of one pixel.
* @param comptypeid  The type ID of one pixel component.
* @param data        The data pointer where the pixel is read.
*
* @return the number of read components.
*/
unsigned int itk_read_pixel_ascii(FILE * fp, const unsigned int ncomp, const int comptypeid, char * data)
{
  const unsigned int compsz = itk_get_datatype_size(comptypeid);
  if (compsz == UINT_MAX) return 0;
  unsigned int i;
  bool readerror = false;
  char compdata[8];
  const char * fmtstr = NULL;
  if ((comptypeid == 1) || (comptypeid == 3) || (comptypeid == 5) || (comptypeid == 7)) fmtstr = "%lu";
  else if ((comptypeid == 2) || (comptypeid == 4) || (comptypeid == 6) || (comptypeid == 8)) fmtstr = "%ld";
  else if ((comptypeid == 9) || (comptypeid == 10)) fmtstr = "%lf";
  for (i = 0; (i < ncomp) && (!readerror); i++)
  {
    if (fscanf(fp, fmtstr, compdata) == 1)
    {
      memcpy(&data[i*compsz], &compdata, compsz*sizeof(char));
      //TODO perform byte swap !!!
    }
    else readerror = true;
  }
  return i;
}



/**
* @brief The ITK image class
*/
class itk_image
{
  public:

  unsigned int npix;    ///< The number of pixels.
  unsigned int ncomp;   ///< The number of components.
  unsigned int compsz;  ///< The size in bytes of one component.
  unsigned int pixsz;   ///< The size in bytes of one pixel.
  unsigned int datasz;  ///< The overall size in bytes of the data.

  triple<unsigned int> dim;   ///< The number of pixel in each dimension.
  int comptypeid;             ///< ID encoding the data type of the data components.
  triple<double> pixspc;      ///< The spacing of the pixels.
  triple<double> orig;        ///< The origin coordinate.
  mt_vector<char> databuff;   ///< The data buffer.

  private:

  /**
  * @brief Since binary ITK data is stored in big endian, byteswap_data() does a
  * byte swap of the whole data array if the machine is little endian.
  */
  void byteswap_data()
  {
    #if __BYTE_ORDER == __LITTLE_ENDIAN
    if(itk_get_datatype_size(this->comptypeid) == 1)
      return;

    switch(this->compsz)
    {
      default:
      case 2:
      {
        short* dp = (short*) this->databuff.data();
        for(unsigned int i=0, idx=0; i<this->npix; i++)
          for(unsigned int j=0; j<this->ncomp; j++, idx++)
            dp[idx] = byte_swap(dp[idx]);
        break;
      }
      case 4:
      {
        int* dp = (int*) this->databuff.data();
        for(unsigned int i=0, idx=0; i<this->npix; i++)
          for(unsigned int j=0; j<this->ncomp; j++, idx++)
            dp[idx] = byte_swap(dp[idx]);
        break;
      }
      case 8:
      {
        double* dp = (double*) this->databuff.data();
        for(unsigned int i=0, idx=0; i<this->npix; i++)
          for(unsigned int j=0; j<this->ncomp; j++, idx++)
            dp[idx] = byte_swap(dp[idx]);
        break;
      }
    }
    #endif
  }

  /**
  * @brief Convert the type of the internal data buffer.
  *
  * The data types are defined through the input pointer types.
  *
  * @tparam T   Output type.
  * @tparam S   Input type.
  * @param out  Output pointer.
  * @param in   Input pointer.
  */
  template<typename T, typename S>
  void convert_dtype(T* out, S* in)
  {
    for(unsigned int c = 0; c < this->npix * this->ncomp; c++)
      out[c] = in[c];
  }

  public:
  /// Empty constructor initializes to 0
  itk_image() :
  npix(0), ncomp(0), compsz(0), pixsz(0), datasz(0),
  dim({0, 0, 0}), comptypeid(0), pixspc({0.0, 0.0, 0.0}),
  orig({0.0, 0.0, 0.0}), databuff()
  {}

  /**
  * @brief Read an itk image from file
  *
  * @param filename The filename.
  */
  void read_file(const char* filename)
  {
    FILE * fp;
    fp = fopen(filename, MT_FOPEN_READ);
    if (fp) {
      int datatype = 0; // 0 ... not set, 1 ... ascii, 2 ... binary

      char line[256], tmpstr[32];
      int readflag = 0;
      while (fgets(line, 256*sizeof(char), fp) && (!readflag))
      {
        if (strncmp(line, "ASCII", 5) == 0)
          datatype = 1;
        else if (strncmp(line, "BINARY", 6) == 0)
          datatype = 2;
        else if (strncmp(line, "DATASET", 7) == 0)
        {
          sscanf(line, "%*s %s", tmpstr);
          if (strncmp(tmpstr, "STRUCTURED_POINTS", 17) != 0) readflag = -1;
        }
        else if (strncmp(line, "DIMENSIONS", 10) == 0)
        {
          sscanf(line, "%*s %u %u %u", &this->dim.v1, &this->dim.v2, &this->dim.v3);
          printf("[ITKIMAGE]: dimensions %u x %u x %u ( %u )\n", this->dim.v1, this->dim.v2, this->dim.v3, this->dim.v1*this->dim.v2*this->dim.v3);
        }
        else if (strncmp(line, "SPACING", 7) == 0)
        {
          sscanf(line, "%*s %lf %lf %lf", &this->pixspc.v1, &this->pixspc.v2, &this->pixspc.v3);
          printf("[ITKIMAGE]: spacing %lf x %lf x %lf\n", this->pixspc.v1, this->pixspc.v2, this->pixspc.v3);
        }
        else if (strncmp(line, "ORIGIN", 6) == 0)
        {
          sscanf(line, "%*s %lf %lf %lf", &this->orig.v1, &this->orig.v2, &this->orig.v3);
          printf("[ITKIMAGE]: origin %lf x %lf x %lf\n", this->orig.v1, this->orig.v2, this->orig.v3);
        }
        else if (strncmp(line, "POINT_DATA", 10) == 0)
        {
          sscanf(line, "%*s %u", &this->npix);
          printf("[ITKIMAGE]: pointdata %u\n", this->npix);
          char * retval = fgets(line, 256*sizeof(char), fp);
          if ((retval) && (strncmp(line, "SCALARS", 7) == 0))
          {
            char comptypestr[16];
            this->ncomp = 1;
            sscanf(line, "%*s %s %s %u", tmpstr, comptypestr, &this->ncomp);
            printf("[ITKIMAGE]: pointdata, name '%s', type '%s', length %u\n", tmpstr, comptypestr, this->ncomp);
            this->comptypeid = itk_get_datatype_id(comptypestr);
            this->compsz = itk_get_datatype_size(this->comptypeid);
            if (this->compsz != UINT_MAX)
            {
              this->pixsz = this->compsz*this->ncomp;
              this->datasz = this->npix*this->pixsz;
              printf("[ITKIMAGE]: pointdata, size %u bytes\n", this->datasz);
              long int fpos = ftell(fp);
              if (fgets(line, 256*sizeof(char), fp) && (strncmp(line, "LOOKUP_TABLE", 12) == 0))
                printf("[ITKIMAGE]: LOOKUP_TABLE found\n");
              else
              {
                printf("[ITKIMAGE]: no LOOKUP_TABLE found\n");
                fseek (fp, fpos, SEEK_SET);
              }
              if (datatype == 1)
              {
                this->databuff.resize(this->datasz);
                unsigned int i = 0;
                bool readerror = false;

                // ascii read data
                for (i = 0; (i < this->npix) && (!readerror); i++)
                  if (itk_read_pixel_ascii(fp, this->ncomp, this->comptypeid, &this->databuff[i * this->pixsz]) != this->ncomp)
                    readerror = true;

                if (readerror) readflag = -5;
              }
              else if (datatype == 2)
              {
                // allocate data
                this->databuff.resize(this->datasz);

                // binary read data
                if (fread(this->databuff.data(), sizeof(char), this->datasz, fp) != this->datasz)
                  readflag = -5;
                else
                  byteswap_data();
              }
              else readflag = -4;
            }
            else readflag = -3;
          }
          else if ((retval) && (strncmp(line, "COLOR_SCALARS", 13) == 0))
          {
            this->ncomp = 1;
            sscanf(line, "%*s %s %u", tmpstr, &this->ncomp);
            printf("[ITKIMAGE]: pointdata, color_scalars, name '%s', length %u\n", tmpstr, this->ncomp);
            this->comptypeid = 1; // unsigned char
            this->compsz = 1;
            if (this->compsz != UINT_MAX)
            {
              this->pixsz = this->compsz*this->ncomp;
              this->datasz = this->npix*this->pixsz;
              printf("[ITKIMAGE]: pointdata, size %u bytes\n", this->datasz);
              if (datatype == 1)
              {
                this->databuff.resize(this->datasz);
                unsigned int i = 0;
                float compval = 0.0;
                unsigned char disccompval = 0;
                bool readerror = false;
                for (i = 0; (i < this->npix*this->ncomp) && (!readerror); i++)
                {
                  if (fscanf(fp, "%f", &compval) == 1)
                  {
                    disccompval = (unsigned char) (compval*255.0);
                    this->databuff[i] = (char) disccompval;
                  }
                  else readerror = true;
                }
                if (readerror) readflag = -5;
              }
              else if (datatype == 2)
              {
                this->databuff.resize(this->datasz);
                if (fread(this->databuff.data(), sizeof(char), this->datasz, fp) != this->datasz)
                  readflag = -5;
                else
                  byteswap_data();
              }
              else readflag = -4;
            }
            else readflag = -3;
          }
          else readflag = -2;
        }
      }
      fclose(fp);

      if (readflag < 0)
      {
        switch (readflag)
        {
          case -1 : printf("[ITKIMAGE]: error while reading file, not supported data set !\n"); break;
          case -2 : printf("[ITKIMAGE]: error while reading file, no scalar point data found !\n"); break;
          case -3 : printf("[ITKIMAGE]: error while reading file, invalid pixel data type id %d !\n", this->comptypeid); break;
          case -4 : printf("[ITKIMAGE]: error while reading file, ascii format not supported yet !\n"); break;
          case -5 : printf("[ITKIMAGE]: error while reading file, pixel data size error !\n"); break;
        }
      }
    }
    else {
      printf("[ITKIMAGE]: failed to open '%s' \n", filename);
    }
  }


  /**
  * @brief Write an itk image file.
  *
  * @param filename  The filename.
  */
  void write_file(const char * filename)
  {
    byteswap_data();

    FILE * fp;
    fp = fopen(filename, "wb");
    if (fp)
    {
      char line[256];
      sprintf(line, "# vtk DataFile Version 3.0\n"
                    "GENERATED BY CFEL\n"
                    "BINARY\n"
                    "DATASET STRUCTURED_POINTS\n"
                    "DIMENSIONS %u %u %u\n", this->dim.v1, this->dim.v2, this->dim.v3);
      fwrite(line, sizeof(char), strlen(line), fp);
      sprintf(line, "SPACING %.16e %.16e %.16e\n"
                    "ORIGIN %.16e %.16e %.16e\n", this->pixspc.v1, this->pixspc.v2, this->pixspc.v3,
                                                  this->orig.v1, this->orig.v2, this->orig.v3);
      fwrite(line, sizeof(char), strlen(line), fp);
      sprintf(line, "POINT_DATA %u\n"
                    "SCALARS scalars %s %u\n"
                    "LOOKUP_TABLE default\n", this->npix, itk_get_datatype_str(this->comptypeid), this->ncomp);
      fwrite(line, sizeof(char), strlen(line), fp);
      fwrite(this->databuff.data(), sizeof(char), this->datasz, fp);
      fclose(fp);
      printf("[ITKIMAGE]: written to '%s' \n", filename);
    }
    else {
      printf("[ITKIMAGE]: failed to open '%s' \n", filename);
    }

    byteswap_data();
  }

  void write_file_color_scalars(const char * filename)
  {
    if (this->comptypeid != 1) {
      printf("[ITKIMAGE]: failed to write color-scalar-file, wrong data type %s\n", itk_get_datatype_str(this->comptypeid));
    }
    else {
      byteswap_data();
  
      FILE * fp;
      fp = fopen(filename, "wb");
      if (fp)
      {
        char line[256];
        sprintf(line, "# vtk DataFile Version 3.0\n"
                      "GENERATED BY CFEL\n"
                      "BINARY\n"
                      "DATASET STRUCTURED_POINTS\n"
                      "DIMENSIONS %u %u %u\n", this->dim.v1, this->dim.v2, this->dim.v3);
        fwrite(line, sizeof(char), strlen(line), fp);
        sprintf(line, "SPACING %.16e %.16e %.16e\n"
                      "ORIGIN %.16e %.16e %.16e\n", this->pixspc.v1, this->pixspc.v2, this->pixspc.v3,
                                                    this->orig.v1, this->orig.v2, this->orig.v3);
        fwrite(line, sizeof(char), strlen(line), fp);
        sprintf(line, "POINT_DATA %u\n"
                      "COLOR_SCALARS scalars %u\n", this->npix, this->ncomp);
        fwrite(line, sizeof(char), strlen(line), fp);
        fwrite(this->databuff.data(), sizeof(char), this->datasz, fp);
        fclose(fp);
        printf("[ITKIMAGE]: written to '%s' \n", filename);
      }
      else {
        printf("[ITKIMAGE]: failed to open '%s' \n", filename);
      }
  
      byteswap_data();  
    }
  }

  void convert(const int dtype)
  {
    if ((dtype < 1) || (dtype > 11))
      return;

    if (dtype == this->comptypeid)
      return;

    itk_image tmp_img;
    tmp_img.assign(*this, itk_get_datatype_str(dtype));

    this->comptypeid = dtype;
    this->compsz     = tmp_img.compsz;
    this->pixsz      = tmp_img.pixsz;
    this->datasz     = tmp_img.datasz;
    this->databuff   = tmp_img.databuff;
  }

  /**
  * @brief Copy the contents of a given itk_image, possibly changing the datatype
  *
  * @param inp   Input itk_image.
  * @param type  The new data type. Use "" to keep the old type.
  *
  * @post Data "inp" has been copied.
  */
  void assign(const itk_image & inp, const char* type = "")
  {
    this->npix   = inp.npix;
    this->ncomp  = inp.ncomp;
    this->dim    = inp.dim;
    this->pixspc = inp.pixspc;
    this->orig   = inp.orig;

    if(strlen(type) == 0) {
      this->compsz     = inp.compsz;
      this->comptypeid = inp.comptypeid;
    }
    else {
      this->comptypeid = itk_get_datatype_id(type);
      this->compsz     = itk_get_datatype_size(this->comptypeid);
    }

    this->pixsz  = this->ncomp * this->compsz;
    this->datasz = this->npix  * this->pixsz;

    this->databuff.assign(this->datasz, char(0));

    const char *input  = inp.databuff.data();
    char       *output = this->databuff.data();

    switch(this->comptypeid) // dtype to convert to
    {
      case 1:  // (unsigned char);
      {
         switch(inp.comptypeid) {
          case 1:  convert_dtype((unsigned char*)output, (unsigned char*) input); break;
          case 2:  convert_dtype((unsigned char*)output, (char*) input); break;
          case 3:  convert_dtype((unsigned char*)output, (unsigned short*) input); break;
          case 4:  convert_dtype((unsigned char*)output, (short*) input); break;
          case 5:  convert_dtype((unsigned char*)output, (unsigned int*) input); break;
          case 6:  convert_dtype((unsigned char*)output, (int*) input); break;
          case 7:  convert_dtype((unsigned char*)output, (unsigned long*) input); break;
          case 8:  convert_dtype((unsigned char*)output, (long*) input); break;
          case 9:  convert_dtype((unsigned char*)output, (float*) input); break;
          case 10: convert_dtype((unsigned char*)output, (double*) input); break;
          default: assert(0);
        }
        break;
      }
      case 2:  // (char);
      {
         switch(inp.comptypeid) {
          case 1:  convert_dtype((char*)output, (unsigned char*) input); break;
          case 2:  convert_dtype((char*)output, (char*) input); break;
          case 3:  convert_dtype((char*)output, (unsigned short*) input); break;
          case 4:  convert_dtype((char*)output, (short*) input); break;
          case 5:  convert_dtype((char*)output, (unsigned int*) input); break;
          case 6:  convert_dtype((char*)output, (int*) input); break;
          case 7:  convert_dtype((char*)output, (unsigned long*) input); break;
          case 8:  convert_dtype((char*)output, (long*) input); break;
          case 9:  convert_dtype((char*)output, (float*) input); break;
          case 10: convert_dtype((char*)output, (double*) input); break;
          default: assert(0);
        }
        break;
      }
      case 3:  // (unsigned short);
      {
        switch(inp.comptypeid) {
          case 1:  convert_dtype((unsigned short*)output, (unsigned char*) input); break;
          case 2:  convert_dtype((unsigned short*)output, (char*) input); break;
          case 3:  convert_dtype((unsigned short*)output, (unsigned short*) input); break;
          case 4:  convert_dtype((unsigned short*)output, (short*) input); break;
          case 5:  convert_dtype((unsigned short*)output, (unsigned int*) input); break;
          case 6:  convert_dtype((unsigned short*)output, (int*) input); break;
          case 7:  convert_dtype((unsigned short*)output, (unsigned long*) input); break;
          case 8:  convert_dtype((unsigned short*)output, (long*) input); break;
          case 9:  convert_dtype((unsigned short*)output, (float*) input); break;
          case 10: convert_dtype((unsigned short*)output, (double*) input); break;
          default: assert(0);
        }
        break;
      }
      case 4:  // (short);
      {
        switch(inp.comptypeid) {
          case 1:  convert_dtype((short*)output, (unsigned char*) input); break;
          case 2:  convert_dtype((short*)output, (char*) input); break;
          case 3:  convert_dtype((short*)output, (unsigned short*) input); break;
          case 4:  convert_dtype((short*)output, (short*) input); break;
          case 5:  convert_dtype((short*)output, (unsigned int*) input); break;
          case 6:  convert_dtype((short*)output, (int*) input); break;
          case 7:  convert_dtype((short*)output, (unsigned long*) input); break;
          case 8:  convert_dtype((short*)output, (long*) input); break;
          case 9:  convert_dtype((short*)output, (float*) input); break;
          case 10: convert_dtype((short*)output, (double*) input); break;
          default: assert(0);
        }
        break;
      }
      case 5:  // (unsigned int);
      {
        switch(inp.comptypeid) {
          case 1:  convert_dtype((unsigned int*)output, (unsigned char*) input); break;
          case 2:  convert_dtype((unsigned int*)output, (char*) input); break;
          case 3:  convert_dtype((unsigned int*)output, (unsigned short*) input); break;
          case 4:  convert_dtype((unsigned int*)output, (short*) input); break;
          case 5:  convert_dtype((unsigned int*)output, (unsigned int*) input); break;
          case 6:  convert_dtype((unsigned int*)output, (int*) input); break;
          case 7:  convert_dtype((unsigned int*)output, (unsigned long*) input); break;
          case 8:  convert_dtype((unsigned int*)output, (long*) input); break;
          case 9:  convert_dtype((unsigned int*)output, (float*) input); break;
          case 10: convert_dtype((unsigned int*)output, (double*) input); break;
          default: assert(0);
        }
        break;
      }
      case 6:  // (int);
      {
        switch(inp.comptypeid) {
          case 1:  convert_dtype((int*)output, (unsigned char*) input); break;
          case 2:  convert_dtype((int*)output, (char*) input); break;
          case 3:  convert_dtype((int*)output, (unsigned short*) input); break;
          case 4:  convert_dtype((int*)output, (short*) input); break;
          case 5:  convert_dtype((int*)output, (unsigned int*) input); break;
          case 6:  convert_dtype((int*)output, (int*) input); break;
          case 7:  convert_dtype((int*)output, (unsigned long*) input); break;
          case 8:  convert_dtype((int*)output, (long*) input); break;
          case 9:  convert_dtype((int*)output, (float*) input); break;
          case 10: convert_dtype((int*)output, (double*) input); break;
          default: assert(0);
        }
        break;
      }
      case 7:  // (unsigned long);
      {
        switch(inp.comptypeid) {
          case 1:  convert_dtype((unsigned long*)output, (unsigned char*) input); break;
          case 2:  convert_dtype((unsigned long*)output, (char*) input); break;
          case 3:  convert_dtype((unsigned long*)output, (unsigned short*) input); break;
          case 4:  convert_dtype((unsigned long*)output, (short*) input); break;
          case 5:  convert_dtype((unsigned long*)output, (unsigned int*) input); break;
          case 6:  convert_dtype((unsigned long*)output, (int*) input); break;
          case 7:  convert_dtype((unsigned long*)output, (unsigned long*) input); break;
          case 8:  convert_dtype((unsigned long*)output, (long*) input); break;
          case 9:  convert_dtype((unsigned long*)output, (float*) input); break;
          case 10: convert_dtype((unsigned long*)output, (double*) input); break;
          default: assert(0);
        }
        break;
      }
      case 8:  // (long);
      {
        switch(inp.comptypeid) {
          case 1:  convert_dtype((long*)output, (unsigned char*) input); break;
          case 2:  convert_dtype((long*)output, (char*) input); break;
          case 3:  convert_dtype((long*)output, (unsigned short*) input); break;
          case 4:  convert_dtype((long*)output, (short*) input); break;
          case 5:  convert_dtype((long*)output, (unsigned int*) input); break;
          case 6:  convert_dtype((long*)output, (int*) input); break;
          case 7:  convert_dtype((long*)output, (unsigned long*) input); break;
          case 8:  convert_dtype((long*)output, (long*) input); break;
          case 9:  convert_dtype((long*)output, (float*) input); break;
          case 10: convert_dtype((long*)output, (double*) input); break;
          default: assert(0);
        }
        break;
      }
      case 9:  // (float);
      {
        switch(inp.comptypeid) {
          case 1:  convert_dtype((float*)output, (unsigned char*) input); break;
          case 2:  convert_dtype((float*)output, (char*) input); break;
          case 3:  convert_dtype((float*)output, (unsigned short*) input); break;
          case 4:  convert_dtype((float*)output, (short*) input); break;
          case 5:  convert_dtype((float*)output, (unsigned int*) input); break;
          case 6:  convert_dtype((float*)output, (int*) input); break;
          case 7:  convert_dtype((float*)output, (unsigned long*) input); break;
          case 8:  convert_dtype((float*)output, (long*) input); break;
          case 9:  convert_dtype((float*)output, (float*) input); break;
          case 10: convert_dtype((float*)output, (double*) input); break;
          default: assert(0);
        }
        break;
      }
      case 10: // (double);
      {
        switch(inp.comptypeid) {
          case 1:  convert_dtype((double*)output, (unsigned char*) input); break;
          case 2:  convert_dtype((double*)output, (char*) input); break;
          case 3:  convert_dtype((double*)output, (unsigned short*) input); break;
          case 4:  convert_dtype((double*)output, (short*) input); break;
          case 5:  convert_dtype((double*)output, (unsigned int*) input); break;
          case 6:  convert_dtype((double*)output, (int*) input); break;
          case 7:  convert_dtype((double*)output, (unsigned long*) input); break;
          case 8:  convert_dtype((double*)output, (long*) input); break;
          case 9:  convert_dtype((double*)output, (float*) input); break;
          case 10: convert_dtype((double*)output, (double*) input); break;
          default: assert(0);
        }
        break;
      }
      default: assert(0);
    }
  }

  void assign(const triple<unsigned int> idim,
              const triple<double>       iorig,
              const triple<double>       ispacing,
              const short incomp,
              const char* datatype)
  {
    dim     = idim;
    orig    = iorig;
    pixspc  = ispacing;
    ncomp   = incomp;

    comptypeid = itk_get_datatype_id(datatype);
    compsz     = itk_get_datatype_size(comptypeid);

    npix   = dim.v1 * dim.v2 * dim.v3;
    pixsz  = ncomp * compsz;
    datasz = npix * pixsz;

    databuff.resize(datasz);
    databuff.zero();
  }



  /**
  * @brief Refine the itk image data.
  *
  * @param scale Each pixel is refined into "scale" many pixels in each dimension.
  */
  void refine(triple<short> scale)
  {
    mt_vector<char> tbuff(this->databuff);
    triple<unsigned int> odim = this->dim;

    this->dim.v1 *= scale.v1;
    this->dim.v2 *= scale.v2;
    this->dim.v3 *= scale.v3;

    this->pixspc.v1 /= scale.v1;
    this->pixspc.v2 /= scale.v2;
    this->pixspc.v3 /= scale.v3;

    this->npix = this->dim.v1 * this->dim.v2 * this->dim.v3;
    this->datasz = this->npix * this->pixsz;
    this->databuff.resize(this->datasz);

    for(unsigned int i = 0; i < odim.v1; i++)
      for(unsigned int j = 0; j < odim.v2; j++)
        for(unsigned int k = 0; k < odim.v3; k++)
        {
          char* input = tbuff.data() + (k*odim.v1*odim.v2 + j*odim.v1 + i) * this->pixsz;

          for(short i2=0; i2<scale.v1; i2++)
            for(short j2 = 0; j2 < scale.v2; j2++)
              for(short k2 = 0; k2 < scale.v3; k2++)
              {
                unsigned int i3 = i*scale.v1 + i2;
                unsigned int j3 = j*scale.v2 + j2;
                unsigned int k3 = k*scale.v3 + k2;
                char* output = this->databuff.data() +
                               (k3*this->dim.v1*this->dim.v2 + j3*this->dim.v1 + i3) * this->pixsz;
                memcpy(output, input, this->pixsz);
              }
        }
  }

  void bounding_box(triple<unsigned int> & from, triple<unsigned int> & to)
  {
    const char itkpixel_zero[8*ITK_MAX_PIXEL_COMP] = {0};
    from = this->dim;
    to.v1 = to.v2 = to.v3 = 0;
    for (unsigned int i = 0; i < this->dim.v1; i++)
    {
      for (unsigned int j = 0; j < this->dim.v2; j++)
      {
        for (unsigned int k = 0; k < this->dim.v3; k++)
        {
          char * pixel = this->databuff.data() + (k*this->dim.v1*this->dim.v2 + j*this->dim.v1 + i) * this->pixsz;
          if (memcmp(pixel, itkpixel_zero, this->pixsz*sizeof(char)) != 0)
          {
            if (i < from.v1) from.v1 = i;
            if (j < from.v2) from.v2 = j;
            if (k < from.v3) from.v3 = k;
            if (i > to.v1) to.v1 = i;
            if (j > to.v2) to.v2 = j;
            if (k > to.v3) to.v3 = k;
          }
        }
      }
    }
    printf("[ITKIMAGE]: bounding box ( %u %u %u ) - ( %u %u %u )\n", from.v1, from.v2, from.v3, to.v1, to.v2, to.v3);
  }

  /**
  * @brief extract slices with respect to a given plane
  */
  void extract_slices(unsigned int plane, const std::set<unsigned int> & idx)
  {
    if ((plane != 3) && (plane != 5) && (plane != 6))
      return;
    
    triple<unsigned int> newdim = this->dim;
    if (plane == 3) 
    {
      newdim.v3 = idx.size();
      if (newdim.v3 > this->dim.v3) return;
    }
    else if (plane == 5) 
    {
      newdim.v2 = idx.size();
      if (newdim.v2 > this->dim.v2) return;
    }
    else 
    {
      newdim.v1 = idx.size();
      if (newdim.v1 > this->dim.v1) return;
    }

    unsigned int newnpix = newdim.v1*newdim.v2*newdim.v3;
    unsigned int newdatasz = newnpix*this->pixsz;
    char *newdata = new char[newdatasz];
     
    if (plane == 3)
    {
      unsigned int ik, ok;
      for (unsigned int i = 0; i < this->dim.v1; i++)
      {
        for (unsigned int j = 0; j < this->dim.v2; j++)
        {
          for (auto it = idx.cbegin(); it != idx.cend(); it++)
          {
            ik = *it;
            ok = std::distance(idx.cbegin(), it);
            if (ik < this->dim.v3)
            {
              char *src  = this->databuff.data() + (ik*this->dim.v1*this->dim.v2 + j*this->dim.v1 + i)*this->pixsz;
              char *dest = newdata + (ok*newdim.v1*newdim.v2 + j*newdim.v1 + i)*this->pixsz;
              memcpy(dest, src, this->pixsz*sizeof(char));
            } 
          }
        }
      }
    }
    else if (plane == 5)
    {
      unsigned int ij, oj;
      for (unsigned int i = 0; i < this->dim.v1; i++)
      {
       for (auto it = idx.cbegin(); it != idx.cend(); it++)
        {
          for (unsigned int k = 0; k < this->dim.v3; k++)
          {
            ij = *it;
            oj = std::distance(idx.cbegin(), it);
            if (ij < this->dim.v2)
            {
              char *src  = this->databuff.data() + (k*this->dim.v1*this->dim.v2 + ij*this->dim.v1 + i)*this->pixsz;
              char *dest = newdata + (k*newdim.v1*newdim.v2 + oj*newdim.v1 + i)*this->pixsz;
              memcpy(dest, src, this->pixsz*sizeof(char));
            } 
          }
        }
      }
    }
    else
    {
      unsigned int ii, oi;
      for (auto it = idx.cbegin(); it != idx.cend(); it++)
      {
        for (unsigned int j = 0; j < this->dim.v2; j++)
        {
          for (unsigned int k = 0; k < this->dim.v3; k++)
          {
            ii = *it;
            oi = std::distance(idx.cbegin(), it);
            if (ii < this->dim.v1)
            {
              char *src  = this->databuff.data() + (k*this->dim.v1*this->dim.v2 + j*this->dim.v1 + ii)*this->pixsz;
              char *dest = newdata + (k*newdim.v1*newdim.v2 + j*newdim.v1 + oi)*this->pixsz;
              memcpy(dest, src, this->pixsz*sizeof(char));             
            } 
          }
        }
      }
    }

    this->dim = newdim;
    this->npix = newnpix;
    this->datasz = newdatasz;
    this->databuff.assign(newdatasz, newdata);
  }

  /**
  * @brief Normalize the voxel spacing and change image dimension to keep the image spacing
  */
  void normalize_spacing()
  {
    const double min_spacing = std::min(std::min(this->pixspc.v1, this->pixspc.v2), std::min(this->pixspc.v3, 1.0));
    this->pixspc.v1 /= min_spacing;
    this->pixspc.v2 /= min_spacing;
    this->pixspc.v3 /= min_spacing;
    printf("[ITKIMAGE]: normalize spacing, new spacing %lf\n", min_spacing);
    triple<unsigned int> newdim;
    unsigned int l, m, n;
    newdim.v1 = static_cast<unsigned int>(ceil(this->dim.v1*this->pixspc.v1));
    newdim.v2 = static_cast<unsigned int>(ceil(this->dim.v2*this->pixspc.v2));
    newdim.v3 = static_cast<unsigned int>(ceil(this->dim.v3*this->pixspc.v3));
    if ((1.0*newdim.v1-0.5) > (this->dim.v1*this->pixspc.v1)) newdim.v1--;
    if ((1.0*newdim.v2-0.5) > (this->dim.v2*this->pixspc.v2)) newdim.v2--;
    if ((1.0*newdim.v3-0.5) > (this->dim.v3*this->pixspc.v3)) newdim.v3--;
    unsigned int newnpix = newdim.v1*newdim.v2*newdim.v3;
    unsigned int newdatasz = newnpix*this->pixsz;
    char *newdata = new char[newdatasz];
    for (unsigned int i = 0; i < newdim.v1; i++)
    {
      for (unsigned int j = 0; j < newdim.v2; j++)
      {
        for (unsigned int k = 0; k < newdim.v3; k++)
        {
          char *dest = newdata + (k*newdim.v1*newdim.v2 + j*newdim.v1 + i)*this->pixsz;
          l = static_cast<unsigned int>(floor((1.0*i+0.5)/this->pixspc.v1));
          m = static_cast<unsigned int>(floor((1.0*j+0.5)/this->pixspc.v2));
          n = static_cast<unsigned int>(floor((1.0*k+0.5)/this->pixspc.v3));
          char *src = this->databuff.data() + (n*this->dim.v1*this->dim.v2 + m*this->dim.v1 + l)*this->pixsz;
          memcpy(dest, src, this->pixsz*sizeof(char));
        }
      }
    }
    this->dim = newdim;
    this->npix = newnpix;
    this->datasz = newdatasz;
    this->pixspc.v1 = this->pixspc.v2 = this->pixspc.v3 = min_spacing;
    this->databuff.assign(newdatasz, newdata);
  }

  /**
   * @brief Add padding to the voxel data
   */
  void padding(const triple<unsigned int> & size)
  {
    triple<unsigned int> newdim;
    newdim.v1 = this->dim.v1 + 2*size.v1;
    newdim.v2 = this->dim.v2 + 2*size.v2;
    newdim.v3 = this->dim.v3 + 2*size.v3;
    unsigned int newnpix = newdim.v1*newdim.v2*newdim.v3;
    unsigned int newdatasz = newnpix * this->pixsz;
    char *newdata = new char[newdatasz];
    for (unsigned int i = 0; i < this->dim.v1; i++)
    {
      for (unsigned int j = 0; j < this->dim.v2; j++)
      {
        for (unsigned int k = 0; k < this->dim.v3; k++)
        {
          char *src = this->databuff.data() + (k*this->dim.v1*this->dim.v2 + j*this->dim.v1 + i)*this->pixsz;
          char *dest = newdata + ((k+size.v3)*newdim.v1*newdim.v2 + (j+size.v2)*newdim.v1 + (i+size.v1))*this->pixsz;
          memcpy(dest, src, this->pixsz*sizeof(char));
        }
      }
    }
    this->dim = newdim;
    this->npix = newnpix;
    this->datasz = newdatasz;
    this->databuff.assign(newdatasz, newdata);
  }

  /**
   * @brief Crop pixel data
   */
  void crop()
  {
    triple<unsigned int> from, to, newdim;
    this->bounding_box(from, to);
    newdim.v1 = to.v1 - from.v1 + 1;
    newdim.v2 = to.v2 - from.v2 + 1;
    newdim.v3 = to.v3 - from.v3 + 1;
    unsigned int newnpix = newdim.v1*newdim.v2*newdim.v3;
    unsigned int newdatasz = newnpix * this->pixsz;
    char *newdata = new char[newdatasz];
    for (unsigned int i = 0; i < newdim.v1; i++)
    {
      for (unsigned int j = 0; j < newdim.v2; j++)
      {
        for (unsigned int k = 0; k < newdim.v3; k++)
        {
          char *src = this->databuff.data() + ((k+from.v3)*this->dim.v1*this->dim.v2 + (j+from.v2)*this->dim.v1 + (i+from.v1))*this->pixsz;
          char *dest = newdata + (k*newdim.v1*newdim.v2 + j*newdim.v2 + i)*this->pixsz;
          memcpy(dest, src, this->pixsz*sizeof(char));
        }
      }
    }
    this->dim = newdim;
    this->npix = newnpix;
    this->datasz = newdatasz;
    this->databuff.assign(newdatasz, newdata);
  }

  void flip(unsigned int axes)
  {
    unsigned int ni, nj, nk;
    char *newdata = new char[this->datasz];
    for (unsigned int i = 0; i < this->dim.v1; i++)
    {
      for (unsigned int j = 0; j < this->dim.v2; j++)
      {
        for (unsigned int k = 0; k < this->dim.v3; k++)
        {
          char *src = this->databuff.data() + (k*this->dim.v1*this->dim.v2 + j*this->dim.v1 + i)*this->pixsz;
          ni = (((axes & ITK_AXIS_X) == ITK_AXIS_X) ? (this->dim.v1-i-1) : i);
          nj = (((axes & ITK_AXIS_Y) == ITK_AXIS_Y) ? (this->dim.v2-j-1) : j);
          nk = (((axes & ITK_AXIS_Z) == ITK_AXIS_Z) ? (this->dim.v3-k-1) : k);
          char *dest = newdata + (nk*this->dim.v1*this->dim.v2 + nj*this->dim.v1 + ni)*this->pixsz;
          memcpy(dest, src, this->pixsz*sizeof(char));
        }
      }
    }
    this->databuff.assign(datasz, newdata);
  }

  void shift(const triple<unsigned int> & offset)
  {
    unsigned int ni, nj, nk;
    char *newdata = new char[this->datasz];
    for (unsigned int i = 0; i < this->dim.v1; i++)
    {
      for (unsigned int j = 0; j < this->dim.v2; j++)
      {
        for (unsigned int k = 0; k < this->dim.v3; k++)
        {
          char *src = this->databuff.data() + (k*this->dim.v1*this->dim.v2 + j*this->dim.v1 + i)*this->pixsz;
          ni = (i + offset.v1) % this->dim.v1;
          nj = (j + offset.v2) % this->dim.v2;
          nk = (k + offset.v3) % this->dim.v3;
          char *dest = newdata + (nk*this->dim.v1*this->dim.v2 + nj*this->dim.v1 + ni)*this->pixsz;
          memcpy(dest, src, this->pixsz*sizeof(char));
        }
      }
    }
    this->databuff.assign(datasz, newdata);
  }

};

/**
* @brief Class used for accessing itk image data.
*
* @tparam T image data type.
*/
template<typename T>
class itk_access
{
  private:
  itk_image & _img;

  public:
  triple<unsigned int> & dim;
  triple<double> & pixspc;
  triple<float> spc_fct;

  /// constructor
  itk_access(itk_image & img) : _img(img), dim(_img.dim), pixspc(_img.pixspc)
  {
    spc_fct.v1 = _img.pixspc.v1 / _img.pixspc.v1;
    spc_fct.v2 = _img.pixspc.v1 / _img.pixspc.v2;
    spc_fct.v3 = _img.pixspc.v1 / _img.pixspc.v3;
  }


  /**
  * @brief Data access function
  *
  * @param [in] i    Discrete coordinate in x.
  * @param [in] j    Discrete coordinate in y.
  * @param [in] k    Discrete coordinate in z.
  * @param [in] c    Optional component index.
  *
  * @return Reference to data value.
  */
  T & operator()(unsigned int i, unsigned int j, unsigned int k, unsigned int c = 0)
  {
    T* data = (T*) _img.databuff.data();
    return (data + (k * dim.v1 * dim.v2 + j * dim.v1 + i) * _img.ncomp)[c];
  }

  /**
  * @brief Compute physical distance for two pixels.
  *
  * @param i1 Pixel 1, discrete coordinate in x.
  * @param j1 Pixel 1, discrete coordinate in y.
  * @param k1 Pixel 1, discrete coordinate in z.
  * @param i2 Pixel 2, discrete coordinate in x.
  * @param j2 Pixel 2, discrete coordinate in y.
  * @param k2 Pixel 2, discrete coordinate in z.
  *
  * @return Distance between the pixels.
  */
  double distance(unsigned int i1, unsigned int j1, unsigned int k1,
                  unsigned int i2, unsigned int j2, unsigned int k2)
  {
    double len_x = i1 * pixspc.v1 - i2 * pixspc.v1;
    double len_y = j1 * pixspc.v2 - j2 * pixspc.v2;
    double len_z = k1 * pixspc.v3 - k2 * pixspc.v3;

    return sqrt(len_x*len_x + len_y*len_y + len_z*len_z);
  }

  bool is_in_ellipsoid(unsigned int i1, unsigned int j1, unsigned int k1,
                       unsigned int i2, unsigned int j2, unsigned int k2,
                       double rada, double radb, double radc)
  {
    double len_x = (i1 * pixspc.v1 - i2 * pixspc.v1)/rada;
    double len_y = (j1 * pixspc.v2 - j2 * pixspc.v2)/radb;
    double len_z = (k1 * pixspc.v3 - k2 * pixspc.v3)/radc;
    return ((len_x*len_x + len_y*len_y + len_z*len_z) < 1.0);
  }


  template<class V>
  mt_point<V> position(unsigned int i, unsigned int j, unsigned int k)
  {
    mt_point<V> pos;
    pos.x = _img.orig.v1 + i * pixspc.v1 + 0.5*pixspc.v1;
    pos.y = _img.orig.v2 + j * pixspc.v2 + 0.5*pixspc.v2;
    pos.z = _img.orig.v3 + k * pixspc.v3 + 0.5*pixspc.v3;

    return pos;
  }

  template<class V>
  vec3i pixel(mt_point<V> pt)
  {
    float i = (pt.x - _img.orig.v1 - 0.5*pixspc.v1) / pixspc.v1;
    float j = (pt.y - _img.orig.v2 - 0.5*pixspc.v2) / pixspc.v2;
    float k = (pt.z - _img.orig.v3 - 0.5*pixspc.v3) / pixspc.v3;

    vec3i ijk = {(mt_int)i, (mt_int)j, (mt_int)k};
    return ijk;
  }
};


/**
* @brief Check if pixel has any non-zero data in a certain radius
*
* @tparam T Image data type.
*
* @param [in] data Data access.
* @param [in] rad  Radius around pixel we check for.
* @param [in] i    Discrete coordinate in x.
* @param [in] j    Discrete coordinate in y.
* @param [in] k    Discrete coordinate in z.
*
* @return Whether there is non-zero data.
*/
template<typename T>
bool pix_has_data(itk_access<T> & data, triple<short> rad, unsigned i, unsigned j, unsigned k)
{
  unsigned int rad_i = rad.v1 * data.spc_fct.v1;
  unsigned int rad_j = rad.v2 * data.spc_fct.v2;
  unsigned int rad_k = rad.v3 * data.spc_fct.v3; 
  const unsigned int i_start = int(i - rad_i) > 0 ? i - rad_i : 0;
  const unsigned int j_start = int(j - rad_j) > 0 ? j - rad_j : 0;
  const unsigned int k_start = int(k - rad_k) > 0 ? k - rad_k : 0;
  rad_i += 1, rad_j += 1, rad_k += 1;
  const unsigned int i_end = i + rad_i < data.dim.v1 ? i + rad_i : data.dim.v1;
  const unsigned int j_end = j + rad_j < data.dim.v2 ? j + rad_j : data.dim.v2;
  const unsigned int k_end = k + rad_k < data.dim.v3 ? k + rad_k : data.dim.v3;

  for(unsigned int i2 = i_start; i2 < i_end; i2++)
    for(unsigned int j2 = j_start; j2 < j_end; j2++)
      for(unsigned int k2 = k_start; k2 < k_end; k2++)
        if(data(i2, j2, k2) != T(0))
          return true;

  return false;
}

template <typename T>
bool pix_has_data(itk_access<T> & data,  unsigned i, unsigned j, unsigned k)
{
  return (data(i,j,k) != T(0));
}

/**
* @brief Check if pixel data in a certain radius consists only of non-zero vals
*
* @tparam T Image data type.
*
* @param [in] data Data access.
* @param [in] rad  Radius around pixel we check for.
* @param [in] i    Discrete coordinate in x.
* @param [in] j    Discrete coordinate in y.
* @param [in] k    Discrete coordinate in z.
*
* @return Whether there is non-zero data.
*/
template<typename T>
bool pix_full_data(itk_access<T> & data, triple<short> rad, unsigned i, unsigned j, unsigned k)
{
  unsigned int rad_i = rad.v1 * data.spc_fct.v1;
  unsigned int rad_j = rad.v2 * data.spc_fct.v2;
  unsigned int rad_k = rad.v3 * data.spc_fct.v3;
  const unsigned int i_start = int(i - rad_i) > 0 ? i - rad_i : 0;
  const unsigned int j_start = int(j - rad_j) > 0 ? j - rad_j : 0;
  const unsigned int k_start = int(k - rad_k) > 0 ? k - rad_k : 0;
  rad_i += 1, rad_j += 1, rad_k += 1;
  const unsigned int i_end = i + rad_i < data.dim.v1 ? i + rad_i : data.dim.v1;
  const unsigned int j_end = j + rad_j < data.dim.v2 ? j + rad_j : data.dim.v2;
  const unsigned int k_end = k + rad_k < data.dim.v3 ? k + rad_k : data.dim.v3;

  for(unsigned int i2 = i_start; i2 < i_end; i2++)
    for(unsigned int j2 = j_start; j2 < j_end; j2++)
      for(unsigned int k2 = k_start; k2 < k_end; k2++)
        if(data(i2, j2, k2) == T(0))
          return false;

  return true;
}

/**
* @brief Smooth (average) the pixel data, taking a certain radius into account.
*
* @tparam T Image data type.
*
* @param [in] data Data access.
* @param [in] rad  Radius around pixel.
* @param [in] smth Smoothing coefficient.
* @param [in] i    Discrete coordinate in x.
* @param [in] j    Discrete coordinate in y.
* @param [in] k    Discrete coordinate in z.
*
* @return
*/
template<typename T>
void pix_smooth(itk_access<T> & data, triple<short> rad, float smth, unsigned i, unsigned j, unsigned k)
{
  unsigned int numadd = 0;
  T avrg = 0.0f;
  T curr = data(i,j,k);

  unsigned int rad_i = rad.v1 * data.spc_fct.v1;
  unsigned int rad_j = rad.v2 * data.spc_fct.v2;
  unsigned int rad_k = rad.v3 * data.spc_fct.v3;
  const unsigned int i_start = int(i - rad_i) > 0 ? i - rad_i : 0;
  const unsigned int j_start = int(j - rad_j) > 0 ? j - rad_j : 0;
  const unsigned int k_start = int(k - rad_k) > 0 ? k - rad_k : 0;
  rad_i += 1, rad_j += 1, rad_k += 1;
  const unsigned int i_end = i + rad_i < data.dim.v1 ? i + rad_i : data.dim.v1;
  const unsigned int j_end = j + rad_j < data.dim.v2 ? j + rad_j : data.dim.v2;
  const unsigned int k_end = k + rad_k < data.dim.v3 ? k + rad_k : data.dim.v3;

  #ifdef ITKSMOOTH_CHECK_RAD
  triple<double> space_rad = {rad_i*data.pixspc.v1, rad_j*data.pixspc.v2, rad_k*data.pixspc.v3};
  #endif

  for(unsigned int k2 = k_start; k2 < k_end; k2++)
    for(unsigned int j2 = j_start; j2 < j_end; j2++)
      for(unsigned int i2 = i_start; i2 < i_end; i2++)
      {
        #ifdef ITKSMOOTH_CHECK_RAD
        double len_x = (i*data.pixspc.v1 - i2*data.pixspc.v1)/space_rad.v1;
        double len_y = (j*data.pixspc.v2 - j2*data.pixspc.v2)/space_rad.v2;
        double len_z = (k*data.pixspc.v3 - k2*data.pixspc.v3)/space_rad.v3;
        if((len_x*len_x + len_y*len_y + len_z*len_z) < 1.0)
        #endif
        {
          avrg += data(i2,j2,k2);
          numadd++;
        }
      }

  avrg /= numadd;
  data(i,j,k) = curr + (avrg - curr)*smth;
}

/**
* @brief Smooth an itk image.
*
* The data has to be of float type.
*
* @param [in, out] img  ITK image.
* @param [in]      rad  Smoothing radius.
* @param [in]      smth Smoothing coefficent.
* @param [in]      iter Number of iterations.
*/
void smooth_itk(itk_image & img, triple<short> rad, float smth, unsigned int iter)
{
  itk_access<float> data(img);
  std::list<triple<unsigned int> > smth_list;

  for(unsigned int k = 0; k < img.dim.v3; k++)
    for(unsigned int j = 0; j < img.dim.v2; j++)
      for(unsigned int i = 0; i < img.dim.v1; i++)
        if(pix_has_data(data, rad, i, j, k)) smth_list.push_back( {i, j, k} );

  PROGRESS<unsigned int> prg(iter, "Smooting progress: ");

  std::vector< triple<unsigned int> > work;
  work.assign(smth_list.begin(), smth_list.end());

  for(unsigned int it = 0; it < iter; it++)
  {
    #ifdef OPENMP
    #pragma omp parallel for schedule(guided, 10)
    #endif
    for(size_t w = 0; w < work.size(); w++) {
      triple<unsigned int> & s = work[w];
      pix_smooth(data, rad, smth, s.v1, s.v2, s.v3);
    }

    #ifdef OPENMP
    #pragma omp parallel for schedule(guided, 10)
    #endif
    for(size_t w = 0; w < work.size(); w++) {
      triple<unsigned int> & s = work[w];
      pix_smooth(data, rad, -smth * SMOOTH_FREQ_SCA, s.v1, s.v2, s.v3);
    }

    prg.next();
  }

  prg.finish();
}

void itk_mask(itk_image & img, float val, float zeroval = 0.0f)
{
  itk_access<float> data(img);

  for(unsigned int k = 0; k < img.dim.v3; k++)
    for(unsigned int j = 0; j < img.dim.v2; j++)
      for(unsigned int i = 0; i < img.dim.v1; i++)
        if(data(i,j,k) != val) data(i,j,k) = zeroval;
}


void itk_extract_vals(itk_image & in,
                      itk_image & out,
                      const std::set<float> & vals,
                      float zeroval = 0.0f)
{
  assert(in.dim.v1 == out.dim.v1 &&
         in.dim.v2 == out.dim.v2 && in.dim.v3 == out.dim.v3);

  itk_access<float> input(in);
  itk_access<float> output(out);

  for(unsigned int k = 0; k < in.dim.v3; k++)
    for(unsigned int j = 0; j < in.dim.v2; j++)
      for(unsigned int i = 0; i < in.dim.v1; i++) {
        float v = input(i,j,k);
        if(vals.count(v)) {
          input(i,j,k)  = zeroval;
          output(i,j,k) = v;
        }
      }
}

void itk_insert_vals(itk_image & in, itk_image & out, const std::set<float> & vals)
{
  assert(in.dim.v1 == out.dim.v1 &&
         in.dim.v2 == out.dim.v2 && in.dim.v3 == out.dim.v3);

  itk_access<float> input(in);
  itk_access<float> output(out);

  for(unsigned int k = 0; k < in.dim.v3; k++)
    for(unsigned int j = 0; j < in.dim.v2; j++)
      for(unsigned int i = 0; i < in.dim.v1; i++) {
        float v = input(i,j,k);
        if(vals.size() == 0 || vals.count(v))
          output(i,j,k) = v;
      }
}

void itk_set_val(itk_image & img, float val)
{
  itk_access<float> data(img);

  for(unsigned int k = 0; k < img.dim.v3; k++)
    for(unsigned int j = 0; j < img.dim.v2; j++)
      for(unsigned int i = 0; i < img.dim.v1; i++)
        data(i,j,k) = val;
}

void itk_dilate(itk_image & img, triple<short> rad, float val)
{
  itk_access<float> data(img);

  #ifdef OPENMP
  #pragma omp parallel
  #endif
  {
    std::vector<triple<unsigned int> > set_list; set_list.reserve(10000);

    #ifdef OPENMP
    #pragma omp for schedule(guided, 10)
    #endif
    for(unsigned int k = 0; k < img.dim.v3; k++)
      for(unsigned int j = 0; j < img.dim.v2; j++)
        for(unsigned int i = 0; i < img.dim.v1; i++)
          if(pix_has_data(data, rad, i, j, k))
            set_list.push_back( {i, j, k} );

    #ifdef OPENMP
    #pragma omp barrier
    #endif

    for(const triple<unsigned int> & t : set_list)
      data(t.v1, t.v2, t.v3) = val;
  }
}

void itk_erode(itk_image & img, triple<short> rad)
{
  itk_access<float> data(img);
  #ifdef OPENMP
  #pragma omp parallel
  #endif
  {
    std::vector<triple<unsigned int> > rem_list; rem_list.reserve(10000);

    #ifdef OPENMP
    #pragma omp for schedule(guided, 10)
    #endif
    for(unsigned int k = 0; k < img.dim.v3; k++)
      for(unsigned int j = 0; j < img.dim.v2; j++)
        for(unsigned int i = 0; i < img.dim.v1; i++)
          if(!pix_full_data(data, rad, i,j,k))
            rem_list.push_back( {i,j,k} );

    #ifdef OPENMP
    #pragma omp barrier
    #endif

    for(const triple<unsigned int> & t : rem_list)
      data(t.v1, t.v2, t.v3) = 0.0f;
  }
}



/**
* @brief Clamp to a value based on minimal distance.
*
* @param [in] curr        Current value.
* @param [in] candidates  Set of value candidates.
*
* @return The value we clamped to.
*/
template<typename T, typename S>
T clamp_min_dist(T curr, const std::set<S> & candidates)
{
  auto it = candidates.begin();

  T dist = fabs(T(*it) - curr);
  T choice = T(*it);
  ++it;

  while(it != candidates.end())
  {
    T f = fabs(T(*it) - curr);
    if(dist > f) {
      dist = f;
      choice = T(*it);
    }
    ++it;
  }

  return choice;
}


/**
* @brief Clamp data of one pixel based on reference data in a radius.
*
* @param [in] cdata  Clamping data.
* @param [in] rdata  Reference data.
* @param [in] rad    Radius of interest.
* @param [in] i      Discrete x coordinate of pixel.
* @param [in] j      Discrete y coordinate of pixel.
* @param [in] k      Discrete z coordinate of pixel.
*/
template<typename T, typename S>
void pix_clamp(itk_access<T> & cdata, itk_access<S> & rdata, triple<short> rad, unsigned i, unsigned j, unsigned k)
{
  std::set<S> candidates;

  unsigned int rad_i = rad.v1 * rdata.spc_fct.v1;
  unsigned int rad_j = rad.v2 * rdata.spc_fct.v2;
  unsigned int rad_k = rad.v3 * rdata.spc_fct.v3;
  const unsigned int i_start = int(i - rad_i) > 0 ? i - rad_i : 0;
  const unsigned int j_start = int(j - rad_j) > 0 ? j - rad_j : 0;
  const unsigned int k_start = int(k - rad_k) > 0 ? k - rad_k : 0;
  rad_i += 1, rad_j += 1, rad_k += 1;
  const unsigned int i_end = i + rad_i < rdata.dim.v1 ? i + rad_i : rdata.dim.v1;
  const unsigned int j_end = j + rad_j < rdata.dim.v2 ? j + rad_j : rdata.dim.v2;
  const unsigned int k_end = k + rad_k < rdata.dim.v3 ? k + rad_k : rdata.dim.v3;

  triple<double> space_rad = {rad_i*rdata.pixspc.v1, rad_j*rdata.pixspc.v2, rad_k*rdata.pixspc.v3};

  // pixels of reference data in a given radius are selected as candidates
  for(unsigned int k2 = k_start; k2 < k_end; k2++)
    for(unsigned int j2 = j_start; j2 < j_end; j2++)
      for(unsigned int i2 = i_start; i2 < i_end; i2++)
        if(rdata.is_in_ellipsoid(i,j,k,i2,j2,k2,space_rad.v1,space_rad.v2,space_rad.v3))
          candidates.insert(rdata(i2,j2,k2));

  T curr = cdata(i,j,k);
  if(candidates.size() > 0)
    cdata(i,j,k) = clamp_min_dist(curr, candidates);
  else {
    std::cerr << "Clamping error: No candidates!" << std::endl;
    exit(1);
  }
}


/**
* @brief Clamp the pixel data of an itk image.
*
* @param [in] cimg   Image of floating point data to clamp.
* @param [in] rimg   Image of integer reference data to clamp to.
* @param [in] rad    Radius of interest.
*/
void clamp_itk(itk_image & cimg, itk_image & rimg, triple<short> rad)
{
  itk_access<float> cdata(cimg);  // data to clamp
  itk_access<short> rdata(rimg);  // reference data

  std::list<triple<unsigned int> > clamp_list;

  // all pixels containing non-zero data are clamped
  for(unsigned int k = 0; k < cimg.dim.v3; k++)
    for(unsigned int j = 0; j < cimg.dim.v2; j++)
      for(unsigned int i = 0; i < cimg.dim.v1; i++)
        if(pix_has_data(cdata, i, j, k)) clamp_list.push_back( {i, j, k} );

  PROGRESS<unsigned int> prg(clamp_list.size(), "Clamping progress: ");

  for(auto c = clamp_list.begin(); c != clamp_list.end(); ++c) {
    pix_clamp(cdata, rdata, rad, c->v1, c->v2, c->v3);
    prg.next();
  }
  prg.finish();
}

#endif
