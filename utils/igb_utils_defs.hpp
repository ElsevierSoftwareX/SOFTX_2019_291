/**
* @file igb_header_utils.h
* @brief IGB format specific IO utils.
* @author Elias Karabelas
* @version
* @date 2016-12-13
*/

#ifndef _IGB_HEADER_UTILS
#define _IGB_HEADER_UTILS

#include <cmath>
#include <limits>
#include <type_traits>
#include <algorithm>


//Check for almost equal of a type

template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
    almost_equal(T x, T y, int ulp = 4)
{
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x-y) <= std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
        // unless the result is subnormal
        || std::abs(x-y) < std::numeric_limits<T>::min();
}

const char* Header_type[IGB_MAX_TYPE+1] =
{
  "",
  "byte",
  "char",
  "short",
  "long",
  "float",
  "double",
  "complex",
  "double_complex",
  "rgba",
  "structure",
  "pointer",
  "list",
  "int",
  "uint",
  "ushort",
  "vec3f",
  "vec3d",
  "vec4f",
  "vec4d",
  "vec9f"
};

//* size of the stored data, not the variable type
size_t Data_Size[IGB_MAX_TYPE+1] =
{
  0,
  sizeof(byte),
  sizeof(char),
  sizeof(short),
  sizeof(long),
  sizeof(float),
  sizeof(double),
  0,
  0,
  0,
  0,
  sizeof(void*),
  0,
  sizeof(int),
  sizeof(unsigned int),
  sizeof(unsigned short),
  3*sizeof(float),
  3*sizeof(double),
  4*sizeof(float),
  4*sizeof(double),
  9*sizeof(float)
};

/** the number of components for each data type */
const long Num_Components[IGB_MAX_TYPE+1] =
{
  0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 3, 3, 4, 4, 9
};


const long Header_Systeme_No[2] =
{
  IGB_BIG_ENDIAN,
  IGB_LITTLE_ENDIAN
};

const char* Header_Systeme[2] =
{
  "big_endian",
  "little_endian"
};

//maximum allowed number of keys, at the moment 36
const char* allowed_keys[] =
{
  "x", "y", "z", "t", "type", "systeme", "taille", "bin", "trame", "num", "comp", "lut", "dim_x", "dim_y",
  "dim_z", "dim_t", "inc_x", "inc_y", "inc_z", "inc_t", "org_x", "org_y", "org_z", "org_t", "unites_x", "unites_y",
  "unites_z", "unites_t", "epais", "unites", "facteur", "zero", "struct", "aut", "transparent"
};

//constant mapping for the keywords in the vector
const long allowed_keys_in_numbers[] =
{
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34
};

void make_igb_consistent(igb_header & igb_head)
{
  if(igb_head.bool_dim_x && igb_head.bool_inc_x && igb_head.bool_x)
    if(!almost_equal(igb_head.v_dim_x, igb_head.v_inc_x * (float)((igb_head.v_x -1))))
    {
      fprintf(stderr, "Adjusting dim_x to make dimensions consistent\n");
      igb_head.v_dim_x = igb_head.v_inc_x * (float)(igb_head.v_x-1);
    }
  if(igb_head.bool_dim_y && igb_head.bool_inc_y && igb_head.bool_y)
    if(!almost_equal(igb_head.v_dim_y,igb_head.v_inc_y * (float)(igb_head.v_y -1)))
    {
      fprintf(stderr, "Adjusting dim_y to make dimensions consistent\n");
      igb_head.v_dim_y = igb_head.v_inc_y * (float)(igb_head.v_y-1);
    }
  if(igb_head.bool_dim_z && igb_head.bool_inc_z && igb_head.bool_z)
    if(!almost_equal(igb_head.v_dim_z,igb_head.v_inc_z * (float)(igb_head.v_z -1)))
    {
      fprintf(stderr, "Adjusting dim_z to make dimensions consistent\n");
      igb_head.v_dim_z = igb_head.v_inc_z * (float)(igb_head.v_z-1);
    }
  if(igb_head.bool_dim_t && igb_head.bool_inc_t && igb_head.bool_t)
    if(!almost_equal(igb_head.v_dim_t, igb_head.v_inc_t * (float)(igb_head.v_t -1)))
    {
      fprintf(stderr, "Adjusting dim_t to make dimensions consistent\n");
      igb_head.v_dim_t = igb_head.v_inc_t * (float)(igb_head.v_t-1);
    }
}

long igb_system_endian_code()
{
  float val=IGB_ENDIAN_VAL;
  // quick fix.
  char  le_val[] = { (char)0,(char)48,(char)242,(char)199 },*pval   = (char *)(&val);
  assert( sizeof(float) == 4 );
  if ( *pval == le_val[0] )
    return IGB_LITTLE_ENDIAN;
  else
    return IGB_BIG_ENDIAN;
}

void igb_treat_file_open_error(std::string file, int errnum){
    fprintf(stderr, "An IO error occured when opening file %s: %s\n\n",
            file.c_str(), strerror(errnum));
    exit(1);
}


long init_igb_header(const std::string filename, igb_header & igb_head)
{
  long ret = 0;
  //Check file extension

  if(!(filename.substr(filename.find_last_of(".") + 1) == "dynpt")) //no dynpt extension
  {
    //maybe igb extension
    if(!(filename.substr(filename.find_last_of(".") + 1) == "igb"))
    {
      ret = -1;
      fprintf(stderr, "Error in %s: Cannot read %s\n", __func__, filename.c_str());
      return ret;
    }
  }

  igb_head.filename = filename;
  igb_head.fileptr  = NULL;

  igb_head.v_x = igb_head.v_y = igb_head.v_type = 0 ;
  igb_head.v_systeme = igb_system_endian_code();
  igb_head.v_comp = igb_head.v_lut = igb_head.v_num = igb_head.v_bin = 0 ;
  igb_head.v_trame = CC8 ;
  igb_head.v_z = igb_head.v_t = 1 ;
  igb_head.v_epais = 0.0 ;
  igb_head.v_inc_x = igb_head.v_inc_y = igb_head.v_inc_z = igb_head.v_inc_t = 1.0 ;
  igb_head.v_org_x = igb_head.v_org_y = igb_head.v_org_z = 1;
  igb_head.v_org_t = 0.0 ;
  igb_head.v_unites_x[0] = igb_head.v_unites_x[40] = '\000' ;
  igb_head.v_unites_y[0] = igb_head.v_unites_y[40] = '\000' ;
  igb_head.v_unites_z[0] = igb_head.v_unites_z[40] = '\000' ;
  igb_head.v_unites_t[0] = igb_head.v_unites_t[40] = '\000' ;
  igb_head.v_unites[0] = igb_head.v_unites[40] = '\000' ;
  igb_head.v_facteur = 1.0 ;
  igb_head.v_zero = 0.0 ;
  igb_head.v_aut_name[0] = igb_head.v_aut_name[40] = '\000' ;
  igb_head.v_struct_desc[0] = igb_head.v_struct_desc[40] = '\000' ;
  igb_head.v_transparent = NULL;

  igb_head.bool_x = igb_head.bool_y = igb_head.bool_type = false;
  igb_head.bool_z = igb_head.bool_t = false;
  igb_head.bool_taille = false;
  igb_head.bool_num = false;
  igb_head.bool_bin = false;
  igb_head.bool_trame = false;
  igb_head.bool_lut = false;
  igb_head.bool_comp = false;
  igb_head.bool_epais = false;
  igb_head.bool_org_x = igb_head.bool_org_y = igb_head.bool_org_z = igb_head.bool_org_t = false;
  igb_head.bool_inc_x = igb_head.bool_inc_y = igb_head.bool_inc_z = igb_head.bool_inc_t = false;
  igb_head.bool_dim_x = igb_head.bool_dim_y = igb_head.bool_dim_z = igb_head.bool_dim_t = false;
  igb_head.bool_unites_x = igb_head.bool_unites_y = igb_head.bool_unites_z = igb_head.bool_unites_t = false;
  igb_head.bool_unites = false;
  igb_head.bool_facteur = igb_head.bool_zero = false;
  igb_head.bool_struct_desc = false;
  igb_head.bool_aut_name = false;
  igb_head.bool_transparent = false;
  igb_head.igb_header_initialized = true;

  return(ret);
}

void set_igb_header_datatype(const std::string datatype, igb_header & igb_head)
{
  long tn = IGB_MIN_TYPE;
  while( tn <= IGB_MAX_TYPE && datatype.compare(Header_type[tn]))
  {
    tn++;
    if(tn <= IGB_MAX_TYPE)
      igb_head.v_type = tn;
    else {
      fprintf(stderr, "ERROR in %s: Illegal data type \"%s\" specified for IGB header\n",__func__, datatype.c_str());
      exit(-1);
    }
    igb_head.bool_type = true;
  }
}

void set_igb_header_systeme(const std::string datatype, igb_header & igb_head)
{
  if(datatype.compare(Header_Systeme[0]) == 0) //BIG ENDIAN
    igb_head.v_systeme = Header_Systeme_No[0];
  else if (datatype.compare(Header_Systeme[1]) == 0) //LITTLE ENDIAN
    igb_head.v_systeme = Header_Systeme_No[1];
  else
  {
    fprintf(stderr, "ERROR in %s: Illegal system \"%s\" specified for IGB header\n", __func__, datatype.c_str());
    exit(-1);
  }
}

long read_igb_header(igb_header & igb_head)
{
  long ret = 0;
  if(!igb_head.igb_header_initialized)
  {
    fprintf(stderr, "ERROR in %s: igb header was not initialized!", __func__);
    ret = -1;
    return ret;
  }

  char * header = new char[MY_IGB_HEADERSIZE+1];
  igb_head.fileptr = fopen(igb_head.filename.c_str(), "rb");
  if(igb_head.fileptr == NULL) igb_treat_file_open_error(igb_head.filename, errno);
  size_t chk = fread(header, sizeof(char), MY_IGB_HEADERSIZE, igb_head.fileptr);

  if(chk != MY_IGB_HEADERSIZE) {
    fprintf(stderr, "%s error: igb header size does not match! Aborting!", __func__);
    exit(1);
  }

  header[MY_IGB_HEADERSIZE] = '\0';
  //parse the header file into a key value string map
  char* pch = strtok (header, " \r\n");
  std::map<std::string, std::string> token_key_list;
  while(pch != nullptr)
  {
    std::string token = std::string(pch);
    std::string key, value;
    std::string::size_type pos = token.find(':');
    if(token.npos != pos) {
      value = token.substr(pos +1);
      key = token.substr(0,pos);
    }
    token_key_list.insert(std::make_pair(key,value));
    pch = strtok(nullptr, " \r\n");
  }

  delete[] header;

  //now parse the map
  //loop over the allowed keys
  for(size_t i=0; i < sizeof(allowed_keys) / sizeof(char*); i++)
  {
    std::map<std::string, std::string>::iterator search = token_key_list.find(allowed_keys[i]);
    if(search != token_key_list.end()) //found the token i in the map
    {
      //parse the stuff
      switch(i)
      {
        case 0 : igb_head.v_x = (long)(std::stoi(search->second)); igb_head.bool_x = true; break;
        case 1 : igb_head.v_y = (long)(std::stoi(search->second)); igb_head.bool_y = true; break;
        case 2 : igb_head.v_z = (long)(std::stoi(search->second)); igb_head.bool_z = true; break;
        case 3 : igb_head.v_t = (long)(std::stoi(search->second)); igb_head.bool_t = true; break;
        case 4 : set_igb_header_datatype(search->second, igb_head); break;
        case 5 : set_igb_header_systeme(search->second, igb_head); break;
        case 6 : igb_head.v_taille = (long)(std::stoi(search->second)); igb_head.bool_taille = true; break;
        case 7 : igb_head.v_bin = (long)(std::stoi(search->second)); igb_head.bool_bin = true; break;
        case 8 :
                 {
                   if(search->second.compare("c8") == 0)
                     igb_head.v_trame = CC8;
                   else if(search->second.compare("c4") == 0)
                     igb_head.v_trame = CC4;
                   else if(search->second.compare("hex") == 0)
                     igb_head.v_trame = HEX;
                   else if(search->second.compare("hexedges") == 0)
                     igb_head.v_trame = HEXEDGES;
                   else if(search->second.compare("hexbridges") == 0)
                     igb_head.v_trame = HEXBRIDGES;
                   else if(search->second.compare("hexlines") == 0)
                     igb_head.v_trame = HEXLINES;
                   else if(search->second.compare("hex2") == 0)
                     igb_head.v_trame = HEX2;
                   igb_head.bool_trame = true;
                   break;
                 }
        case 9 : igb_head.v_num   = (long)(std::stoi(search->second)); igb_head.bool_num = true; break;
        case 10 : igb_head.v_comp = (long)(std::stoi(search->second)); igb_head.bool_comp = true; break;
        case 11 : igb_head.v_lut  = (long)(std::stoi(search->second)); igb_head.bool_lut = true; break;
        case 12 : igb_head.v_dim_x = std::stof(search->second); igb_head.bool_dim_x = true; break;
        case 13 : igb_head.v_dim_y = std::stof(search->second); igb_head.bool_dim_y = true; break;
        case 14 : igb_head.v_dim_z = std::stof(search->second); igb_head.bool_dim_z = true; break;
        case 15 : igb_head.v_dim_t = std::stof(search->second); igb_head.bool_dim_t = true; break;
        case 16 : igb_head.v_inc_x = std::stof(search->second); igb_head.bool_inc_x = true; break;
        case 17 : igb_head.v_inc_y = std::stof(search->second); igb_head.bool_inc_y = true; break;
        case 18 : igb_head.v_inc_z = std::stof(search->second); igb_head.bool_inc_z = true; break;
        case 19 : igb_head.v_inc_t = std::stof(search->second); igb_head.bool_inc_t = true; break;
        case 20 : igb_head.v_org_x = std::stof(search->second); igb_head.bool_org_x = true; break;
        case 21 : igb_head.v_org_y = std::stof(search->second); igb_head.bool_org_y = true; break;
        case 22 : igb_head.v_org_z = std::stof(search->second); igb_head.bool_org_z = true; break;
        case 23 : igb_head.v_org_t = std::stof(search->second); igb_head.bool_org_t = true; break;
        case 24 : strncpy(igb_head.v_unites_x, search->second.c_str(), 40); igb_head.bool_unites_x = true; break;
        case 25 : strncpy(igb_head.v_unites_y, search->second.c_str(), 40); igb_head.bool_unites_y = true; break;
        case 26 : strncpy(igb_head.v_unites_z, search->second.c_str(), 40); igb_head.bool_unites_z = true; break;
        case 27 : strncpy(igb_head.v_unites_t, search->second.c_str(), 40); igb_head.bool_unites_t = true; break;
        case 28 : igb_head.v_epais = std::stof(search->second); igb_head.bool_epais = true; break;
        case 29 : strncpy(igb_head.v_unites, search->second.c_str(), 40); igb_head.bool_unites = true; break;
        case 30 : igb_head.v_facteur = std::stof(search->second); igb_head.bool_facteur = true; break;
        case 31 : igb_head.v_zero = std::stof(search->second); igb_head.bool_zero = true; break;
        case 32 : strncpy(igb_head.v_struct_desc, search->second.c_str(), 40); igb_head.bool_struct_desc = true; break;
        case 33 : strncpy(igb_head.v_aut_name, search->second.c_str(), 40); igb_head.bool_aut_name = true; break;
        case 34 : strcpy(igb_head.transstr, search->second.c_str()); igb_head.bool_transparent = true; break;
        default:
                  {
                    fprintf(stderr, "ERROR in %s: Unrecognized Keyword \"%s\" found while parsing the IGB header!\n", __func__, search->first.c_str());
                    ret = -1;
                    break;
                  }
      }
    }
  }

  if(igb_head.bool_transparent)
  {
    if( strlen(igb_head.transstr) != 2 * Data_Size[igb_head.v_type])
      fprintf(stderr, "ATTENTION: ignoring invalid transparent value !\n");
    else
    {
      char s[3], *p, *v;
      s[2] = '\0';
      v = (char *)(igb_head.v_transparent = calloc( Data_Size[igb_head.v_type], 1 ));
      for (size_t i=0; i<Data_Size[igb_head.v_type]; i++ ) {
        s[0] = igb_head.transstr[i*2];
        s[1] = igb_head.transstr[i*2+1];
        v[i] = strtol( s, &p, 16 );
      }
    }
  }

  /* --- l'info x y et type est obligatoire --- */
  if ( !igb_head.bool_x || !igb_head.bool_y || !igb_head.bool_type ) {
    fprintf(stderr, "\nERROR in %s: x, y or type not defined\n", __func__) ;
    ret = -1;
    return(ret);
  }
  /* --- calcul des inc et dim --- */
  if ( igb_head.bool_dim_x ) {
    if ( igb_head.bool_inc_x ) {
      float dim_x = igb_head.v_inc_x * (float)(igb_head.v_x) ;
      if ( !almost_equal(dim_x,igb_head.v_dim_x)) {
        fprintf(stderr, "\nATTENTION:\n") ;
        fprintf(stderr,
            "conflit entre x (%ld) * inc_x (%.3g) = %.3g et dim_x (%.3g)\n",
            igb_head.v_x, igb_head.v_inc_x, dim_x, igb_head.v_dim_x) ;
        ret = -1 ;
       }
     }
    else {
      igb_head.v_inc_x = igb_head.v_dim_x / (float)(igb_head.v_x) ;
      igb_head.bool_inc_x = true;
    }
  }
  else
  {
    igb_head.v_dim_x = (float)(igb_head.v_x) * igb_head.v_inc_x ;
    if ( igb_head.bool_inc_x ) igb_head.bool_dim_x = true;
  }

  if ( igb_head.bool_dim_y ) {
    if ( igb_head.bool_inc_y ) {
      float dim_y = igb_head.v_inc_y * (float)(igb_head.v_y) ;
      if ( !almost_equal(dim_y,igb_head.v_dim_y) )
       {
        fprintf(stderr, "\nATTENTION:\n") ;
        fprintf(stderr,
            "conflit entre y (%ld) * inc_y (%.3g) = %.3g et dim_y (%.3g)\n",
            igb_head.v_y, igb_head.v_inc_y, dim_y, igb_head.v_dim_y) ;
        ret = -1 ;
      }
    }
    else {
      igb_head.v_inc_y = igb_head.v_dim_y / (float)(igb_head.v_y) ;
      igb_head.bool_inc_y = true;
    }
  }
  else {
    igb_head.v_dim_y = (float)(igb_head.v_y) * igb_head.v_inc_y ;
    if ( igb_head.bool_inc_y ) igb_head.bool_dim_y = true;
  }

  if ( igb_head.bool_dim_z ) {
    if ( igb_head.bool_inc_z ) {
      float dim_z = igb_head.v_inc_z * (float)(igb_head.v_z) ;
      if ( !almost_equal(dim_z, igb_head.v_dim_z)) {
        fprintf(stderr, "\nATTENTION:\n") ;
        fprintf(stderr,
            "conflit entre z (%ld) * inc_z (%.3g) = %.3g et dim_z (%.3g)\n",
            igb_head.v_z, igb_head.v_inc_z, dim_z, igb_head.v_dim_z) ;
        ret = -1 ;
      }
    }
    else {
      igb_head.v_inc_z = igb_head.v_dim_z / (float)(igb_head.v_z) ;
      igb_head.bool_inc_y = true;
    }
  }
  else {
    igb_head.v_dim_z = (float)(igb_head.v_z) * igb_head.v_inc_z ;
    if ( igb_head.bool_inc_z ) igb_head.bool_dim_z = true;
  }


  if ( igb_head.bool_dim_t ) {
    if ( igb_head.bool_inc_t ) {
      float dim_t = igb_head.v_inc_t * (float)((igb_head.v_t-1)) ;
      if ( !almost_equal(dim_t, igb_head.v_dim_t) ) {
        fprintf(stderr, "\nATTENTION:\n") ;
        fprintf(stderr,
            "conflit entre t (%ld) * inc_t (%.3g) = %.3g et dim_t (%.3g)\n",
            igb_head.v_t, igb_head.v_inc_t, dim_t, igb_head.v_dim_t) ;
        ret = -1;
      }
    }
    else {
      igb_head.v_inc_t = igb_head.v_dim_t / (float)((igb_head.v_t - 1)) ;
      igb_head.bool_inc_t = true;
    }
  }
  else {
    igb_head.v_dim_t = (float)((igb_head.v_t-1)) * igb_head.v_inc_t ;
    if ( igb_head.bool_inc_t )
      igb_head.bool_dim_t = true;
  }

  if ( igb_head.bool_taille ) {
    if (igb_head.v_type!=IGB_STRUCTURE) {
      fprintf(stderr, "\nERREUR taille redefinie pour type autre que structure\n") ;
      ret = -1;
    }
  }
  else {
    if (igb_head.v_type==IGB_STRUCTURE) {
      fprintf(stderr,
          "\nERREUR taille non definie pour type structure\n") ;
      ret = -1;
    }
    else
      igb_head.v_taille = Data_Size[igb_head.v_type];
  }
  return ret;
}

long write_igb_header(igb_header & igb_head)
{
  long ret = 0;
  igb_head.fileptr = fopen(igb_head.filename.c_str(), "wb");
  if(igb_head.fileptr == NULL) igb_treat_file_open_error(igb_head.filename, errno);

  if (igb_head.v_type<IGB_MIN_TYPE || igb_head.v_type>IGB_MAX_TYPE) {
    std::cerr<< "\nHeader_Write: unknown data type: "<< igb_head.v_type;
    return (-1);
  }

  //do again this fancy routine since the amount of points in the header may have changed
  /* --- l'info x y et type est obligatoire --- */
  if ( !igb_head.bool_x || !igb_head.bool_y || !igb_head.bool_type )
  {
    fprintf(stderr, "\nERROR in %s: x, y or type not defined\n", __func__) ;
    ret = -1;
    return(ret);
  }
  /* --- calcul des inc et dim --- */
  if ( igb_head.bool_dim_x ) {
    if ( igb_head.bool_inc_x )
    {
      float dim_x = igb_head.v_inc_x * (float)(igb_head.v_x) ;
      if ( !almost_equal(dim_x, igb_head.v_dim_x) )
      {
        fprintf(stderr, "\nATTENTION:\n") ;
        fprintf(stderr,
            "conflit entre x (%ld) * inc_x (%.3g) = %.3g et dim_x (%.3g)\n",
            igb_head.v_x, igb_head.v_inc_x, dim_x, igb_head.v_dim_x) ;
        ret = -1 ;
      }
    }
    else {
      igb_head.v_inc_x = igb_head.v_dim_x / (float)(igb_head.v_x) ;
      igb_head.bool_inc_x = true;
    }
  }
  else {
    igb_head.v_dim_x = (float)(igb_head.v_x) * igb_head.v_inc_x ;
    if ( igb_head.bool_inc_x ) igb_head.bool_dim_x = true;
  }

  if ( igb_head.bool_dim_y ) {
    if ( igb_head.bool_inc_y ) {
      float dim_y = igb_head.v_inc_y * (float)(igb_head.v_y) ;
      if ( !almost_equal(dim_y, igb_head.v_dim_y) )
      {
        fprintf(stderr, "\nATTENTION:\n") ;
        fprintf(stderr,
            "conflit entre y (%ld) * inc_y (%.3g) = %.3g et dim_y (%.3g)\n",
            igb_head.v_y, igb_head.v_inc_y, dim_y, igb_head.v_dim_y) ;
        ret = -1 ;
      }
    }
    else {
      igb_head.v_inc_y = igb_head.v_dim_y / (float)(igb_head.v_y) ;
      igb_head.bool_inc_y = true;
    }
  }
  else {
    igb_head.v_dim_y = (float)(igb_head.v_y) * igb_head.v_inc_y ;
    if ( igb_head.bool_inc_y ) igb_head.bool_dim_y = true;
  }

  if ( igb_head.bool_dim_z ) {
    if ( igb_head.bool_inc_z ) {
      float dim_z = igb_head.v_inc_z * (float)(igb_head.v_z) ;
      if ( !almost_equal(dim_z, igb_head.v_dim_z) )
      {
        fprintf(stderr, "\nATTENTION:\n") ;
        fprintf(stderr,
            "conflit entre z (%ld) * inc_z (%.3g) = %.3g et dim_z (%.3g)\n",
            igb_head.v_z, igb_head.v_inc_z, dim_z, igb_head.v_dim_z) ;
        ret = -1 ;
      }
    }
    else {
      igb_head.v_inc_z = igb_head.v_dim_z / (float)(igb_head.v_z) ;
      igb_head.bool_inc_y = true;
    }
  }
  else {
    igb_head.v_dim_z = (float)(igb_head.v_z) * igb_head.v_inc_z ;
    if ( igb_head.bool_inc_z ) igb_head.bool_dim_z = true;
  }


  if ( igb_head.bool_dim_t ) {
    if ( igb_head.bool_inc_t ) {
      float dim_t = igb_head.v_inc_t * (float)((igb_head.v_t-1)) ;
      if ( !almost_equal(dim_t, igb_head.v_dim_t) )
      {
        fprintf(stderr, "\nATTENTION:\n") ;
        fprintf(stderr,
            "conflit entre t (%ld) * inc_t (%.3g) = %.3g et dim_t (%.3g)\n",
            igb_head.v_t, igb_head.v_inc_t, dim_t, igb_head.v_dim_t) ;
        ret = -1;
      }
    }
    else {
      igb_head.v_inc_t = igb_head.v_dim_t / (float)((igb_head.v_t - 1)) ;
      igb_head.bool_inc_t = true;
    }
  }
  else {
    igb_head.v_dim_t = (float)((igb_head.v_t-1)) * igb_head.v_inc_t ;
    if ( igb_head.bool_inc_t )
      igb_head.bool_dim_t = true;
  }

  std::string type = Header_type[igb_head.v_type];

  if (igb_head.v_type==IGB_STRUCTURE && igb_head.v_taille<1) {
    std::cerr << "\nHeader_Write: taille invalide:" << igb_head.v_taille << std::endl;
    return (-1);
  }

  if (igb_head.v_trame<MIN_TRAME || igb_head.v_trame>MAX_TRAME) {
    fprintf(stderr, "\nHeader_Write: trame inconnue: %ld\n", igb_head.v_trame);
    return (-1);
  }

  make_igb_consistent(igb_head);

  // we will now only allow writing of big or little endian */
  const char* systeme=(igb_system_endian_code()==IGB_BIG_ENDIAN)?"big_endian":"little_endian";
  char ligne[1024];
  if (igb_head.bool_t) {
    if (igb_head.bool_z) {
      sprintf(ligne, "x:%ld y:%ld z:%ld t:%ld type:%s systeme:%s ",
          igb_head.v_x, igb_head.v_y, igb_head.v_z, igb_head.v_t, type.c_str(), systeme);
    } else {
      sprintf(ligne, "x:%ld y:%ld t:%ld type:%s systeme:%s ",
          igb_head.v_x, igb_head.v_y, igb_head.v_t, type.c_str(), systeme);
    }
  } else {
    if (igb_head.bool_z) {
      sprintf(ligne, "x:%ld y:%ld z:%ld type:%s systeme:%s ",
          igb_head.v_x, igb_head.v_y, igb_head.v_z, type.c_str(), systeme);
    } else {
      sprintf(ligne, "x:%ld y:%ld type:%s systeme:%s ",
          igb_head.v_x, igb_head.v_y, type.c_str(), systeme);
    }
  }

  int n_car = strlen(ligne);

  int n_lignes = 1;
  int n_items = 0;
  int l_item[N_MAX_ITEMS+1];
  char items[N_MAX_ITEMS+1][L_MAX_ITEM];
  /*
     Le mot-clef "taille" n'est ecrit que pour le type STRUCTURE mais il est
     obligatoire pour ce cas.
     */
  if (igb_head.v_type==IGB_STRUCTURE) {
    sprintf(&items[n_items][0], "taille:%ld ", igb_head.v_taille);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_org_x) {
    sprintf(&items[n_items][0], "org_x:%f ", igb_head.v_org_x);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_org_y) {
    sprintf(&items[n_items][0], "org_y:%f ", igb_head.v_org_y);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_org_z) {
    sprintf(&items[n_items][0], "org_z:%f ", igb_head.v_org_z);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_org_t) {
    sprintf(&items[n_items][0], "org_t:%f ", igb_head.v_org_t);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_dim_x) {
    sprintf(&items[n_items][0], "dim_x:%f ", igb_head.v_dim_x);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_inc_x) {
    sprintf(&items[n_items][0], "dim_x:%f ", igb_head.v_inc_x);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_dim_y) {
    sprintf(&items[n_items][0], "dim_y:%f ", igb_head.v_dim_y);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_inc_y) {
    sprintf(&items[n_items][0], "dim_y:%f ", igb_head.v_inc_y);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_dim_z) {
    sprintf(&items[n_items][0], "dim_z:%f ", igb_head.v_dim_z);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_inc_z) {
    sprintf(&items[n_items][0], "dim_z:%f ", igb_head.v_inc_z);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_dim_t) {
    sprintf(&items[n_items][0], "dim_t:%f ", igb_head.v_dim_t);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_inc_t) {
    sprintf(&items[n_items][0], "inc_t:%f ", igb_head.v_inc_t);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_unites_x) {
    sprintf(&items[n_items][0], "unites_x:%.40s ", igb_head.v_unites_x);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_unites_y) {
    sprintf(&items[n_items][0], "unites_y:%.40s ", igb_head.v_unites_y);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_unites_z) {
    sprintf(&items[n_items][0], "unites_z:%.40s ", igb_head.v_unites_z);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_unites_t) {
    sprintf(&items[n_items][0], "unites_t:%.40s ", igb_head.v_unites_t);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_num) {
    sprintf(&items[n_items][0], "num:%ld ", igb_head.v_num);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_bin) {
    sprintf(&items[n_items][0], "bin:%ld ", igb_head.v_bin);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_trame) {
    switch (igb_head.v_trame) {
      case CC8:
        sprintf(&items[n_items][0], "trame:c8 ");
        break;
      case CC4:
        sprintf(&items[n_items][0], "trame:c4 ");
        break;
      case HEX:
        sprintf(&items[n_items][0], "trame:hex ");
        break;
      case HEXEDGES:
        sprintf(&items[n_items][0], "trame:hexedges ");
        break;
      case HEXBRIDGES:
        sprintf(&items[n_items][0], "trame:hexbridges ");
        break;
      case HEXLINES:
        sprintf(&items[n_items][0], "trame:hexlines ");
        break;
      case HEX2:
        sprintf(&items[n_items][0], "trame:hex2 ");
        break;
    }
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_lut) {
    sprintf(&items[n_items][0], "lut:%ld ", igb_head.v_lut);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_comp) {
    sprintf(&items[n_items][0], "comp:%ld ", igb_head.v_comp);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_epais) {
    sprintf(&items[n_items][0], "epais:%f ", igb_head.v_epais);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_unites) {
    sprintf(&items[n_items][0], "unites:%.40s ", igb_head.v_unites);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_facteur) {
    sprintf(&items[n_items][0], "facteur:%f ", igb_head.v_facteur);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_zero) {
    sprintf(&items[n_items][0], "zero:%f ", igb_head.v_zero);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if ( igb_head.v_transparent != NULL ) {
    char *p=(char *)igb_head.v_transparent, value[MAXL];
    for (size_t a=0; a<Data_Size[igb_head.v_type]; a++ )
      sprintf( value+a*2, "%x", *(p++) );
    value[2*Data_Size[igb_head.v_type]] = '\0';
    sprintf(&items[n_items][0], "transparent:%s ", value );
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_struct_desc) {
    sprintf(&items[n_items][0], "struct:%.40s ", igb_head.v_struct_desc);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }
  if (igb_head.bool_aut_name) {
    sprintf(&items[n_items][0], "aut:%.40s ", igb_head.v_aut_name);
    l_item[n_items] = strlen(&items[n_items][0]);
    n_items++;
  }

  int n_car_total = 0;
  /* Ecrit tous les items, sauf les commentaires */
  for (int i=0;i<n_items;i++) {
    if (n_car+l_item[i]<71) {	/*  Ajoute a la ligne courante s'il reste de la place */
      strcat(ligne, &items[i][0]);
      n_car += l_item[i];
    } else {		/*  Sinon, ecrit cette ligne et commence-en une autre */
      ligne[n_car++] = '\r';
      ligne[n_car++] = '\n';
      ligne[n_car]   = '\000';
      n_car_total += n_car;
      if (fputs(ligne, igb_head.fileptr)==-1)
      {
        fprintf(stderr, "\nHeader_Write: Erreur a l'ecriture \n");
        perror("\n *** ");
        fprintf(stderr,  "\n");
        return (-1);
      }
      strcpy(ligne, &items[i][0]);
      n_car = l_item[i];
      n_lignes++;
    }
  }
  /*
     Termine la derniere ligne
     */
  ligne[n_car++] = '\r';
  ligne[n_car++] = '\n';
  ligne[n_car]   = '\000';
  n_car_total += n_car;
  if (fputs(ligne, igb_head.fileptr)==-1)
  {
    fprintf(stderr, "\nHeader_Write: Erreur a l'ecriture \n");
    perror("\n *** ");
    fprintf(stderr,  "\n");
    return (-1);
  }
  n_lignes++;

  /*
     Determine le nombre de caracteres et de lignes supplementaires
     necessaires
     */
  int n_blocs   = 1 + (n_car_total-1)/1024;
  int n_car_sup = n_blocs*1024 - n_car_total;
  int n_lig_sup;
  if (n_car_sup>0) {
    n_lig_sup = 1 + (n_car_sup-1)/72;
  } else {
    n_lig_sup = 0;
  }
  int n_car_dl  = 1 + (n_car_sup-1)%72;

  /*
     Complete l'entete a un multiple de 1024 caracteres
     */
  for (int i=0;i<70;i++) ligne[i] = ' ';
  ligne[70] = '\r';
  ligne[71] = '\n';
  ligne[72] = '\000';
  for (int i=0;i<n_lig_sup-1;i++) {
    if (fputs(ligne, igb_head.fileptr)==-1)
    {
      fprintf(stderr, "\nHeader_Write: Erreur a l'ecriture \n");
      perror("\n *** ");
      fprintf(stderr,  "\n");
      return(-1);
    }
  }

  /*
     La derniere ligne se termine par un saut de page (FF)
     */
  for (int i=0;i<n_car_dl-2;i++) ligne[i] = ' ';
  if (n_car_dl>2) ligne[n_car_dl-3] = '\r';
  if (n_car_dl>1) ligne[n_car_dl-2] = '\n';
  ligne[n_car_dl-1] = FF;
  ligne[n_car_dl] = '\000';
  if (fputs( ligne , igb_head.fileptr)==-1)
  {
    fprintf(stderr, "\nHeader_Write: Erreur a l'ecriture \n");
    perror("\n *** ");
    fprintf(stderr,  "\n");
    return (-1);
  }

  if (n_car_total>1024)
  {
    fprintf(stderr,
          "\nHeader_Write ATTENTION: etiquette de grandeur non-standard \n");
  }
  return ret;
}

bool igb_system_is_big_endian()
{
  int i = 1;
  char *p = (char *)&i;

  if (p[0] == 1) return false;
  else           return true;
}

void set_igb_header_dim(int x, int y, int z, int t, igb_header & igb)
{
  igb.v_x = x; igb.bool_x = true;
  igb.v_y = y; igb.bool_y = true;
  igb.v_z = z; igb.bool_z = true;
  igb.v_t = t; igb.bool_t = true;
}

void set_igb_header_units(const std::string space_unit,
                   const std::string time_unit,
                   const std::string val_unit,
                   igb_header & igb)
{
  sprintf(igb.v_unites, "%.40s", val_unit.c_str()); igb.bool_unites = true;

  if(igb.bool_x)
    sprintf(igb.v_unites_x, "%.40s", space_unit.c_str()), igb.bool_unites_x = true;
  if(igb.bool_y)
    sprintf(igb.v_unites_y, "%.40s", space_unit.c_str()), igb.bool_unites_y = true;
  if(igb.bool_z)
    sprintf(igb.v_unites_z, "%.40s", space_unit.c_str()), igb.bool_unites_z = true;
  if(igb.bool_t)
    sprintf(igb.v_unites_t, "%.40s", time_unit.c_str()), igb.bool_unites_t = true;
}
#endif


