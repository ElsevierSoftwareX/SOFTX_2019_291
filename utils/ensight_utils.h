#ifndef ENSIGHT_UTILS_H
#define ENSIGHT_UTILS_H

struct ens_typeinfo {
  mt_vector<elem_t> types;
  mt_vector<mt_int> eidx_cnt;
  mt_vector<mt_int> eidx_dsp;
  mt_vector<mt_int> eidx;
};

inline const char*
elemtype_to_ensight_string(const elem_t type)
{
  const char* ret = NULL;
  switch(type) {
    case Line:    ret = "bar2";     break;
    case Tri:     ret = "tria3";    break;
    case Quad:    ret = "quad4";    break;
    case Tetra:   ret = "tetra4";   break;
    case Pyramid: ret = "pyramid5"; break;
    case Prism:   ret = "penta6";   break;
    case Hexa:    ret = "hexa8";    break;
  }

  return ret;
}

void get_ensight_type_info(const mt_meshdata & mesh, ens_typeinfo & ti);

void write_ensight_mesh(const mt_meshdata & mesh, const ens_typeinfo & ti,
                        bool binary, std::string filename);

void write_ensight_fibers(const mt_vector<mt_real> & lon, const ens_typeinfo & ti,
                          const size_t nelem, const int numfib,
                          const std::string basename);

void write_ensight_meshdata(const mt_meshdata & mesh, const ens_typeinfo & ti,
                            const std::string basename);

void write_ensight_case_header(const std::string casefile,
                               const std::string geometry,
                               bool have_sheets);

template<typename T>
void write_ensight_nodedata(const T* data,
                            const size_t npts,
                            const short dpn,
                            const std::string desc,
                            const std::string filename)
{
  FILE* fd = fopen(filename.c_str(), MT_FOPEN_WRITE);
  if(fd == NULL) treat_file_open_error(filename, errno);

  const int linelen = 80;
  int       intbuff = 1;
  mt_vector<float> fbuff;

  mt_vector<char> lbuff(linelen, ' ');
  char* linebuff = lbuff.data();
  linebuff[linelen-1] = '\n';

  snprintf(linebuff, linelen-1, "%s", desc.c_str()); fwrite(linebuff, 1, linelen, fd);
  sprintf(linebuff, "part                        "); fwrite(linebuff, 1, linelen, fd);
  intbuff = 1; fwrite(&intbuff, sizeof(int), 1, fd);
  sprintf(linebuff, "coordinates                 "); fwrite(linebuff, 1, linelen, fd);

  fbuff.resize(npts);

  for(short j=0; j<dpn; j++) {
    for(size_t i=0; i<npts; i++)
      fbuff[i] = data[i*dpn+j];
    fwrite(fbuff.data(), sizeof(float), npts, fd);
  }

  fclose(fd);
}

template<typename T>
void write_ensight_elemdata(const ens_typeinfo & ti,
                            const T* data,
                            const short dpn,
                            const std::string desc,
                            const std::string filename)
{
  FILE* fd = fopen(filename.c_str(), MT_FOPEN_WRITE);
  if(fd == NULL) treat_file_open_error(filename, errno);

  const int linelen = 80;
  int       intbuff = 1;
  mt_vector<float> fbuff;

  mt_vector<char> lbuff(linelen, ' ');
  char* linebuff = lbuff.data();
  linebuff[linelen-1] = '\n';

  snprintf(linebuff, linelen-1, "%s", desc.c_str()); fwrite(linebuff, 1, linelen, fd);
  sprintf(linebuff, "part                        "); fwrite(linebuff, 1, linelen, fd);
  intbuff = 1; fwrite(&intbuff, sizeof(int), 1, fd);

  // now the elements by type
  for(size_t tidx = 0; tidx < ti.types.size(); tidx++) {
    elem_t curtype  = ti.types   [tidx];
    long   cur_nele = ti.eidx_cnt[tidx];
    long   cur_dsp  = ti.eidx_dsp[tidx];

    // element type
    sprintf(linebuff, "%s  ", elemtype_to_ensight_string(curtype));
    fwrite(linebuff, 1, linelen, fd);
    fbuff.resize(cur_nele);

    // per element-type values, sorted by component
    for(short j=0; j<dpn; j++) {
      for(long i=0; i<cur_nele; i++) {
        mt_int eidx = ti.eidx[cur_dsp + i];
        fbuff[i] = data[eidx*dpn+j];
      }
      fwrite(fbuff.data(), sizeof(float), cur_nele, fd);
    }
  }

  fclose(fd);
}


#endif
