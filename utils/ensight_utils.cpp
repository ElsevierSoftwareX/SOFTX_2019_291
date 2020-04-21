#include "mt_utils.h"

void write_ensight_mesh(const mt_meshdata & mesh, const ens_typeinfo & ti,
                        bool binary, std::string filename)
{
  FILE* fd = fopen(filename.c_str(), MT_FOPEN_WRITE);
  if(fd == NULL) treat_file_open_error(filename, errno);

  long npts = mesh.xyz.size() / 3;

  mt_vector<mt_int> e2n_dsp;
  bucket_sort_offset(mesh.e2n_cnt, e2n_dsp);

  if(binary == false) {
    PROGRESS<int> prg((2 + 2*ti.types.size()), "Writing mesh in ensight txt format");

    fprintf(fd, "Ensight ASCII geometry file\n");
    fprintf(fd, "Written by meshtool\n");
    // we provide both element and node IDs
    fprintf(fd, "node id given\n");
    fprintf(fd, "element id given\n");
    // we currently write only one part
    int part_id = 1;
    fprintf(fd, "part\n%10d\n", part_id);
    fprintf(fd, "======================================\n");

    // first the coordinates
    fprintf(fd, "coordinates\n%10ld\n", npts);
    // index part
    for(long i=0; i<npts; i++)
      fprintf(fd, "%10ld\n", i+1);
    prg.next();
    // coordinates by component
    for(long i=0; i<npts; i++)
      fprintf(fd, "%12.5e\n", mesh.xyz[i*3+0]);
    for(long i=0; i<npts; i++)
      fprintf(fd, "%12.5e\n", mesh.xyz[i*3+1]);
    for(long i=0; i<npts; i++)
      fprintf(fd, "%12.5e\n", mesh.xyz[i*3+2]);
    prg.next();

    // now the elements by type
    for(size_t tidx = 0; tidx < ti.types.size(); tidx++) {
      elem_t curtype  = ti.types   [tidx];
      long   cur_nele = ti.eidx_cnt[tidx];
      long   cur_dsp  = ti.eidx_dsp[tidx];

      fprintf(fd, "%s\n", elemtype_to_ensight_string(curtype));
      fprintf(fd, "%10ld\n", cur_nele);

      // write element indices
      for(long i=0; i<cur_nele; i++) {
        mt_int eidx = ti.eidx[cur_dsp + i];
        fprintf(fd, "%10ld\n", eidx+1);
      }
      prg.next();

      // write element connectivities
      for(long i=0; i<cur_nele; i++) {
        mt_int eidx = ti.eidx[cur_dsp + i];
        const mt_int* con = mesh.e2n_con.data() + e2n_dsp[eidx];

        for(mt_int j=0; j<mesh.e2n_cnt[eidx]; j++)
          fprintf(fd, "%10ld ", con[j]+1);
        fprintf(fd, "\n");
      }

      prg.next();
    }
    prg.finish();
  }
  else {
    PROGRESS<int> prg((2 + 2*ti.types.size()), "Writing mesh in ensight binary format");
    const int linelen = 80;
    int intbuff = 1;

    mt_vector<char> lbuff(linelen, ' ');
    char* linebuff = lbuff.data();
    linebuff[linelen-1] = '\n';

    mt_vector<int> idxbuff(npts);
    mt_vector<float> fbuff(npts);

    sprintf(linebuff, "C Binary");                     fwrite(linebuff, 1, linelen, fd);
    sprintf(linebuff, "Ensight binary geometry file"); fwrite(linebuff, 1, linelen, fd);
    sprintf(linebuff, "Written by meshtool         "); fwrite(linebuff, 1, linelen, fd);
    sprintf(linebuff, "node id given               "); fwrite(linebuff, 1, linelen, fd);
    sprintf(linebuff, "element id given            "); fwrite(linebuff, 1, linelen, fd);
    sprintf(linebuff, "part                        "); fwrite(linebuff, 1, linelen, fd);
    intbuff = 1;    fwrite(&intbuff, sizeof(int), 1, fd);
    sprintf(linebuff, "============================"); fwrite(linebuff, 1, linelen, fd);
    sprintf(linebuff, "coordinates                 "); fwrite(linebuff, 1, linelen, fd);
    intbuff = npts; fwrite(&intbuff, sizeof(int), 1, fd);

    prg.next();

    // vertex indices
    for(int i=0; i<npts; i++) idxbuff[i] = i+1;
    fwrite(idxbuff.data(), sizeof(int), npts, fd);
    prg.next();
    // coords x components
    for(int i=0; i<npts; i++) fbuff[i] = mesh.xyz[i*3+0];
    fwrite(fbuff.data(), sizeof(float), npts, fd);
    // coords y components
    for(int i=0; i<npts; i++) fbuff[i] = mesh.xyz[i*3+1];
    fwrite(fbuff.data(), sizeof(float), npts, fd);
    // coords z components
    for(int i=0; i<npts; i++) fbuff[i] = mesh.xyz[i*3+2];
    fwrite(fbuff.data(), sizeof(float), npts, fd);
    prg.next();

    // now the elements by type
    for(size_t tidx = 0; tidx < ti.types.size(); tidx++) {
      elem_t curtype  = ti.types   [tidx];
      long   cur_nele = ti.eidx_cnt[tidx];
      long   cur_dsp  = ti.eidx_dsp[tidx];

      // element type
      sprintf(linebuff, "%s  ", elemtype_to_ensight_string(curtype));
      fwrite(linebuff, 1, linelen, fd);
      // number of elements
      intbuff = cur_nele; fwrite(&intbuff, sizeof(int), 1, fd);

      // element indices
      idxbuff.resize(cur_nele);
      for(long i=0; i<cur_nele; i++) {
        mt_int eidx = ti.eidx[cur_dsp + i];
        idxbuff[i] = eidx+1;
      }

      fwrite(idxbuff.data(), sizeof(int), cur_nele, fd);
      prg.next();

      idxbuff.resize(0);
      idxbuff.reserve(mesh.e2n_con.size());

      // write element connectivities
      for(long i=0; i<cur_nele; i++) {
        mt_int eidx = ti.eidx[cur_dsp + i];
        const mt_int* con = mesh.e2n_con.data() + e2n_dsp[eidx];

        for(mt_int j=0; j<mesh.e2n_cnt[eidx]; j++)
          idxbuff.push_back(int(con[j]+1));
      }

      fwrite(idxbuff.data(), sizeof(int), idxbuff.size(), fd);
      prg.next();
    }
    prg.finish();
  }
  fclose(fd);

  // we also write the mesh data
  // we extract the name without extension by using the mt_filename functionality
  mt_filename geofile(filename, "");
  write_ensight_meshdata(mesh, ti, geofile.base);
}

void write_ensight_fibers(const mt_vector<mt_real> & lon, const ens_typeinfo & ti,
                          const size_t nelem, const int numfib,
                          const std::string basename)
{
  check_condition(numfib == 1 || numfib == 2, "number of fibers is 1 or 2", __func__);
  check_condition((nelem * numfib * 3) == lon.size(), "fiber values and number of elems match",
                  __func__);

  mt_vector<mt_real> fib(nelem * 3);

  for(size_t i=0; i<nelem; i++) {
    fib[i*3+0] = lon[i*3*numfib+0];
    fib[i*3+1] = lon[i*3*numfib+1];
    fib[i*3+2] = lon[i*3*numfib+2];
  }

  std::string filename = basename + ".fiber" + ENSIGHT_DATA_EXT;
  std::string pathless_basename = mt_basename(basename);
  std::string desc = pathless_basename + " fiber directions";

  printf("Writing elemdata \"%s\" to: %s\n", desc.c_str(), filename.c_str());
  write_ensight_elemdata(ti, fib.data(), 3, desc, filename);

  if(numfib == 2) {
    for(size_t i=0; i<nelem; i++) {
      fib[i*3+0] = lon[i*3*numfib+3];
      fib[i*3+1] = lon[i*3*numfib+4];
      fib[i*3+2] = lon[i*3*numfib+5];
    }

    filename = basename + ".sheet" + ENSIGHT_DATA_EXT;
    desc     = pathless_basename + " sheet directions";

    printf("Writing elemdata \"%s\" to: %s\n", desc.c_str(), filename.c_str());
    write_ensight_elemdata(ti, fib.data(), 3, desc, filename);
  }
}


void write_ensight_meshdata(const mt_meshdata & mesh, const ens_typeinfo & ti,
                            const std::string basename)
{
  short numfib = mesh.lon.size() == mesh.e2n_cnt.size() * 6 ? 2 : 1;

  write_ensight_fibers(mesh.lon, ti, mesh.e2n_cnt.size(), numfib, basename);

  std::string pathless_basename = mt_basename(basename);
  std::string filename = basename + ".tags" + ENSIGHT_DATA_EXT;
  std::string desc     = pathless_basename + " element tags";

  printf("Writing elemdata \"%s\" to: %s\n", desc.c_str(), filename.c_str());
  write_ensight_elemdata(ti, mesh.etags.data(), 1, desc, filename);
}

void write_ensight_case_header(const std::string casefile,
                               const std::string geometry,
                               bool have_sheets)
{
  FILE* fd = fopen(casefile.c_str(), MT_FOPEN_WRITE);
  if(fd == NULL) treat_file_open_error(casefile, errno);

  printf("Writing ensight case file %s \n", casefile.c_str());

  fprintf(fd, "FORMAT\n");
  fprintf(fd, "type: ensight gold\n");
  fprintf(fd, "GEOMETRY\n");
  fprintf(fd, "model: %s\n", geometry.c_str());
  fprintf(fd, "\nVARIABLE\n\n");

  mt_filename geofile(geometry, "");
  std::string datastr = geofile.base + ".fiber" + ENSIGHT_DATA_EXT;
  fprintf(fd, "vector per element: fiber %s\n", datastr.c_str());

  if(have_sheets) {
    datastr = geofile.base + ".sheet" + ENSIGHT_DATA_EXT;
    fprintf(fd, "vector per element: sheet %s\n", datastr.c_str());
  }

  datastr = geofile.base + ".tags" + ENSIGHT_DATA_EXT;
  fprintf(fd, "scalar per element: tags %s\n", datastr.c_str());

  fclose(fd);
}

void get_ensight_type_info(const mt_meshdata & mesh, ens_typeinfo & ti)
{
  // ensight sorts elements by type therefore we first compute the set of types
  MT_USET<int> types_set;
  for(elem_t t : mesh.etype) types_set.insert(int(t));
  types_set.sort();
  long nele = mesh.e2n_cnt.size();

  ti.types   .resize(types_set.size());
  ti.eidx_cnt.resize(types_set.size());
  ti.eidx    .reserve(nele);
  int widx = 0;

  for(int type_int : types_set) {
    elem_t curtype  = elem_t(type_int);
    ti.types   [widx] = curtype;
    ti.eidx_cnt[widx] = 0;

    // count number of elements of type curtype
    for(long i=0; i<nele; i++)
      if(mesh.etype[i] == curtype) ti.eidx_cnt[widx]++;

    // element indices
    for(long i=0; i<nele; i++)
      if(mesh.etype[i] == curtype) ti.eidx.push_back(i);

    widx++;
  }

  bucket_sort_offset(ti.eidx_cnt, ti.eidx_dsp);
}

