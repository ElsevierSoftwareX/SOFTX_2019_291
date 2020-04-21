/**
* @file collect_mode.cpp
* @brief Merge a mesh with datasets into a collection
* @author Aurel Neic, Elias Karabelas, Matthias Gsell, Andrew Crosier
* @version 
* @date 2019-11-07
*/

#include "mt_modes_base.h"

static const std::string fib_par = "-fib=";
static const std::string dsp_par = "-dsp=";
static const std::string nod_par = "-nod=";
static const std::string ele_par = "-ele=";
static const std::string tstart_par = "-tstart=";
static const std::string tend_par   = "-tend=";
static const std::string inc_par    = "-inc=";
static const std::string trim_par   = "-trim=";
static const std::string extend_par = "-extend=";
static const std::string periodic_par      = "-prd=";
static const std::string periodic_init_par = "-prd_init=";

enum file_format {CARP_TXT, CARP_BIN, IGB, DAT, VTX, VEC, VTU_ASCII, VTU_BIN, ENS_TXT, ENS_BIN, FMT_UNSET};

struct data_handler {
  std::string filename = "";
  bool        exists = false;
  file_format fmt = FMT_UNSET;

  std::vector<float> data_buff;
  std::vector<float> backup_buff;
  igb_header hdr;
  fpos_t     rewind_pos;

  int       tsteps;
  int       tsteps_initial_phase;
  int       tsteps_period;
  int       nodes;
  int       dim;
  float     dt;
  float     Torg;
  size_t    timeslice_sz;
  size_t    read_cnt;
};

struct collect_options
{
  // files
  mt_filename inp_file;
  mt_filename out_file;

  // string options
  std::string out_fmt;
  std::string dsp;
  mt_vector<std::string> fib;
  mt_vector<std::string> vtx;
  mt_vector<std::string> nod;
  mt_vector<std::string> ele;

  // other type options
  mt_real scale;
  bool write_idx;
  bool absolute_dsp;
  bool trim_name;
  bool extend_data;
  int periodic;
  int periodic_init;
  int tstart;
  int tend;
  int inc;
};

void print_collect_help()
{
  fprintf(stderr, "collect: merge a mesh with datasets \n");
  fprintf(stderr, "parameters:\n");
  fprintf(stderr, "%s<path>\t (input) path to basename of the input mesh\n", inp_mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (output) path to basename of the output mesh\n", out_mesh_par.c_str());
  fprintf(stderr, "%s<path>\t (optional) Additional lon file to load.\n", fib_par.c_str());
  fprintf(stderr, "%s<path>\t (optional) .dynpt / .igb displacement file.\n", dsp_par.c_str());
  fprintf(stderr, "%s<path>\t (optional) .vtx vertex file to load.\n", vtx_par.c_str());
  fprintf(stderr, "%s<path>\t (optional) nodal data to include. May be .dat, .vec or .igb.\n", nod_par.c_str());
  fprintf(stderr, "%s<path>\t (optional) element data to include. May be .dat, .vec or .igb.\n", ele_par.c_str());
  fprintf(stderr, "%s<int>\t (optional) trim filenames? 0 = no, 1 = yes, default is 1.\n", trim_par.c_str());
  fprintf(stderr, "%s<int>\t (optional) include node and element indices? 0 = no, 1 = yes, default is 0.\n", idx_par.c_str());
  fprintf(stderr, "%s<int>\t (optional) extend incomplete data? 0 = no, 1 = yes, default is 0.\n", extend_par.c_str());
  fprintf(stderr, "%s<int>\t (optional) extend the output periodically by given periods? 0 = off. default is 0.\n", periodic_par.c_str());
  fprintf(stderr, "%s<int>\t (optional) How many timesteps describe the initial non periodic phase.\n", periodic_init_par.c_str());
  fprintf(stderr, "%s<int>\t (optional) start time.\n", tstart_par.c_str());
  fprintf(stderr, "%s<int>\t (optional) end time.\n", tend_par.c_str());
  fprintf(stderr, "%s<int>\t (optional) time increment.\n", inc_par.c_str());
  fprintf(stderr, "%s<float>\t (optional) Vertex scaling applied to output mesh\n", scale_par.c_str());
  fprintf(stderr, "%s<path>\t (optional) format of the input mesh\n", inp_format_par.c_str());
  fprintf(stderr, "%s<path>\t (optional) format of the output mesh\n", out_format_par.c_str());
  fprintf(stderr, "\nThe supported input formats are:\n%s\n", input_formats.c_str());
  fprintf(stderr, "\nThe supported output formats are:\nvtu_txt, vtu_bin, ens_txt, ens_bin\n");
  fprintf(stderr, "\n");
}


int collect_parse_options(int argc, char** argv, collect_options & opts)
{
  if(argc < 3) {
    print_collect_help();
    return 1;
  }

  std::string inp_msh_base;
  std::string out_msh_base;
  std::string inp_format;
  std::string strbuff;

  // to reduce complexity in the main mode code, we want to convert the parameter strings
  // already here. Also, we set default values here.
  opts.scale = 1.0;
  opts.write_idx = false;
  opts.absolute_dsp  = false;
  opts.trim_name     = true;
  opts.extend_data   = false;
  opts.periodic      = 0;
  opts.periodic_init = 0;
  opts.tstart        = 0;
  opts.tend          = -1;
  opts.inc           = 1;

  for(int i=2; i<argc; i++) {
    std::string param = argv[i];
    bool match = false;

    if(!match) match = parse_param(param, inp_mesh_par, inp_msh_base);
    if(!match) match = parse_param(param, out_mesh_par, out_msh_base);
    if(!match) match = parse_param(param, dsp_par, opts.dsp);
    if(!match) match = parse_param(param, inp_format_par, inp_format);
    if(!match) match = parse_param(param, out_format_par, opts.out_fmt);

    // scale
    if(!match && (match = parse_param(param, scale_par, strbuff)))
      opts.scale = atof(strbuff.c_str());
    // trim names
    if(!match && (match = parse_param(param, trim_par, strbuff)))
      opts.trim_name = atoi(strbuff.c_str()) == 1;
    // write indices
    if(!match && (match = parse_param(param, idx_par, strbuff)))
      opts.write_idx = atoi(strbuff.c_str()) == 1;
    // extend
    if(!match && (match = parse_param(param, extend_par, strbuff)))
      opts.extend_data = atoi(strbuff.c_str()) == 1;
    // periodic
    if(!match && (match = parse_param(param, periodic_par, strbuff)))
      opts.periodic = atoi(strbuff.c_str());
    // periodic init
    if(!match && (match = parse_param(param, periodic_init_par, strbuff)))
      opts.periodic_init = atoi(strbuff.c_str());
    // tstart
    if(!match && (match = parse_param(param, tstart_par, strbuff)))
      opts.tstart = atoi(strbuff.c_str());
    // tend
    if(!match && (match = parse_param(param, tend_par, strbuff)))
      opts.tend = atoi(strbuff.c_str());
    // inc
    if(!match && (match = parse_param(param, inc_par, strbuff)))
      opts.inc = atoi(strbuff.c_str());

    if(!match && (match = parse_param(param, fib_par, strbuff)))
      opts.fib.push_back(strbuff);
    if(!match && (match = parse_param(param, vtx_par, strbuff)))
      opts.vtx.push_back(strbuff);
    if(!match && (match = parse_param(param, nod_par, strbuff)))
      opts.nod.push_back(strbuff);
    if(!match && (match = parse_param(param, ele_par, strbuff)))
      opts.ele.push_back(strbuff);

    if(!match) {
      std::cerr << "Error: Cannot parse parameter " << param << std::endl;
      return 2;
    }
  }

  opts.inp_file.assign(inp_msh_base, inp_format);
  opts.out_file.assign(out_msh_base, "");
  opts.out_file.format = "vtk_bin";

  if(! (opts.inp_file.isSet() && opts.out_file.isSet()) ) {
    std::cerr << "Mesh convert error: Insufficient parameters provided." << std::endl;
    print_collect_help();
    return 3;
  }

  return 0;
}

void filedesc_check_extension_data(data_handler & file)
{
  // Get extension
  const char *dot = strrchr(file.filename.c_str(), '.');

  if (dot && !strcmp(dot, ".igb"))      file.fmt = IGB;
  else if (dot && !strcmp(dot, ".dat")) file.fmt = DAT;
  else if (dot && !strcmp(dot, ".vec")) file.fmt = VEC;
  else {
    fprintf(stderr, "ERROR: Unsupported data file type, file '%s'.\n", file.filename.c_str());
    exit(EXIT_FAILURE);
  }
}

void filedesc_check_exists(data_handler & file, bool required)
{
  file.exists = file_exists(file.filename.c_str());
  if (!file.exists && required) {
    fprintf(stderr, "ERROR: File '%s' does not exist\n", file.filename.c_str());
    exit(EXIT_FAILURE);
  }
}

#define HEADER_SIZE  1024             //!< header size of binary mesh files
#define BUFSIZE      5000

void parse_data(data_handler & ndat, const int nvals, const bool extend,
                const bool extend_periodic, const int periodic_initialsteps, const int periods)
{
  mt_vector<float>  fbuff;
  mt_vector<mt_int> ibuff;

  switch (ndat.fmt) {
    case DAT: {
      ndat.dim = 1;
      ndat.tsteps = 1;
      read_vector_ascii(fbuff, ndat.filename);
      ndat.nodes = fbuff.size();
      ndat.data_buff.assign(fbuff.begin(), fbuff.end());
      break;
    }

    case VEC: {
      ndat.dim = 3;
      ndat.tsteps = 1;
      read_vector_ascii(fbuff, ndat.filename);
      ndat.nodes = fbuff.size() / 3;
      ndat.data_buff.assign(fbuff.begin(), fbuff.end());
      break;
    }

    case VTX: {
      ndat.dim = 1;
      ndat.tsteps = 1;
      read_vtx(ibuff, ndat.filename);
      ndat.data_buff.assign(nvals, 0.0);
      ndat.nodes = nvals;
      for(const int v : ibuff) ndat.data_buff[v] = 1.0;
      break;
    }

    case IGB: {
      //Parse IGB Header
      igb_header & igbhead = ndat.hdr;
      init_igb_header(ndat.filename, igbhead);
      read_igb_header(igbhead);

      int num_pts = igbhead.v_x, tsteps = igbhead.v_t;
      int dim = Num_Components[igbhead.v_type];
      float Torg = igbhead.v_org_t, dt = igbhead.v_inc_t;

      //Figure out how many timesteps we really write when doing periodic continuation
      const int tsteps_period = tsteps - periodic_initialsteps - 1;
      const int total_tsteps  = extend_periodic ? periodic_initialsteps + periods * tsteps_period + 1 : tsteps;

      // Allocate data array
      ndat.nodes                = num_pts;
      ndat.dim                  = dim;
      ndat.tsteps               = total_tsteps;
      ndat.tsteps_period        = tsteps_period;
      ndat.tsteps_initial_phase = periodic_initialsteps;
      ndat.dt                   = dt;
      ndat.Torg                 = Torg;
      ndat.data_buff  .assign(num_pts * dim, 0.);
      ndat.backup_buff.assign(num_pts * dim, 0.);
      ndat.timeslice_sz = num_pts * dim * sizeof(float);
      ndat.read_cnt = 0;
      break;
    }

    default:
      fprintf(stderr, "ERROR: Unsupported file format, file '%s'.\n", ndat.filename.c_str());
      exit(EXIT_FAILURE);
  }
}

void read_igb(data_handler & ndat,
              const bool extend,
              const bool extend_periodic)
{
  //backup last full data
  if(ndat.data_buff.size())
    ndat.backup_buff.assign(ndat.data_buff.begin(), ndat.data_buff.end());

  //return backup buffer in case we have only one time step
  if(ndat.tsteps == 1 && ndat.read_cnt > 0) {
    ndat.data_buff.assign(ndat.backup_buff.begin(), ndat.backup_buff.end());
    return;
  }

  // set filepos to be able to extend data periodically
  if(extend_periodic && ndat.read_cnt == size_t(ndat.tsteps_initial_phase)) {
    fgetpos(ndat.hdr.fileptr, &ndat.rewind_pos);
  }

  // Load data
  size_t nvals = ndat.nodes * ndat.dim;
  size_t count = read_igb_slice(ndat.data_buff, ndat.hdr);

  if(count < nvals) {
    fprintf(stderr, "Warning: Number of values read does not match igb header, file '%s'\n",
            ndat.hdr.filename.c_str());
    fprintf(stderr, "         read:     %zd\n", count);
    fprintf(stderr, "         expected: %zd\n", nvals);

    if(extend) {
      fprintf(stdout, "Warning: Truncated data in file '%s' at step %zd! Extend with last full read!",
              ndat.hdr.filename.c_str(), ndat.read_cnt);
      ndat.data_buff.assign(ndat.backup_buff.begin(), ndat.backup_buff.end());
    }
    else {
      fprintf(stdout, "Warning: Truncated data in file '%s' at step %zd! Extend with zero!",
              ndat.hdr.filename.c_str(), ndat.read_cnt);
      ndat.data_buff.assign(nvals, 0.);
    }
  }

  ndat.read_cnt++;

  // rewind file
  if (extend_periodic &&
      ndat.read_cnt == size_t(ndat.tsteps_initial_phase + ndat.tsteps_period + 1))
  {
    fsetpos(ndat.hdr.fileptr, &ndat.rewind_pos);
    ndat.read_cnt = ndat.tsteps_initial_phase;
  }
}

void read_dynpt(data_handler & ndat,
                const bool extend,
                const bool extend_periodic)
{

  // Read data from a .dynpt file
  // This function is not much more than a glorified alias of read_igb. All it
  // does extra is to check that the read file contains vector data.

  read_igb(ndat, extend, extend_periodic);
  if (ndat.dim != 3) {
    fprintf(stderr, "ERROR: .dynpt file must be of dimension 3\n");
    exit(EXIT_FAILURE);
  }
}

void read_data(data_handler & ndat,
               const int nvals,
               const bool extend,
               const bool extend_periodic)
{
  // Read the data from a suitable file as a NodeDataArray
  //   The nvals parameter is only needed for formats like .dat and .vec that do
  //   not have a header containing this information
  switch (ndat.fmt) {
    case IGB:
      read_igb(ndat, extend, extend_periodic);
      break;
    case DAT:
    case VEC:
    case VTX:
      break;
    default:
      fprintf(stderr, "ERROR: Unsupported file format, file '%s'.\n", ndat.filename.c_str());
      exit(EXIT_FAILURE);
  }
}

void write_VTU_indices (const size_t ninds, const char* name, bool binary, FILE* output_file)
{
  mt_vector<int> inds(ninds);
  for(size_t k=0; k < inds.size(); k++) inds[k] = k;

  bool is_float = false;
  write_VTU_dataset(inds.data(), inds.size(), 1, name, is_float, binary, output_file);
}

class output_data_container
{
  public:
  mt_vector<data_handler>     * nodedata;
  mt_vector<data_handler>     * elemdata;
  data_handler * dynpoint;
  mt_vector<mt_vector<mt_real>> * fibres;
  mt_vector<data_handler> * fibre_files;
  mt_vector<int>          * fibre_axes;

  output_data_container(mt_vector<data_handler>     & nod,
                        mt_vector<data_handler>     & elem,
                        data_handler & dyn,
                        mt_vector<mt_vector<mt_real>> & fib,
                        mt_vector<data_handler> & fib_files,
                        mt_vector<int>         & fib_axes) :
                        nodedata(&nod),
                        elemdata(&elem),
                        dynpoint(&dyn),
                        fibres(&fib),
                        fibre_files(&fib_files),
                        fibre_axes(&fib_axes)
  {}
};


class basic_writer
{
  public:
  double dur_write = 0;
  double dur_write_lon = 0;
  double dur_write_dyn = 0;
  double dur_write_nodedat = 0;
  double dur_write_elemdat = 0;

  virtual void output_data(double time) = 0;
  virtual void output_collection() = 0;
  virtual ~basic_writer() {}
};

class vtu_writer : public basic_writer
{
  public:
  const std::string basename;
  const mt_meshdata & imesh;
  const output_data_container & out_data;

  const bool write_idx;
  const bool binary_fmt;
  const bool trim_name;
  const int  tsteps;

  int output_idx;
  mt_vector<double>      times;
  mt_vector<std::string> filenames;

  vtu_writer(const std::string inp_basename,
             const mt_meshdata & inp_mesh,
             const output_data_container & inp_out_data,
             const bool inp_write_idx,
             const bool inp_binary_fmt,
             const bool inp_trim_name,
             const int inp_tsteps) :
             basename(inp_basename),
             imesh(inp_mesh),
             out_data(inp_out_data),
             write_idx(inp_write_idx),
             binary_fmt(inp_binary_fmt),
             trim_name(inp_trim_name),
             tsteps(inp_tsteps),
             output_idx(0)
  {}

  void output_data(double time)
  {
    char oname[BUFSIZE];
    if (tsteps == 1) {
      sprintf(oname, "%s.vtu", basename.c_str());
    } else {
      sprintf(oname, "%s_%04d.vtu", basename.c_str(), output_idx);
    }

    times    .push_back(time);
    filenames.push_back(mt_basename(oname));

    FILE* out = fopen(oname,MT_FOPEN_WRITE);
    if (out == NULL) treat_file_open_error(oname, errno);

    struct timeval t0, t_end;
    size_t nelems = imesh.e2n_cnt.size();
    size_t npts   = imesh.xyz.size() / 3;

    mt_vector<data_handler>     & nodedata       = *out_data.nodedata;
    mt_vector<data_handler>     & elemdata       = *out_data.elemdata;
    data_handler                & dynpoint       = *out_data.dynpoint;
    mt_vector<mt_vector<mt_real>> & fibres      = *out_data.fibres;
    mt_vector<data_handler>       & fibre_files = *out_data.fibre_files;
    mt_vector<int>                & fibre_axes  = *out_data.fibre_axes;

    // Write the main mesh data ======================================================
    gettimeofday(&t0, NULL);
    write_VTU_mesh(imesh, binary_fmt, false, out);
    gettimeofday(&t_end,NULL);
    dur_write += timediff_sec(t0, t_end);

    // Write cell data ===============================================================
    fprintf(out, "      <CellData>\n");

    // write elem indices
    if (write_idx)
      write_VTU_indices(nelems, "elemIndices", binary_fmt, out);

    // write elem IDs (e.g. material)
    // TODO: make ID name adjustable
    const char* IDname = "elemTags";
    write_VTU_elemtags(imesh.etags, IDname, binary_fmt, out);

    gettimeofday(&t0,NULL);
    int num_axes = imesh.lon.size() == imesh.e2n_cnt.size() * 6 ? 2 : 1;
    write_VTU_londata(imesh.lon, nelems, num_axes, "", binary_fmt, out);

    for (size_t j=0; j<fibre_files.size(); j++) {
      if (fibre_files[j].exists) {
        std::string prefix = mt_basename(fibre_files[j].filename) + "-";
        write_VTU_londata(fibres[j], nelems, fibre_axes[j], prefix.c_str(), binary_fmt, out);
      }
    }
    gettimeofday(&t_end,NULL);
    dur_write_lon += timediff_sec(t0, t_end);

    std::string dataname;
    //write element data
    for(size_t j=0; j < elemdata.size(); j++) {
      dataname = elemdata[j].filename;

      if(trim_name) {
        dataname = mt_basename(dataname);
        remove_extension(dataname);
        dataname += "_" + std::to_string(j);
      }

      gettimeofday(&t0,NULL);
      write_VTU_dataset(elemdata[j].data_buff.data(), nelems, elemdata[j].dim,
                        dataname.c_str(), true, binary_fmt, out);
      gettimeofday(&t_end,NULL);
      dur_write_elemdat += timediff_sec(t0, t_end);
    }
    fprintf(out, "      </CellData>\n");

    // Write point data ==============================================================
    bool points_started = false;
    if (write_idx) {
      fprintf(out, "      <PointData>\n");
      points_started = true;

      write_VTU_indices(npts, "nodeIndices", binary_fmt, out);

      if(nodedata.size() == 0 && dynpoint.hdr.fileptr != NULL)
        fprintf(out, "      </PointData>\n");
    }

    if (dynpoint.exists) {
      if (!points_started) {
        fprintf(out, "      <PointData>\n");
        points_started = true;
      }

      gettimeofday(&t0,NULL);
      write_VTU_dataset(dynpoint.data_buff.data(), npts, 3, "displacement", true, binary_fmt, out);
      gettimeofday(&t_end,NULL);
      dur_write_dyn += timediff_sec(t0, t_end);

      if(nodedata.size() == 0)
        fprintf(out, "      </PointData>\n");
    }

    for(size_t j=0; j < nodedata.size(); j++) {
      if (!points_started) {
        fprintf(out, "      <PointData>\n");
        points_started = true;
      }

      dataname = nodedata[j].filename;
      if(trim_name) {
        dataname = mt_basename(dataname);
        remove_extension(dataname);
        dataname += "_" + std::to_string(j);
      }

      gettimeofday(&t0,NULL);
      write_VTU_dataset(nodedata[j].data_buff.data(), npts, nodedata[j].dim,
                        dataname.c_str(), true, binary_fmt, out);
      gettimeofday(&t_end,NULL);
      dur_write_nodedat += timediff_sec(t0, t_end);
    }

    if(points_started)
      fprintf(out, "      </PointData>\n");

    finish_VTU(out);
    output_idx++;
  }

  void output_collection()
  {
    std::string pvd_file_name = basename + "_collection.pvd";
    FILE *pvd_file = NULL;

    if(tsteps > 1)
    {
#if MT_ENDIANNESS == 0
      const char* endian = "LittleEndian";
#else
      const char* endian = "BigEndian";
#endif

      pvd_file = fopen(pvd_file_name.c_str(), MT_FOPEN_WRITE);

      fprintf(pvd_file, "<?xml version=\"1.0\"?>\n");
      fprintf(pvd_file, "<VTKFile type=\"Collection\" version=\"0.1\"");
      fprintf(pvd_file, " byte_order=\"%s\">\n", endian);
      fprintf(pvd_file, "  <Collection>\n");

      for(size_t i=0; i<size_t(tsteps); i++) {
        fprintf(pvd_file, "    <DataSet timestep=\"%f\"", times[i]);
        fprintf(pvd_file, " group=\"\" part=\"0\" \n");
        fprintf(pvd_file, "             file=\"%s\"/>\n", filenames[i].c_str());
      }

      fprintf(pvd_file, "  </Collection>\n");
      fprintf(pvd_file, "</VTKFile>\n");
      fclose(pvd_file);
    }
  }
};

class ensight_writer : public basic_writer
{
  public:
  const std::string basename;
  const mt_meshdata & imesh;
  const output_data_container & out_data;

  const bool binary_fmt;
  const bool trim_name;
  const int  tsteps, tsteps_integers;

  bool mesh_outputted;
  int output_idx;

  mt_vector<double>      times;
  mt_vector<std::string> filenames;
  mt_vector<std::string> datatype_descr;
  char name_buffer[BUFSIZE];

  std::string generate_desc_string(std::string dataname, int dpn, bool elementdata)
  {
    std::string descr,  datatype;

    switch(dpn) {
      case 1:
        datatype = std::string("scalar");
        break;
      case 3:
        datatype = std::string("vector");
        break;
      case 9:
        datatype = std::string("tensor asym");
        break;
      default:
        fprintf(stderr, "%s error: data dimension %d unsupported yet!\n", __func__, dpn);
        exit(EXIT_FAILURE);
    }

    if(elementdata)
      descr = datatype + " per element: " + dataname;
    else
      descr = datatype + " per node: " + dataname;

    return descr;
  }

  void generate_names(const std::string & inp_filename,
                      const int idx,
                      std::string & out_filename,
                      std::string & description,
                      std::string & case_filename)
  {
    description = inp_filename;
    if(trim_name) {
      description = mt_basename(description);
      remove_extension(description);
    }
    else {
      swap_char(description, '/', '_');
    }

    std::string idx_string = std::to_string(idx);
    mt_vector<char> wildcards;
    wildcards.reserve(36);

    int req_zeros = tsteps_integers - int(idx_string.length());
    if(req_zeros > 0)
      wildcards.assign(req_zeros, '0');

    wildcards.push_back(0);
    sprintf(name_buffer, "%s_%s_%s%s" ENSIGHT_DATA_EXT, basename.c_str(),
            description.c_str(), wildcards.data(), idx_string.c_str());

    out_filename = name_buffer;

    wildcards.assign(tsteps_integers, '*');
    wildcards.push_back(0);
    sprintf(name_buffer, "%s_%s_%s" ENSIGHT_DATA_EXT, basename.c_str(),
            description.c_str(), wildcards.data());
    case_filename = mt_basename(name_buffer);
  }

  void generate_names(const std::string & inp_filename,
                      std::string & out_filename,
                      std::string & description,
                      std::string & case_filename)
  {
    description = inp_filename;
    if(trim_name) {
      description = mt_basename(description);
      remove_extension(description);
    }
    else {
      swap_char(description, '/', '_');
    }

    out_filename = basename + "_" + description + std::string(ENSIGHT_DATA_EXT);
    case_filename = mt_basename(out_filename);
  }


  ensight_writer(const std::string inp_basename,
                 const mt_meshdata & inp_mesh,
                 const output_data_container & inp_out_data,
                 const bool inp_binary_fmt,
                 const bool inp_trim_name,
                 const int inp_tsteps) :
                 basename(inp_basename),
                 imesh(inp_mesh),
                 out_data(inp_out_data),
                 binary_fmt(inp_binary_fmt),
                 trim_name(inp_trim_name),
                 tsteps(inp_tsteps),
                 tsteps_integers(std::to_string(tsteps).length()),
                 mesh_outputted(false),
                 output_idx(0)
  {}

  void output_data(double time)
  {
    struct timeval t0, t_end;
    size_t npts   = imesh.xyz.size() / 3;
    size_t nelems = imesh.e2n_cnt.size();

    times.push_back(time);
    bool populate_filenames = filenames.size() == 0;

    mt_vector<data_handler>     & nodedata      = *out_data.nodedata;
    mt_vector<data_handler>     & elemdata      = *out_data.elemdata;
    data_handler                & dynpoint      = *out_data.dynpoint;
    mt_vector<mt_vector<mt_real>> & fibres      = *out_data.fibres;
    mt_vector<data_handler>       & fibre_files = *out_data.fibre_files;
    mt_vector<int>                & fibre_axes  = *out_data.fibre_axes;

    std::string meshname = basename + ENSIGHT_MESH_EXT;
    ens_typeinfo ti;
    get_ensight_type_info(imesh, ti);

    mt_vector<char> strbuff;
    std::string output_filename, description, case_filename;

    // Write the main mesh data only once
    if(mesh_outputted == false) {
      gettimeofday(&t0, NULL);
      write_ensight_mesh(imesh, ti, binary_fmt, meshname);
      gettimeofday(&t_end,NULL);
      dur_write += timediff_sec(t0, t_end);

      // Write stationary data only once here. ====================
      gettimeofday(&t0,NULL);
      for(size_t j=0; j < elemdata.size(); j++) {
        if(elemdata[j].fmt != IGB) {
          generate_names(elemdata[j].filename, output_filename, description, case_filename);

          write_ensight_elemdata(ti, elemdata[j].data_buff.data(), elemdata[j].dim,
                                 description, output_filename);

          filenames.push_back(case_filename);
          datatype_descr.push_back(generate_desc_string(description, elemdata[j].dim, true));
        }
      }
      gettimeofday(&t_end,NULL);
      dur_write_elemdat += timediff_sec(t0, t_end);

      gettimeofday(&t0,NULL);
      for(size_t j=0; j < nodedata.size(); j++) {
        if(nodedata[j].fmt != IGB) {
          generate_names(nodedata[j].filename, output_filename, description, case_filename);

          write_ensight_nodedata(nodedata[j].data_buff.data(), npts, nodedata[j].dim,
                                 description, output_filename);

          filenames.push_back(case_filename);
          datatype_descr.push_back(generate_desc_string(description, nodedata[j].dim, false));
        }
      }
      gettimeofday(&t_end,NULL);
      dur_write_nodedat += timediff_sec(t0, t_end);

      // write additional fibers =====================================
      gettimeofday(&t0,NULL);
      for (size_t j=0; j<fibre_files.size(); j++) {
        if (fibre_files[j].exists) {
          std::string pathless_fiber_name = mt_basename(fibre_files[j].filename);
          remove_extension(pathless_fiber_name);

          std::string fiber_basename = basename + "." + pathless_fiber_name;
          pathless_fiber_name = mt_basename(fiber_basename);

          write_ensight_fibers(fibres[j], ti, nelems, fibre_axes[j], fiber_basename);

          // we also need to store the files that will be written by write_ensight_fibers
          // for the case file generation
          filenames.push_back(pathless_fiber_name + ".fiber" + ENSIGHT_DATA_EXT);
          datatype_descr.push_back(generate_desc_string(pathless_fiber_name+".fiber", 3, true));
          if(fibre_axes[j] == 2) {
            filenames.push_back(pathless_fiber_name + ".sheet" + ENSIGHT_DATA_EXT);
            datatype_descr.push_back(generate_desc_string(pathless_fiber_name+".sheet", 3, true));
          }
        }
      }
      gettimeofday(&t_end,NULL);
      dur_write_lon += timediff_sec(t0, t_end);

      mesh_outputted = true;
    }

    // Write cell data ===============================================================
    gettimeofday(&t0,NULL);
    for(size_t j=0; j < elemdata.size(); j++) {
      if(elemdata[j].fmt == IGB)
      {
        generate_names(elemdata[j].filename, output_idx, output_filename, description,
                       case_filename);

        write_ensight_elemdata(ti, elemdata[j].data_buff.data(), elemdata[j].dim,
                               description, output_filename);

        if(populate_filenames) {
          filenames.push_back(case_filename);
          datatype_descr.push_back(generate_desc_string(description, elemdata[j].dim, true));
        }
      }
    }
    gettimeofday(&t_end,NULL);
    dur_write_elemdat += timediff_sec(t0, t_end);

    // Write point data ==============================================================
    if (dynpoint.exists) {
      gettimeofday(&t0,NULL);

      sprintf(name_buffer, "%s_displacement_%04d" ENSIGHT_DATA_EXT, basename.c_str(), output_idx);
      write_ensight_nodedata(dynpoint.data_buff.data(), npts, 3, "displacement", name_buffer);

      if(populate_filenames) {
        sprintf(name_buffer, "%s_displacement_0***" ENSIGHT_DATA_EXT, basename.c_str());
        filenames.push_back(mt_basename(name_buffer, strbuff));
        datatype_descr.push_back(generate_desc_string("displacement", 3, false));
      }

      gettimeofday(&t_end,NULL);
      dur_write_dyn += timediff_sec(t0, t_end);
    }

    gettimeofday(&t0,NULL);
    for(size_t j=0; j < nodedata.size(); j++) {
      if(nodedata[j].fmt == IGB) {
        generate_names(nodedata[j].filename, output_idx, output_filename, description,
                       case_filename);

        write_ensight_nodedata(nodedata[j].data_buff.data(), npts, nodedata[j].dim,
                               description, output_filename);

        if(populate_filenames) {
          filenames.push_back(case_filename);
          datatype_descr.push_back(generate_desc_string(description, nodedata[j].dim, false));
        }
      }
    }
    gettimeofday(&t_end,NULL);
    dur_write_nodedat += timediff_sec(t0, t_end);

    output_idx++;
  }

  void output_collection()
  {
    size_t steps = times.size();

    check_condition(filenames.size() == datatype_descr.size(),
                    "filenames and datatype_descr have same size", __func__);

    std::string casefile = basename + ENSIGHT_CASE_EXT;

    write_ensight_case_header(casefile,
                              mt_basename(basename + ENSIGHT_MESH_EXT),
                              imesh.lon.size() == imesh.e2n_cnt.size() * 6);

    FILE* fd = fopen(casefile.c_str(), MT_FOPEN_APPEND);
    if(fd == NULL) treat_file_open_error(casefile, errno);

    for(size_t i=0; i<filenames.size(); i++)
      fprintf(fd, "%s %s\n", datatype_descr[i].c_str(), filenames[i].c_str());

    fprintf(fd, "\nTIME\n\n");
    fprintf(fd, "time set: 1\n");
    fprintf(fd, "number of steps: %d\n", int(steps));
    fprintf(fd, "filename start number: 0\n");
    fprintf(fd, "filename increment: 1\n");
    fprintf(fd, "time values:\n");
    for(size_t i=0; i<times.size(); i++)
      fprintf(fd, "%lf\n", times[i]);
  }
};

void collect_mode(int argc, char** argv)
{
  struct collect_options opts;
  int ret = collect_parse_options(argc, argv, opts);
  if(ret != 0) return;

  struct timeval t0, t_end;  // time measures
  std::string out_basename = opts.out_file.base;

  // additional fibre storage
  size_t n_fibres = opts.fib.size();
  mt_vector<data_handler> fibre_files(n_fibres);  // the file handles to multiple input fibre files
  mt_vector<int>          fibre_axes(n_fibres);   // the number of axes of each input fibre file

  // Were fibre files manually set?
  for (size_t i=0; i<opts.fib.size(); i++) {
    fibre_files[i].filename = opts.fib[i];
    // Check extension
    if (endswith(opts.fib[i], ".lon"))
      fibre_files[i].fmt = CARP_TXT;
    else if (endswith(opts.fib[i], ".blon"))
      fibre_files[i].fmt = CARP_BIN;
    else {
      fprintf(stderr, "ERROR: fibre file '%s' has invalid extension.\n",
          opts.fib[i].c_str());
      exit(EXIT_FAILURE);
    }
  }

  // Was a displacement file set?
  data_handler dynpoint;
  if (opts.dsp.size()) {
    dynpoint.filename = opts.dsp;
    dynpoint.fmt = IGB;
  }

  // Were node data or vertex files set?
  // Vertex files are loaded as node data for visualisation purposes
  size_t n_ndat = opts.nod.size(), n_vtx = opts.vtx.size();
  size_t n_nodedata_files = n_ndat + n_vtx;

  mt_vector<data_handler> nodedata(n_nodedata_files);
  for (size_t i=0; i<n_ndat; i++) {
    nodedata[i].filename = opts.nod[i];
    filedesc_check_extension_data(nodedata[i]);
  }
  for (size_t i=0; i<n_vtx; i++) {
    nodedata[n_ndat+i].filename = opts.vtx[i];
    nodedata[n_ndat+i].fmt      = VTX;
  }

  // Were element data files set?
  size_t n_elemdata_files = opts.ele.size();
  mt_vector<data_handler> elemdata(n_elemdata_files);

  for (size_t i=0; i<n_elemdata_files; i++) {
    elemdata[i].filename = opts.ele[i];
    filedesc_check_extension_data(elemdata[i]);
  }

  // Check files exist
  for (size_t i=0; i<n_fibres; i++)
    filedesc_check_exists(fibre_files[i], true);
  if (opts.dsp.size())
    filedesc_check_exists(dynpoint, true);
  for (size_t i=0; i<n_nodedata_files; i++)
    filedesc_check_exists(nodedata[i], true);
  for (size_t i=0; i<n_elemdata_files; i++)
    filedesc_check_exists(elemdata[i], true);

  // Define data structures
  mt_vector<mt_vector<mt_real>> fibres(n_fibres);

  //=========================================================================
  // READ
  //=========================================================================
  bool extend                       = opts.extend_data;
  bool extend_periodic              = opts.periodic > 0;
  int extend_periodic_periods       = opts.periodic;
  int extend_periodic_initial_steps = opts.periodic_init;

  double dur_mesh = 0, dur_dynpts = 0, dur_lon = 0, dur_nodedat = 0, dur_elemdat = 0;
  mt_meshdata imesh; // the mesh storage struct

  // read mesh
  gettimeofday(&t0,NULL);
  read_mesh_selected(imesh, opts.inp_file.format, opts.inp_file.base);

  // number of points and elements
  const size_t npts = imesh.xyz.size() / 3, nelems = imesh.e2n_cnt.size();

  gettimeofday(&t_end,NULL);
  dur_mesh = timediff_sec(t0, t_end);

  // read additional fibres
  // ------------------------------------------------------------------------
  gettimeofday(&t0,NULL);
  for (size_t i=0; i<n_fibres; i++) {
    switch(fibre_files[i].fmt) {
      case CARP_TXT:
        fibre_axes[i] = readFibers(fibres[i], nelems, fibre_files[i].filename);
        break;

      case CARP_BIN:
        fibre_axes[i] = readFibersBinary(fibres[i], nelems, fibre_files[i].filename);
        break;

      default: assert(0);
    }
  }
  gettimeofday(&t_end,NULL);
  dur_lon = timediff_sec(t0, t_end);

  // parse dynpt header
  // ------------------------------------------------------------------------
  if (dynpoint.exists) {
    gettimeofday(&t0,NULL);
    parse_data(dynpoint, npts, extend, extend_periodic,
               extend_periodic_initial_steps, extend_periodic_periods);
    gettimeofday(&t_end,NULL);
    dur_dynpts = timediff_sec(t0, t_end);
  }

  // read / parse node data
  // ------------------------------------------------------------------------
  gettimeofday(&t0,NULL);
  for (size_t i=0; i<n_nodedata_files; i++)
    parse_data(nodedata[i], npts, extend, extend_periodic,
               extend_periodic_initial_steps, extend_periodic_periods);
  gettimeofday(&t_end,NULL);
  dur_nodedat = timediff_sec(t0, t_end);

  // read / parse element data
  // ------------------------------------------------------------------------
  gettimeofday(&t0,NULL);
  for (size_t i=0; i<n_elemdata_files; i++)
    parse_data(elemdata[i], nelems, extend, extend_periodic,
               extend_periodic_initial_steps, extend_periodic_periods);
  gettimeofday(&t_end,NULL);
  dur_elemdat = timediff_sec(t0, t_end);

  //=========================================================================
  // SCALE AND OFFSET
  //=========================================================================

  if (fabs(opts.scale - 1.0) > 1e-5)
    for(mt_real & v : imesh.xyz) v *= opts.scale;

  if (dynpoint.exists && opts.absolute_dsp)
  {
    //read in first entry in displacement
    fgetpos(dynpoint.hdr.fileptr, &dynpoint.rewind_pos);
    gettimeofday(&t0,NULL);
    read_dynpt(dynpoint, extend, extend_periodic);
    gettimeofday(&t_end,NULL);
    dur_dynpts += timediff_sec(t0, t_end);

    check_size(dynpoint.data_buff.size(), imesh.xyz.size(), "offsetting mesh coords by dynpoints");

    for(size_t i=0; i<imesh.xyz.size(); i++)
      imesh.xyz[i] = dynpoint.data_buff[i] * opts.scale;

    fsetpos(dynpoint.hdr.fileptr, &dynpoint.rewind_pos);
  }

  //=========================================================================
  // VALIDATE
  //=========================================================================

  int tsteps = 1;

  if (dynpoint.exists)
    tsteps = dynpoint.tsteps;

  // Validate node file times are consistent
  for (size_t i=0; i<n_nodedata_files; i++) {
    if (size_t(nodedata[i].nodes) != npts) {
       fprintf(stderr, "Warning: Num nodes in '%s' does not match points file\n",
               nodedata[i].filename.c_str());
    }
    if (nodedata[i].tsteps != tsteps) {
      if (tsteps == 1) {
        // Update the tstep
        tsteps = nodedata[i].tsteps;
      } else if (nodedata[i].tsteps != 1) {
        // File tstep mismatch
        fprintf(stderr, "ERROR: Mismatch between number of time samples in data files\n");
        exit(EXIT_FAILURE);
      }
    }
  }

  // Validate element file times are consistent
  for (size_t i=0; i<n_elemdata_files; i++) {
    if (size_t(elemdata[i].nodes) != nelems) {
       fprintf(stderr, "ERROR: Num elements in '%s' does not match elem file\n",
               elemdata[i].filename.c_str());
       exit(EXIT_FAILURE);
    }
    if (elemdata[i].tsteps != tsteps) {
      if (tsteps == 1) {
        // Update the tstep
        tsteps = elemdata[i].tsteps;
      } else if (elemdata[i].tsteps != 1) {
        // File tstep mismatch
        fprintf(stderr, "ERROR: Mismatch between number of time samples in data files\n");
        exit(EXIT_FAILURE);
      }
    }
  }

  fprintf(stdout, "Data read:\n");
  fprintf(stdout, " - Points:     %ld\n", (long int) npts);
  fprintf(stdout, " - Elements:   %ld\n", (long int) nelems);
  fprintf(stdout, " - Time steps: %d\n\n", tsteps);

  //=========================================================================
  // WRITE
  //=========================================================================

  bool binary_fmt = true;
  bool ensight_fmt = true;
  if(opts.out_fmt.size()) {
    if (opts.out_fmt.compare("vtu_txt") == 0) {
      binary_fmt  = false;
      ensight_fmt = false;
    }
    else if(opts.out_fmt.compare("vtu_bin") == 0) {
      binary_fmt  = true;
      ensight_fmt = false;
    }
    else if(opts.out_fmt.compare("ens_bin") == 0) {
      binary_fmt  = true;
      ensight_fmt = true;
    }
    else if(opts.out_fmt.compare("ens_txt") == 0) {
      binary_fmt  = false;
      ensight_fmt = true;
    }
    else {
      fprintf(stderr, "%s warning: unrecognized output format %s.", __func__, opts.out_fmt.c_str());
    }
  }

  int start = opts.tstart;
  int stop  = opts.tend > -1 ? opts.tend : tsteps;

  if(opts.periodic && (start > 0 || stop < tsteps)) {
    fprintf(stderr, "WARNING: periodic continuation doesn't work with start and stop set manually!\n");
    start = 0;
    stop = tsteps;
  }

  /*
  // TODO: Set up time properly
  float dt = 0.;
  float Torg = 0.;
  if(tsteps > 1) {
    if(n_ndat > 0) {
      dt = nodedata[0].dt;
      Torg = nodedata[0].Torg;
    }
    else if(n_elemdata_files > 0) {
      dt = elemdata[0].dt;
      Torg = elemdata[0].Torg;
    }
  }
  */

  //Reset fileptrs in case of nonzero start point
  if(start > 0) {
    for(size_t j=0; j < n_nodedata_files; j++) {
      if(nodedata[j].fmt == IGB) {
        if(fseek(nodedata[j].hdr.fileptr, start * nodedata[j].timeslice_sz, SEEK_CUR) != 0)
        {
          if(ferror(nodedata[j].hdr.fileptr)) {
            perror("fseek()");
            fprintf(stderr,"fseek() failed in file %s at line # %d\n", __FILE__,__LINE__);
            exit(EXIT_FAILURE);
          }
        }
      }
    }
    for(size_t j=0; j < n_elemdata_files; j++) {
      if(elemdata[j].fmt == IGB) {
        if(fseek(elemdata[j].hdr.fileptr, start * elemdata[j].timeslice_sz, SEEK_CUR) != 0)
        {
          if(ferror(elemdata[j].hdr.fileptr)) {
            perror("fseek()");
            fprintf(stderr,"fseek() failed in file %s at line # %d\n", __FILE__,__LINE__);
            exit(EXIT_FAILURE);
          }
        }
      }
    }
  }

  // pointers to the output data are stored in this container to be easier to pass
  // to output functions
  output_data_container out_data(nodedata, elemdata, dynpoint, fibres, fibre_files, fibre_axes);

  basic_writer* writer = NULL;

  if(ensight_fmt)
    writer = new ensight_writer(out_basename, imesh, out_data, binary_fmt, opts.trim_name, tsteps);
  else
    writer = new vtu_writer(out_basename, imesh, out_data, opts.write_idx, binary_fmt, opts.trim_name, tsteps);

  int cnt = 0;
  for(int i=start; i < stop; i++, cnt++) {
    // All reading happens here =========================================================
    // element data
    for(size_t j=0; j < n_elemdata_files; j++) {
      gettimeofday(&t0,NULL);
      read_data(elemdata[j], nelems, extend, extend_periodic);
      gettimeofday(&t_end,NULL);
      dur_elemdat += timediff_sec(t0, t_end);
    }

    // displacement data
    if (dynpoint.exists) {
      gettimeofday(&t0,NULL);
      read_dynpt(dynpoint, extend, extend_periodic);
      check_size(dynpoint.data_buff.size(), imesh.xyz.size(), "offset dynpoints by mesh coords");

      for(size_t j=0; j<imesh.xyz.size(); j++)
        dynpoint.data_buff[j] = (dynpoint.data_buff[j] * opts.scale) - imesh.xyz[j];

      gettimeofday(&t_end,NULL);
      dur_dynpts += timediff_sec(t0, t_end);
    }

    // nodal data
    for(size_t j=0; j < n_nodedata_files; j++) {
      gettimeofday(&t0,NULL);
      read_data(nodedata[j], npts, extend, extend_periodic);
      gettimeofday(&t_end,NULL);
      dur_nodedat += timediff_sec(t0, t_end);
    }

    // writing =========================================================================
    if(!(i % opts.inc)) {
      fprintf(stdout, "\rWriting timestep: %d of %d", i+1, stop);
      fflush(stdout);
      // Determine output filename

      // write into PVD file
      // TODO: Set up time properly
      // const double time = tsteps > 1 ? static_cast<double>(Torg + i * dt) : 0.;
      const double time = cnt;

      // output
      writer->output_data(time);
    }
  }
  printf("\n");

  // start the collection file
  writer->output_collection();

  // Print some I/O time statistics
  fprintf(stdout, "\nMesh I/O Statistics:\n");

  fprintf(stdout, "\nRead times (%s format mesh):\n", opts.inp_file.format.c_str());
  fprintf(stdout, " %.5f s\n", dur_mesh);
  for (size_t i=0; i<n_fibres; i++) {
    if (fibre_files[i].exists) {
      fprintf(stdout, " - fibres:         %.5f s\n", dur_lon);
      break;
    }
  }
  if (dynpoint.exists)
    fprintf(stdout, " - dynpts:         %.5f s\n", dur_dynpts);
  if (n_nodedata_files != 0)
    fprintf(stdout, " - node/vtx data:  %.5f s\n", dur_nodedat);
  if (n_elemdata_files != 0)
    fprintf(stdout, " - elem data:      %.5f s\n", dur_elemdat);

  fprintf(stdout, "\nWrite times (%s vtk format):\n", opts.out_fmt.c_str());
  fprintf(stdout, " - mesh :      %.5f s\n", writer->dur_write);
  for (size_t i=0; i<n_fibres; i++) {
    if (fibre_files[i].exists) {
      fprintf(stdout, " - fibres:         %.5f s\n", writer->dur_write_lon);
      break;
    }
  }
  if (dynpoint.exists)
    fprintf(stdout, " - dyn:            %.5f s\n", writer->dur_write_dyn);
  if (n_nodedata_files != 0)
    fprintf(stdout, " - node/vtx data:  %.5f s\n", writer->dur_write_nodedat);
  if (n_elemdata_files != 0)
    fprintf(stdout, " - elem data:      %.5f s\n", writer->dur_write_elemdat);

  //Cleanup
  for(size_t i=0; i < n_nodedata_files; i++) {
    if(nodedata[i].fmt == IGB)
      fclose(nodedata[i].hdr.fileptr);
  }

  for(size_t i=0; i < n_elemdata_files; i++) {
    if(elemdata[i].fmt == IGB)
      fclose(elemdata[i].hdr.fileptr);
  }

  delete writer;
}

