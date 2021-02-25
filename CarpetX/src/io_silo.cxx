#include "io_silo.hxx"

#include "driver.hxx"
#include "timer.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <AMReX.H>
#include <AMReX_IntVect.H>

#include <mpi.h>

#ifdef HAVE_CAPABILITY_Silo
#include <silo.hxx>
#endif

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#ifdef HAVE_CAPABILITY_Silo
namespace CarpetX {
using namespace std;

constexpr bool io_verbose = true;

struct mesh_props_t {
  amrex::IntVect ngrow;

  auto to_tuple() const { return make_tuple(ngrow); }
  friend bool operator==(const mesh_props_t &p, const mesh_props_t &q) {
    return p.to_tuple() == q.to_tuple();
  }
  friend bool operator<(const mesh_props_t &p, const mesh_props_t &q) {
    return p.to_tuple() < q.to_tuple();
  }
};

string make_subdirname(const string &simulation_name, const int iteration) {
  ostringstream buf;
  buf << simulation_name                               //
      << ".it" << setw(8) << setfill('0') << iteration //
      << ".silo.dir";
  return buf.str();
}

string make_filename(const string &simulation_name, const int iteration,
                     const int ioserver = -1) {
  ostringstream buf;
  buf << simulation_name //
      << ".it" << setw(8) << setfill('0') << iteration;
  if (ioserver >= 0)
    buf << ".p" << setw(6) << setfill('0') << ioserver;
  buf << ".silo";
  return buf.str();
}

string make_meshname(const int reflevel = -1, const int component = -1) {
  assert((reflevel == -1) == (component == -1));
  ostringstream buf;
  if (reflevel < 0)
    buf << "gh";
  else
    buf << "box"                                        //
        << ".rl" << setw(2) << setfill('0') << reflevel //
        << ".c" << setw(8) << setfill('0') << component;
  return DB::legalize_name(buf.str());
}

string make_varname(const int gi, const int vi, const int reflevel = -1,
                    const int component = -1) {
  assert((reflevel == -1) == (component == -1));
  string varname;
  if (vi < 0) {
    assert(0);
    char *const groupname = CCTK_GroupName(gi);
    varname = groupname;
    free(groupname);
  } else {
    const int v0 = CCTK_FirstVarIndexI(gi);
    varname = CCTK_FullVarName(v0 + vi);
    varname = regex_replace(varname, regex("::"), "-");
    for (auto &ch : varname)
      ch = tolower(ch);
  }
  ostringstream buf;
  buf << varname;
  if (reflevel >= 0)
    buf << ".rl" << setw(2) << setfill('0') << reflevel //
        << ".c" << setw(8) << setfill('0') << component;
  return DB::legalize_name(buf.str());
}

const string driver_name = "CarpetX";

string make_fabarraybasename(const int reflevel) {
  ostringstream buf;
  buf << "FabArrayBase.rl" << setw(2) << setfill('0') << reflevel;
  return DB::legalize_name(buf.str());
}

////////////////////////////////////////////////////////////////////////////////

void InputSilo(const cGH *restrict const cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int ierr;

  // Set up timers
  static Timer timer("InputSilo");
  Interval interval(timer);

  static Timer timer_setup("InputSilo.setup");
  auto interval_setup = make_unique<Interval>(timer_setup);

  const MPI_Comm mpi_comm = MPI_COMM_WORLD;
  const int myproc = CCTK_MyProc(cctkGH);
  const int nprocs = CCTK_nProcs(cctkGH);

  const int ioproc_every = [&]() {
    if (CCTK_EQUALS(out_mode, "proc"))
      return 1;
    if (CCTK_EQUALS(out_mode, "np"))
      return out_proc_every;
    if (CCTK_EQUALS(out_mode, "onefile"))
      return nprocs;
    assert(0);
  }();
  assert(ioproc_every > 0);

  const int myioproc = myproc / ioproc_every * ioproc_every;
  assert(myioproc <= myproc);
  assert(myproc < myioproc + ioproc_every);
  const bool read_file = myproc % ioproc_every == 0;
  assert((myioproc == myproc) == read_file);

  const int metafile_ioproc = 0;
  const bool read_metafile = myproc == metafile_ioproc;

  // Configure Silo library
  DBShowErrors(DB_ALL_AND_DRVR, nullptr);
  // DBSetAllowEmptyObjects(1);
  DBSetCompression("METHOD=GZIP");
  DBSetEnableChecksums(1);

  // Determine input file name
  const string parfilename = [&]() {
    string buf(1024, '\0');
    const int len =
        CCTK_ParameterFilename(buf.length(), const_cast<char *>(buf.data()));
    buf.resize(len);
    return buf;
  }();

  // TODO: Check wether this parameter contains multiple file names
  // const string simulation_name = filereader_ID_files;
  const string simulation_name = in_file;

  // TODO
  // const string pathname = string(filereader_ID_dir);
  const string pathname = string(in_dir);

  // TODO: directories instead of carefully chosen names

  // Find input groups
  const vector<bool> group_enabled = [&]() {
    vector<bool> enabled(CCTK_NumGroups(), false);
    const auto callback{
        [](const int index, const char *const optstring, void *const arg) {
          vector<bool> &enabled = *static_cast<vector<bool> *>(arg);
          enabled.at(CCTK_GroupIndexFromVarI(index)) = true;
        }};
    CCTK_TraverseString(in_silo_vars, callback, &enabled, CCTK_GROUP_OR_VAR);
    if (io_verbose) {
      const bool have_silo_input = any_of(
          group_enabled.begin(), group_enabled.end(), [](bool b) { return b; });
      if (have_silo_input) {
        CCTK_VINFO("Silo input for groups:");
        for (int gi = 0; gi < CCTK_NumGroups(); ++gi) {
          if (group_enabled.at(gi)) {
            char *const groupname = CCTK_GroupName(gi);
            CCTK_VINFO("  %s", groupname);
            free(groupname);
          }
        }
      } else {
        CCTK_VINFO("No Silo input");
      }
    }
    return enabled;
  }();
  const auto num_in_vars =
      count(group_enabled.begin(), group_enabled.end(), true);
  if (num_in_vars == 0)
    return;

  constexpr int ndims = dim;

  interval_setup = nullptr;

  static Timer timer_meta("InputSilo.meta");
  auto interval_meta = make_unique<Interval>(timer_meta);

  // Read metadata
  if (read_metafile) {

    const string metafilename =
        pathname + "/" + make_filename(simulation_name, cctk_iteration);
    // We could use DB_UNKNOWN instead of DB_HDF5
    const DB::ptr<DBfile> metafile =
        DB::make(DBOpen(metafilename.c_str(), DB_HDF5, DB_READ));
    assert(metafile);

    // Read internal driver state
    // TODO: Only do this when recovering. But then, we should also do
    // this as a file reader option!
    {
      const string dirname = DB::legalize_name(driver_name);

      // Read number of levels
      const int nlevels = [&] {
        const string varname = dirname + "/" + DB::legalize_name("nlevels");
        const int vartype = DBGetVarType(metafile.get(), varname.c_str());
        assert(vartype == DB_INT);
        const int varlength = DBGetVarLength(metafile.get(), varname.c_str());
        assert(varlength == 1);
        int result;
        ierr = DBReadVar(metafile.get(), varname.c_str(), &result);
        assert(!ierr);
        return result;
      }();

      // Read FabArrayBase (component positions and shapes)
      vector<vector<amrex::Box> > boxes(nlevels);
      for (int level = 0; level < nlevels; ++level) {
        const string varname = dirname + "/" + make_fabarraybasename(level);
        const int vartype = DBGetVarType(metafile.get(), varname.c_str());
        assert(vartype == DB_INT);
        int vardims[2];
        const int varndims =
            DBGetVarDims(metafile.get(), varname.c_str(), 2, vardims);
        assert(varndims >= 0);
        assert(varndims == 2);
        assert(vardims[1] == 2 * ndims);
        const int nfabs = vardims[0];
        assert(nfabs >= 0);
        vector<int> data(2 * ndims * nfabs);
        ierr = DBReadVar(metafile.get(), varname.c_str(), data.data());

        auto &levboxes = boxes.at(level);
        levboxes.resize(nfabs);
        assert(!ierr);
        for (int component = 0; component < nfabs; ++component) {
          const amrex::IntVect small(&data.at(2 * ndims * component));
          const amrex::IntVect big(&data.at(ndims + 2 * ndims * component));
          levboxes.at(component) = amrex::Box(small, big);
        }
      }
    }
  }

  interval_meta = nullptr;

  static Timer timer_data("InputSilo.data");
  auto interval_data = make_unique<Interval>(timer_data);

  // Read data
  {
    DB::ptr<DBfile> file;
    if (read_file) {
      const string subdirname =
          make_subdirname(simulation_name, cctk_iteration);
      const string filename =
          pathname + "/" + subdirname + "/" +
          make_filename(simulation_name, cctk_iteration, myproc / ioproc_every);
      // We could use DB_UNKNOWN instead of DB_HDF5
      file = DB::make(DBOpen(filename.c_str(), DB_HDF5, DB_READ));
      assert(file);
    }

    // Loop over levels
    for (const auto &leveldata : ghext->leveldata) {
      if (io_verbose)
        CCTK_VINFO("Reading level %d", leveldata.level);

      // Loop over groups
      set<mesh_props_t> have_meshes;
      for (int gi = 0; gi < CCTK_NumGroups(); ++gi) {
        if (!group_enabled.at(gi))
          continue;
        if (CCTK_GroupTypeI(gi) != CCTK_GF)
          continue;
        if (io_verbose) {
          char *const groupname = CCTK_GroupName(gi);
          CCTK_VINFO("  Reading group %s", groupname);
          free(groupname);
        }

        auto &groupdata = *leveldata.groupdata.at(gi);
        const int numvars = groupdata.numvars;
        const int tl = 0;
        amrex::MultiFab &mfab = *groupdata.mfab[tl];
        const amrex::IndexType &indextype = mfab.ixType();
        const amrex::IntVect &ngrow = mfab.nGrowVect();
        const amrex::DistributionMapping &dm = mfab.DistributionMap();

        const mesh_props_t mesh_props{ngrow};
        const bool have_mesh = have_meshes.count(mesh_props);

        // Loop over components (AMReX boxes)
        const int nfabs = dm.size();
        for (int component = 0; component < nfabs; ++component) {
          if (io_verbose)
            CCTK_VINFO("    Reading component %d", component);

          const int proc = dm[component];
          const int ioproc = proc / ioproc_every * ioproc_every;
          const bool recv_this_fab = proc == myproc;
          const bool read_this_fab = ioproc == myproc;
          if (!(recv_this_fab || read_this_fab))
            continue;

          const amrex::Box &fabbox = mfab.fabbox(component); // exterior

          array<int, ndims> dims;
          for (int d = 0; d < ndims; ++d)
            dims[d] = fabbox.length(d);
          ptrdiff_t zonecount = 1;
          for (int d = 0; d < ndims; ++d)
            zonecount *= dims[d];
          assert(zonecount >= 0 && zonecount <= INT_MAX);

          // Communicate variable, part 1
          static Timer timer_mpi("InputSilo.mpi");
          auto interval_mpi = make_unique<Interval>(timer_mpi);
          const int mpi_tag = 22901; // randomly chosen
          vector<double> buffer;
          MPI_Request mpi_req;
          double *data = nullptr;
          if (recv_this_fab && read_this_fab) {
            amrex::FArrayBox &fab = mfab[component];
            data = fab.dataPtr();
          } else if (recv_this_fab) {
            amrex::FArrayBox &fab = mfab[component];
            assert(numvars * zonecount <= INT_MAX);
            MPI_Irecv(fab.dataPtr(), numvars * zonecount, MPI_DOUBLE, ioproc,
                      mpi_tag, mpi_comm, &mpi_req);
          } else {
            buffer.resize(numvars * zonecount);
            assert(numvars * zonecount <= INT_MAX);
            data = buffer.data();
          }
          interval_mpi = nullptr;

          // Read variable
          if (read_file) {
            static Timer timer_var("InputSilo.var");
            Interval interval_var(timer_var);

            const string meshname = make_meshname(leveldata.level, component);

            const int centering = [&]() {
              if (indextype.nodeCentered())
                return DB_NODECENT;
              if (indextype.cellCentered())
                return DB_ZONECENT;
              assert(0);
            }();

            for (int vi = 0; vi < numvars; ++vi) {
              const string varname =
                  make_varname(gi, vi, leveldata.level, component);
              if (io_verbose)
                CCTK_VINFO("      Reading variable %s", varname.c_str());

              const DB::ptr<DBquadvar> quadvar =
                  DB::make(DBGetQuadvar(file.get(), varname.c_str()));
              assert(quadvar);

              assert(quadvar->ndims == ndims);
              assert(ndims <= 3);
              for (int d = 0; d < ndims; ++d)
                assert(quadvar->dims[d] == dims[d]);
              assert(quadvar->datatype == DB_DOUBLE);
              assert(quadvar->centering == centering);
              assert(quadvar->nvals == 1);

              // TODO: check DBOPT_COORDSYS: int cartesian = DB_CARTESIAN;
              const int column_major = 0;
              assert(quadvar->major_order == column_major);

              const void *const read_ptr = quadvar->vals[0];

              void *const data_ptr = data + vi * zonecount;
              memcpy(data_ptr, read_ptr, zonecount * sizeof(double));
            } // for vi
          }   // if read_file

          // Communicate variable, part 2
          static Timer timer_wait("InputSilo.wait");
          auto interval_wait = make_unique<Interval>(timer_wait);
          if (recv_this_fab && read_this_fab) {
            // do nothing
          } else if (recv_this_fab) {
            MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
          } else {
            buffer.resize(numvars * zonecount);
            assert(numvars * zonecount <= INT_MAX);
            MPI_Send(buffer.data(), numvars * zonecount, MPI_DOUBLE, proc,
                     mpi_tag, mpi_comm);
          }
          if (recv_this_fab)
            for (int vi = 0; vi < numvars; ++vi)
              groupdata.valid.at(tl).at(vi).set(
                  valid_t(true), []() { return "read from file"; });
          interval_wait = nullptr;

        } // for component

      } // for gi
    }   // for leveldata
  }     // write data

  interval_data = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

void OutputSilo(const cGH *restrict const cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int ierr;

  // Set up timers
  static Timer timer("OutputSilo");
  Interval interval(timer);

  static Timer timer_setup("OutputSilo.setup");
  auto interval_setup = make_unique<Interval>(timer_setup);

  const MPI_Comm mpi_comm = MPI_COMM_WORLD;
  const int myproc = CCTK_MyProc(cctkGH);
  const int nprocs = CCTK_nProcs(cctkGH);

  const int ioproc_every = [&]() {
    if (CCTK_EQUALS(out_mode, "proc"))
      return 1;
    if (CCTK_EQUALS(out_mode, "np"))
      return out_proc_every;
    if (CCTK_EQUALS(out_mode, "onefile"))
      return nprocs;
    assert(0);
  }();
  assert(ioproc_every > 0);

  const int myioproc = myproc / ioproc_every * ioproc_every;
  assert(myioproc <= myproc);
  assert(myproc < myioproc + ioproc_every);
  const bool write_file = myproc % ioproc_every == 0;
  assert((myioproc == myproc) == write_file);
  // If process 1 exists and if it is not an I/O process, then output
  // metadata there. Else output metadata on process 0.

  const int metafile_ioproc = nprocs == 1 || ioproc_every == 1 ? 0 : 1;
  const bool write_metafile = myproc == metafile_ioproc;

  // Configure Silo library
  DBShowErrors(DB_ALL_AND_DRVR, nullptr);
  // DBSetAllowEmptyObjects(1);
  DBSetCompression("METHOD=GZIP");
  DBSetEnableChecksums(1);

  // Determine output file name
  const string parfilename = [&]() {
    string buf(1024, '\0');
    const int len =
        CCTK_ParameterFilename(buf.length(), const_cast<char *>(buf.data()));
    buf.resize(len);
    return buf;
  }();

  const string simulation_name = [&]() {
    string name = parfilename;
    const size_t last_slash = name.rfind('/');
    if (last_slash != string::npos && last_slash < name.length())
      name = name.substr(last_slash + 1);
    const size_t last_dot = name.rfind('.');
    if (last_dot != string::npos && last_dot > 0)
      name = name.substr(0, last_dot);
    return name;
  }();

  // TODO: directories instead of carefully chosen names

  // Find output groups
  const vector<bool> group_enabled = [&] {
    vector<bool> enabled(CCTK_NumGroups(), false);
    const auto callback{
        [](const int index, const char *const optstring, void *const arg) {
          vector<bool> &enabled = *static_cast<vector<bool> *>(arg);
          enabled.at(CCTK_GroupIndexFromVarI(index)) = true;
        }};
    CCTK_TraverseString(out_silo_vars, callback, &enabled, CCTK_GROUP_OR_VAR);
    if (io_verbose) {
      CCTK_VINFO("Silo output for groups:");
      for (int gi = 0; gi < CCTK_NumGroups(); ++gi) {
        if (group_enabled.at(gi)) {
          char *const groupname = CCTK_GroupName(gi);
          CCTK_VINFO("  %s", groupname);
          free(groupname);
        }
      }
    }
    return enabled;
  }();
  const auto num_out_vars =
      count(group_enabled.begin(), group_enabled.end(), true);
  if (num_out_vars == 0)
    return;

  constexpr int ndims = dim;

  interval_setup = nullptr;

  static Timer timer_data("OutputSilo.data");
  auto interval_data = make_unique<Interval>(timer_data);

  // Write data
  {
    static Timer timer_mkdir("OutputSilo.mkdir");
    auto interval_mkdir = make_unique<Interval>(timer_mkdir);
    const string subdirname = make_subdirname(simulation_name, cctk_iteration);
    const string pathname = string(out_dir) + "/" + subdirname;
    const int mode = 0755;
    static bool did_create_directory = false;
    if (!did_create_directory) {
      ierr = CCTK_CreateDirectory(mode, out_dir);
      assert(ierr >= 0);
      did_create_directory = true;
    }
    ierr = CCTK_CreateDirectory(mode, pathname.c_str());
    assert(ierr >= 0);
    interval_mkdir = nullptr;

    DB::ptr<DBfile> file;
    if (write_file) {
      const string filename =
          pathname + "/" +
          make_filename(simulation_name, cctk_iteration, myproc / ioproc_every);
      file = DB::make(DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL,
                               simulation_name.c_str(), DB_HDF5));
      assert(file);
    }

    {
      // Tell VisIt that the mesh structure may change over time
      const int dims = 1;
      const int value = 1;
      ierr = DBWrite(file.get(), "MetadataIsTimeVarying", &value, &dims, 1,
                     DB_INT);
      assert(!ierr);
    }

    // Loop over levels
    for (const auto &leveldata : ghext->leveldata) {

      // Loop over groups
      set<mesh_props_t> have_meshes;
      for (int gi = 0; gi < CCTK_NumGroups(); ++gi) {
        if (!group_enabled.at(gi))
          continue;
        if (CCTK_GroupTypeI(gi) != CCTK_GF)
          continue;

        const auto &groupdata = *leveldata.groupdata.at(gi);
        const int numvars = groupdata.numvars;
        const int tl = 0;
        const amrex::MultiFab &mfab = *groupdata.mfab[tl];
        const amrex::IndexType &indextype = mfab.ixType();
        const amrex::IntVect &ngrow = mfab.nGrowVect();
        const amrex::DistributionMapping &dm = mfab.DistributionMap();

        const mesh_props_t mesh_props{ngrow};
        const bool have_mesh = have_meshes.count(mesh_props);

        // Loop over components (AMReX boxes)
        const int nfabs = dm.size();
        for (int component = 0; component < nfabs; ++component) {
          const int proc = dm[component];
          const int ioproc = proc / ioproc_every * ioproc_every;
          const bool send_this_fab = proc == myproc;
          const bool write_this_fab = ioproc == myproc;
          if (!(send_this_fab || write_this_fab))
            continue;

          // TODO: Check whether data are valid
          const amrex::Box &fabbox = mfab.fabbox(component); // exterior

          array<int, ndims> dims;
          for (int d = 0; d < ndims; ++d)
            dims[d] = fabbox.length(d);
          ptrdiff_t zonecount = 1;
          for (int d = 0; d < ndims; ++d)
            zonecount *= dims[d];
          assert(zonecount >= 0 && zonecount <= INT_MAX);

          if (write_file && !have_mesh) {
            static Timer timer_mesh("OutputSilo.mesh");
            Interval interval_mesh(timer_mesh);

            const string meshname = make_meshname(leveldata.level, component);

            array<int, ndims> dims_vc;
            for (int d = 0; d < ndims; ++d)
              dims_vc[d] = dims[d] + int(indextype.cellCentered(d));

            const amrex::Geometry &geom = ghext->amrcore->Geom(leveldata.level);
            const double *const x0 = geom.ProbLo();
            const double *const dx = geom.CellSize();
            array<vector<double>, ndims> coords;
            for (int d = 0; d < ndims; ++d) {
              coords[d].resize(dims_vc[d]);
              for (int i = 0; i < dims_vc[d]; ++i)
                coords[d][i] = x0[d] + (fabbox.smallEnd(d) + i) * dx[d];
            }
            array<void *, ndims> coord_ptrs;
            for (int d = 0; d < ndims; ++d)
              coord_ptrs[d] = coords[d].data();

            const DB::ptr<DBoptlist> optlist = DB::make(DBMakeOptlist(10));
            assert(optlist);

            int cartesian = DB_CARTESIAN;
            ierr = DBAddOption(optlist.get(), DBOPT_COORDSYS, &cartesian);
            assert(!ierr);

            int cycle = cctk_iteration;
            ierr = DBAddOption(optlist.get(), DBOPT_CYCLE, &cycle);
            assert(!ierr);

            array<int, ndims> min_index, max_index;
            for (int d = 0; d < ndims; ++d) {
              min_index[d] = ngrow[d];
              max_index[d] = ngrow[d];
            }
            ierr =
                DBAddOption(optlist.get(), DBOPT_LO_OFFSET, min_index.data());
            assert(!ierr);
            ierr =
                DBAddOption(optlist.get(), DBOPT_HI_OFFSET, max_index.data());
            assert(!ierr);

            int column_major = 0;
            ierr = DBAddOption(optlist.get(), DBOPT_MAJORORDER, &column_major);
            assert(!ierr);

            // float time = cctk_time;
            // ierr = DBAddOption(optlist.get(), DBOPT_TIME, &time);
            // assert(!ierr);
            double dtime = cctk_time;
            ierr = DBAddOption(optlist.get(), DBOPT_DTIME, &dtime);
            assert(!ierr);

            int hide_from_gui = 1;
            ierr =
                DBAddOption(optlist.get(), DBOPT_HIDE_FROM_GUI, &hide_from_gui);
            assert(!ierr);

            ierr = DBPutQuadmesh(file.get(), meshname.c_str(), nullptr,
                                 coord_ptrs.data(), dims_vc.data(), ndims,
                                 DB_DOUBLE, DB_COLLINEAR, optlist.get());
            assert(!ierr);
          } // if write mesh

          // Communicate variable
          static Timer timer_mpi("OutputSilo.mpi");
          auto interval_mpi = make_unique<Interval>(timer_mpi);
          const int mpi_tag = 22900; // randomly chosen
          vector<double> buffer;
          const double *data = nullptr;
          if (send_this_fab && write_this_fab) {
            const amrex::FArrayBox &fab = mfab[component];
            data = fab.dataPtr();
          } else if (send_this_fab) {
            const amrex::FArrayBox &fab = mfab[component];
            assert(numvars * zonecount <= INT_MAX);
            MPI_Send(fab.dataPtr(), numvars * zonecount, MPI_DOUBLE, ioproc,
                     mpi_tag, mpi_comm);
          } else {
            buffer.resize(numvars * zonecount);
            assert(numvars * zonecount <= INT_MAX);
            MPI_Recv(buffer.data(), numvars * zonecount, MPI_DOUBLE, proc,
                     mpi_tag, mpi_comm, MPI_STATUS_IGNORE);
            data = buffer.data();
          }
          interval_mpi = nullptr;

          // Write variable
          if (write_file) {
            static Timer timer_var("OutputSilo.var");
            Interval interval_var(timer_var);

            const string meshname = make_meshname(leveldata.level, component);

            const int centering = [&]() {
              if (indextype.nodeCentered())
                return DB_NODECENT;
              if (indextype.cellCentered())
                return DB_ZONECENT;
              assert(0);
            }();

            const DB::ptr<DBoptlist> optlist = DB::make(DBMakeOptlist(10));
            assert(optlist);

            int cartesian = DB_CARTESIAN;
            ierr = DBAddOption(optlist.get(), DBOPT_COORDSYS, &cartesian);
            assert(!ierr);

            int cycle = cctk_iteration;
            ierr = DBAddOption(optlist.get(), DBOPT_CYCLE, &cycle);
            assert(!ierr);

            int column_major = 0;
            ierr = DBAddOption(optlist.get(), DBOPT_MAJORORDER, &column_major);
            assert(!ierr);

            // float time = cctk_time;
            // ierr = DBAddOption(optlist.get(), DBOPT_TIME, &time);
            // assert(!ierr);
            double dtime = cctk_time;
            ierr = DBAddOption(optlist.get(), DBOPT_DTIME, &dtime);
            assert(!ierr);

            int hide_from_gui = 1;
            ierr =
                DBAddOption(optlist.get(), DBOPT_HIDE_FROM_GUI, &hide_from_gui);
            assert(!ierr);

            for (int vi = 0; vi < numvars; ++vi) {
              const string varname =
                  make_varname(gi, vi, leveldata.level, component);

              const void *const data_ptr = data + vi * zonecount;

              ierr =
                  DBPutQuadvar1(file.get(), varname.c_str(), meshname.c_str(),
                                data_ptr, dims.data(), ndims, nullptr, 0,
                                DB_DOUBLE, centering, optlist.get());
              assert(!ierr);
            } // for vi
          }   // if write_file

        } // for component

      } // for gi
    }   // for leveldata
  }     // write data

  interval_data = nullptr;

  static Timer timer_meta("OutputSilo.meta");
  auto interval_meta = make_unique<Interval>(timer_meta);

  // Write metadata
  if (write_metafile) {

    const string metafilename =
        string(out_dir) + "/" + make_filename(simulation_name, cctk_iteration);
    const DB::ptr<DBfile> metafile =
        DB::make(DBCreate(metafilename.c_str(), DB_CLOBBER, DB_LOCAL,
                          simulation_name.c_str(), DB_HDF5));
    assert(metafile);

    {
      // Tell VisIt that the mesh structure may change over time
      const int dims = 1;
      const int value = 1;
      ierr = DBWrite(metafile.get(), "MetadataIsTimeVarying", &value, &dims, 1,
                     DB_INT);
      assert(!ierr);
    }

    // Loop over groups
    set<mesh_props_t> have_meshes;
    for (int gi = 0; gi < CCTK_NumGroups(); ++gi) {
      if (!group_enabled.at(gi))
        continue;
      if (CCTK_GroupTypeI(gi) != CCTK_GF)
        continue;

      const auto &leveldata0 = ghext->leveldata.at(0);
      const auto &groupdata0 = *leveldata0.groupdata.at(gi);
      const int numvars = groupdata0.numvars;
      const int tl = 0;
      const amrex::MultiFab &mfab0 = *groupdata0.mfab[tl];
      const amrex::IndexType &indextype = mfab0.ixType();
      const amrex::IntVect &ngrow = mfab0.nGrowVect();

      const mesh_props_t mesh_props{ngrow};
      const bool have_mesh = have_meshes.count(mesh_props);

      if (!have_mesh) {

        const string multimeshname = make_meshname();

        // Count components per level
        const int nlevels = ghext->leveldata.size();
        vector<int> ncomps_level;
        for (const auto &leveldata : ghext->leveldata) {
          const auto &groupdata = *leveldata.groupdata.at(gi);
          const amrex::MultiFab &mfab = *groupdata.mfab[tl];
          const amrex::DistributionMapping &dm = mfab.DistributionMap();
          const int nfabs = dm.size();
          ncomps_level.push_back(nfabs);
        }
        vector<int> firstcomp_level;
        int ncomps_total = 0;
        for (const int ncomps : ncomps_level) {
          firstcomp_level.push_back(ncomps_total);
          ncomps_total += ncomps;
        }

        // Describe which components belong to which level
        // Question: Can this name be changed?
        const string levelmaps_name = multimeshname + "_wmrgtree_lvlMaps";
        {
          vector<int> segment_types;
          vector<vector<int> > segment_data;
          segment_types.reserve(nlevels);
          segment_data.reserve(nlevels);
          for (int l = 0; l < nlevels; ++l) {
            const int comp0 = firstcomp_level.at(l);
            const int ncomps = ncomps_level.at(l);
            vector<int> data;
            data.reserve(ncomps);
            for (int c = 0; c < ncomps; ++c)
              data.push_back(comp0 + c);
            segment_types.push_back(DB_BLOCKCENT);
            segment_data.push_back(move(data));
          }
          vector<int> segment_lengths;
          vector<const int *> segment_data_ptrs;
          segment_lengths.reserve(segment_data.size());
          segment_data_ptrs.reserve(segment_data.size());
          for (const auto &data : segment_data) {
            segment_lengths.push_back(data.size());
            segment_data_ptrs.push_back(data.data());
          }

          ierr = DBPutGroupelmap(metafile.get(), levelmaps_name.c_str(),
                                 nlevels, segment_types.data(),
                                 segment_lengths.data(), nullptr,
                                 segment_data_ptrs.data(), nullptr, 0, nullptr);
          assert(!ierr);
        }

        // Describe which components are children of which other components
        // Question: Can this name be changed?
        const string childmaps_name = multimeshname + "_wmrgtree_chldMaps";
        vector<int> num_children;
        {
          vector<int> segment_types;
          vector<vector<int> > segment_data;
          segment_types.reserve(ncomps_total);
          segment_data.reserve(ncomps_total);
          for (const auto &leveldata : ghext->leveldata) {
            const int level = leveldata.level;
            const int fine_level = level + 1;
            if (fine_level < nlevels) {
              const auto &groupdata = *leveldata.groupdata.at(gi);
              const amrex::MultiFab &mfab = *groupdata.mfab[tl];
              const auto &fine_leveldata = ghext->leveldata.at(fine_level);
              const auto &fine_groupdata = *fine_leveldata.groupdata.at(gi);
              const amrex::MultiFab &fine_mfab = *fine_groupdata.mfab[tl];

              const int ncomps = ncomps_level.at(level);
              const int fine_comp0 = firstcomp_level.at(fine_level);
              const amrex::BoxArray &fine_boxarray = fine_mfab.boxarray;

              for (int component = 0; component < ncomps; ++component) {
                const amrex::Box &box = mfab.box(component); // interior
                amrex::Box refined_box(box);
                refined_box.refine(2);
                const vector<pair<int, amrex::Box> > child_boxes =
                    fine_boxarray.intersections(refined_box);
                vector<int> children;
                children.reserve(child_boxes.size());
                for (const auto &ib : child_boxes) {
                  const int fine_component = ib.first;
                  children.push_back(fine_comp0 + fine_component);
                }

                segment_types.push_back(DB_BLOCKCENT);
                segment_data.push_back(move(children));
              }

            } else {
              // no finer level, hence no children
              const int ncomps = ncomps_level.at(level);
              for (int component = 0; component < ncomps; ++component) {
                segment_types.push_back(DB_BLOCKCENT);
                segment_data.emplace_back();
              }
            }
          }

          vector<int> &segment_lengths = num_children;
          vector<const int *> segment_data_ptrs;
          segment_lengths.reserve(segment_data.size());
          segment_data_ptrs.reserve(segment_data.size());
          for (const auto &data : segment_data) {
            segment_lengths.push_back(data.size());
            segment_data_ptrs.push_back(data.data());
          }

          ierr = DBPutGroupelmap(metafile.get(), childmaps_name.c_str(),
                                 ncomps_total, segment_types.data(),
                                 segment_lengths.data(), nullptr,
                                 segment_data_ptrs.data(), nullptr, 0, nullptr);
          assert(!ierr);
        }

        // Create Mrgtree
        {
          const int max_mgrtree_children = 2;
          const DB::ptr<DBmrgtree> mrgtree = DB::make(
              DBMakeMrgtree(DB_MULTIMESH, 0, max_mgrtree_children, nullptr));
          assert(mrgtree);

          // Describe AMR configuration
          const int max_amr_decomp_children = 2;
          ierr = DBAddRegion(mrgtree.get(), "amr_decomp", 0,
                             max_amr_decomp_children, nullptr, 0, nullptr,
                             nullptr, nullptr, nullptr);
          assert(!ierr);
          ierr = DBSetCwr(mrgtree.get(), "amr_decomp");
          assert(ierr >= 0);

          // Describe AMR levels
          {
            ierr = DBAddRegion(mrgtree.get(), "levels", 0, nlevels, nullptr, 0,
                               nullptr, nullptr, nullptr, nullptr);
            assert(!ierr);

            ierr = DBSetCwr(mrgtree.get(), "levels");
            assert(ierr >= 0);

            const vector<string> region_names{"@level%d@n"};
            vector<const char *> region_name_ptrs;
            region_name_ptrs.reserve(region_names.size());
            for (const string &name : region_names)
              region_name_ptrs.push_back(name.c_str());

            vector<int> segment_ids;
            vector<int> segment_types;
            segment_ids.reserve(nlevels);
            segment_types.reserve(nlevels);
            for (int l = 0; l < nlevels; ++l) {
              segment_ids.push_back(l);
              segment_types.push_back(DB_BLOCKCENT);
            }

            ierr = DBAddRegionArray(
                mrgtree.get(), nlevels, region_name_ptrs.data(), 0,
                levelmaps_name.c_str(), 1, segment_ids.data(),
                ncomps_level.data(), segment_types.data(), nullptr);
            assert(!ierr);

            ierr = DBSetCwr(mrgtree.get(), "..");
            assert(ierr >= 0);
          }

          // Describe AMR children
          {
            ierr = DBAddRegion(mrgtree.get(), "patches", 0, ncomps_total,
                               nullptr, 0, nullptr, nullptr, nullptr, nullptr);
            assert(ierr >= 0);

            ierr = DBSetCwr(mrgtree.get(), "patches");
            assert(ierr >= 0);

            const vector<string> region_names{"@patch%d@n"};
            vector<const char *> region_name_ptrs;
            region_name_ptrs.reserve(region_names.size());
            for (const string &name : region_names)
              region_name_ptrs.push_back(name.c_str());

            vector<int> segment_types;
            vector<int> segment_ids;
            segment_types.reserve(ncomps_total);
            segment_ids.reserve(ncomps_total);
            for (int c = 0; c < ncomps_total; ++c) {
              segment_ids.push_back(c);
              segment_types.push_back(DB_BLOCKCENT);
            }

            ierr = DBAddRegionArray(
                mrgtree.get(), ncomps_total, region_name_ptrs.data(), 0,
                childmaps_name.c_str(), 1, segment_ids.data(),
                num_children.data(), segment_types.data(), nullptr);

            ierr = DBSetCwr(mrgtree.get(), "..");
            assert(ierr >= 0);
          }

          {
            const vector<string> mrgv_onames{
                multimeshname + "_wmrgtree_lvlRatios",
                multimeshname + "_wmrgtree_ijkExts",
                multimeshname + "_wmrgtree_xyzExts", "rank"};
            vector<const char *> mrgv_oname_ptrs;
            mrgv_oname_ptrs.reserve(mrgv_onames.size() + 1);
            for (const string &name : mrgv_onames)
              mrgv_oname_ptrs.push_back(name.c_str());
            mrgv_oname_ptrs.push_back(nullptr);

            const DB::ptr<DBoptlist> optlist = DB::make(DBMakeOptlist(10));
            assert(optlist);

            ierr = DBAddOption(optlist.get(), DBOPT_MRGV_ONAMES,
                               mrgv_oname_ptrs.data());
            assert(!ierr);

            ierr = DBPutMrgtree(metafile.get(), "mrgTree", "amr_mesh",
                                mrgtree.get(), optlist.get());
            assert(!ierr);
          }
        }

        // Write refinement ratios
        {
          const string levelrationame = multimeshname + "_wmrgtree_lvlRatios";

          const vector<string> compnames{"iRatio", "jRatio", "kRatio"};
          vector<const char *> compname_ptrs;
          compname_ptrs.reserve(compnames.size());
          for (const string &name : compnames)
            compname_ptrs.push_back(name.c_str());

          const vector<string> regionnames{"@level%d@n"};
          vector<const char *> regionname_ptrs;
          regionname_ptrs.reserve(regionnames.size());
          for (const string &name : regionnames)
            regionname_ptrs.push_back(name.c_str());

          array<vector<int>, ndims> data;
          for (int d = 0; d < ndims; ++d) {
            data[d].reserve(1);
            data[d].push_back(2);
          }
          array<const void *, ndims> data_ptrs;
          for (int d = 0; d < ndims; ++d)
            data_ptrs[d] = data[d].data();

          ierr = DBPutMrgvar(metafile.get(), levelrationame.c_str(), "mrgTree",
                             ndims, compname_ptrs.data(), nlevels,
                             regionname_ptrs.data(), DB_INT, data_ptrs.data(),
                             nullptr);
          assert(!ierr);
        }

        typedef array<array<int, ndims>, 2> iextent_t;
        typedef array<array<double, ndims>, 2> extent_t;
        vector<iextent_t> iextents;
        vector<extent_t> extents;
        iextents.reserve(ncomps_total);
        extents.reserve(ncomps_total);
        for (const auto &leveldata : ghext->leveldata) {
          const auto &groupdata = *leveldata.groupdata.at(gi);
          const int tl = 0;
          const amrex::MultiFab &mfab = *groupdata.mfab[tl];
          const amrex::Geometry &geom = ghext->amrcore->Geom(leveldata.level);
          const double *const x0 = geom.ProbLo();
          const double *const dx = geom.CellSize();
          const int nfabs = mfab.size();
          for (int c = 0; c < nfabs; ++c) {
            const amrex::Box &fabbox = mfab.fabbox(c); // exterior
            iextent_t iextent;
            extent_t extent;
            for (int d = 0; d < ndims; ++d) {
              iextent[0][d] = fabbox.smallEnd(d);
              iextent[1][d] = fabbox.bigEnd(d);
              extent[0][d] = x0[d] + fabbox.smallEnd(d) * dx[d];
              extent[1][d] = x0[d] + fabbox.bigEnd(d) * dx[d];
            }
            iextents.push_back(iextent);
            extents.push_back(extent);
          }
        }

        // Write extents
        {
          const string iextentsname = multimeshname + "_wmrgtree_ijkExts";
          const string extentsname = multimeshname + "_wmrgtree_xyzExts";

          const vector<string> icompnames{"iMin", "iMax", "jMin",
                                          "jMax", "kMin", "kMax"};
          const vector<string> compnames{"xMin", "xMax", "yMin",
                                         "yMax", "zMin", "zMax"};
          vector<const char *> icompname_ptrs;
          icompname_ptrs.reserve(icompnames.size());
          for (const string &name : icompnames)
            icompname_ptrs.push_back(name.c_str());
          vector<const char *> compname_ptrs;
          compname_ptrs.reserve(compnames.size());
          for (const string &name : compnames)
            compname_ptrs.push_back(name.c_str());

          const vector<string> regionnames{"@patch%d@n"};
          vector<const char *> regionname_ptrs;
          regionname_ptrs.reserve(regionnames.size());
          for (const string &name : regionnames)
            regionname_ptrs.push_back(name.c_str());

          array<array<vector<int>, 2>, ndims> idata;
          array<array<vector<double>, 2>, ndims> data;
          for (int d = 0; d < ndims; ++d) {
            for (int f = 0; f < 2; ++f) {
              idata[d][f].reserve(ncomps_total);
              data[d][f].reserve(ncomps_total);
            }
          }
          for (int c = 0; c < ncomps_total; ++c) {
            for (int d = 0; d < ndims; ++d) {
              for (int f = 0; f < 2; ++f) {
                idata[d][f].push_back(iextents[c][f][d]);
                data[d][f].push_back(extents[c][f][d]);
              }
            }
          }
          array<array<const void *, 2>, ndims> idata_ptrs;
          array<array<const void *, 2>, ndims> data_ptrs;
          for (int d = 0; d < ndims; ++d) {
            for (int f = 0; f < 2; ++f) {
              idata_ptrs[d][f] = idata[d][f].data();
              data_ptrs[d][f] = data[d][f].data();
            }
          }

          ierr = DBPutMrgvar(metafile.get(), iextentsname.c_str(), "mrgTree",
                             2 * ndims, icompname_ptrs.data(), ncomps_total,
                             regionname_ptrs.data(), DB_INT, idata_ptrs.data(),
                             nullptr);
          assert(!ierr);

          ierr = DBPutMrgvar(metafile.get(), extentsname.c_str(), "mrgTree",
                             2 * ndims, compname_ptrs.data(), ncomps_total,
                             regionname_ptrs.data(), DB_DOUBLE,
                             data_ptrs.data(), nullptr);
          assert(!ierr);

          // Write rank
          const vector<int> ranks(ncomps_total, ndims);
          const vector<const void *> rank_ptrs{ranks.data()};
          ierr = DBPutMrgvar(metafile.get(), "rank", "mrgTree", 1, nullptr,
                             ncomps_total, regionname_ptrs.data(), DB_INT,
                             rank_ptrs.data(), nullptr);
          assert(!ierr);
        }

        // Write multimesh

        vector<string> meshnames;
        for (const auto &leveldata : ghext->leveldata) {
          const auto &groupdata = *leveldata.groupdata.at(gi);
          const amrex::MultiFab &mfab = *groupdata.mfab[tl];
          const amrex::DistributionMapping &dm = mfab.DistributionMap();
          const int nfabs = dm.size();
          for (int c = 0; c < nfabs; ++c) {
            const int proc = dm[c];
            const string proc_filename =
                make_subdirname(simulation_name, cctk_iteration) + "/" +
                make_filename(simulation_name, cctk_iteration,
                              proc / ioproc_every);
            const string meshname =
                proc_filename + ":" + make_meshname(leveldata.level, c);
            meshnames.push_back(meshname);
          }
        }
        vector<const char *> meshname_ptrs;
        meshname_ptrs.reserve(meshnames.size());
        for (const auto &meshname : meshnames)
          meshname_ptrs.push_back(meshname.c_str());

        const DB::ptr<DBoptlist> optlist = DB::make(DBMakeOptlist(10));
        assert(optlist);

        int cycle = cctk_iteration;
        ierr = DBAddOption(optlist.get(), DBOPT_CYCLE, &cycle);
        assert(!ierr);

        // float time = cctk_time;
        // ierr = DBAddOption(optlist.get(), DBOPT_TIME, &time);
        // assert(!ierr);
        double dtime = cctk_time;
        ierr = DBAddOption(optlist.get(), DBOPT_DTIME, &dtime);
        assert(!ierr);

        int quadmesh = DB_QUADMESH;
        ierr = DBAddOption(optlist.get(), DBOPT_MB_BLOCK_TYPE, &quadmesh);
        assert(!ierr);

        int extents_size = 2 * ndims;
        // typedef array<array<double, ndims>, 2> extent_t;
        // vector<extent_t> extents;
        // extents.reserve(meshnames.size());
        // for (const auto &leveldata : ghext->leveldata) {
        //   const auto &groupdata = *leveldata.groupdata.at(gi);
        //   const int tl = 0;
        //   const amrex::MultiFab &mfab = *groupdata.mfab[tl];
        //   const amrex::Geometry &geom =
        //   ghext->amrcore->Geom(leveldata.level); const double *const x0 =
        //   geom.ProbLo(); const double *const dx = geom.CellSize(); const int
        //   nfabs = mfab.size(); for (int c = 0; c < nfabs; ++c) {
        //     const amrex::Box &fabbox = mfab.fabbox(c); // exterior
        //     extent_t extent;
        //     for (int d = 0; d < ndims; ++d) {
        //       extent[0][d] = x0[d] + fabbox.smallEnd(d) * dx[d];
        //       extent[1][d] = x0[d] + fabbox.bigEnd(d) * dx[d];
        //     }
        //     extents.push_back(extent);
        //   }
        // }
        ierr = DBAddOption(optlist.get(), DBOPT_EXTENTS_SIZE, &extents_size);
        assert(!ierr);
        assert(extents.size() == meshname_ptrs.size());
        ierr = DBAddOption(optlist.get(), DBOPT_EXTENTS, extents.data());
        assert(!ierr);

        vector<int> zonecounts;
        zonecounts.reserve(meshnames.size());
        for (const auto &leveldata : ghext->leveldata) {
          const auto &groupdata = *leveldata.groupdata.at(gi);
          const int tl = 0;
          const amrex::MultiFab &mfab = *groupdata.mfab[tl];
          const int nfabs = mfab.size();
          for (int c = 0; c < nfabs; ++c) {
            const amrex::Box &fabbox = mfab.fabbox(c); // exterior
            array<int, ndims> dims_vc;
            for (int d = 0; d < ndims; ++d)
              dims_vc[d] = fabbox.length(d) + int(indextype.cellCentered(d));
            int zonecount = 1;
            for (int d = 0; d < ndims; ++d)
              zonecount *= dims_vc[d];
            zonecounts.push_back(zonecount);
          }
        }
        assert(zonecounts.size() == meshname_ptrs.size());
        ierr = DBAddOption(optlist.get(), DBOPT_ZONECOUNTS, zonecounts.data());
        assert(!ierr);

        const string mrgtreename = "mrgtree";
        ierr = DBAddOption(optlist.get(), DBOPT_MRGTREE_NAME,
                           const_cast<char *>(mrgtreename.c_str()));
        assert(!ierr);

        ierr = DBPutMultimesh(metafile.get(), multimeshname.c_str(),
                              meshname_ptrs.size(), meshname_ptrs.data(),
                              nullptr, optlist.get());
        assert(!ierr);

        have_meshes.insert(mesh_props);
      } // if write multimesh

      // Write multivar
      {
        const string multimeshname = make_meshname();

        const DB::ptr<DBoptlist> optlist = DB::make(DBMakeOptlist(10));
        assert(optlist);

        int cycle = cctk_iteration;
        ierr = DBAddOption(optlist.get(), DBOPT_CYCLE, &cycle);
        assert(!ierr);

        // float time = cctk_time;
        // ierr = DBAddOption(optlist.get(), DBOPT_TIME, &time);
        // assert(!ierr);
        double dtime = cctk_time;
        ierr = DBAddOption(optlist.get(), DBOPT_DTIME, &dtime);
        assert(!ierr);

        ierr = DBAddOption(optlist.get(), DBOPT_MMESH_NAME,
                           const_cast<char *>(multimeshname.c_str()));
        assert(!ierr);

        int vartype_scalar = DB_VARTYPE_SCALAR;
        ierr = DBAddOption(optlist.get(), DBOPT_TENSOR_RANK, &vartype_scalar);
        assert(!ierr);

        int quadvar = DB_QUADVAR;
        ierr = DBAddOption(optlist.get(), DBOPT_MB_BLOCK_TYPE, &quadvar);
        assert(!ierr);

        for (int vi = 0; vi < numvars; ++vi) {
          const string multivarname = make_varname(gi, vi);

          vector<string> varnames;
          for (const auto &leveldata : ghext->leveldata) {
            const auto &groupdata = *leveldata.groupdata.at(gi);
            const int tl = 0;
            const amrex::MultiFab &mfab = *groupdata.mfab[tl];
            const amrex::DistributionMapping &dm = mfab.DistributionMap();
            const int nfabs = dm.size();
            for (int c = 0; c < nfabs; ++c) {
              const int proc = dm[c];
              const string proc_filename =
                  make_subdirname(simulation_name, cctk_iteration) + "/" +
                  make_filename(simulation_name, cctk_iteration,
                                proc / ioproc_every);
              const string varname = proc_filename + ":" +
                                     make_varname(gi, vi, leveldata.level, c);
              varnames.push_back(varname);
            }
          }
          vector<const char *> varname_ptrs;
          varname_ptrs.reserve(varnames.size());
          for (const auto &varname : varnames)
            varname_ptrs.push_back(varname.c_str());

          ierr = DBPutMultivar(metafile.get(), multivarname.c_str(),
                               varname_ptrs.size(), varname_ptrs.data(),
                               nullptr, optlist.get());
          assert(!ierr);
        } // for vi
      }   // write multivar

    } // for gi

    // Write internal driver state
    {
      const string dirname = DB::legalize_name(driver_name);
      ierr = DBMkDir(metafile.get(), dirname.c_str());
      assert(!ierr);

      // Write number of levels
      {
        const int dims = 1;
        const int value = ghext->leveldata.size();
        const string varname = dirname + "/" + DB::legalize_name("nlevels");
        ierr =
            DBWrite(metafile.get(), varname.c_str(), &value, &dims, 1, DB_INT);
        assert(!ierr);
      }

      // Write FabArrayBase (component positions and shapes)
      for (const auto &leveldata : ghext->leveldata) {
        const amrex::FabArrayBase &fab = *leveldata.fab;
        const int nfabs = fab.size();
        vector<int> boxes(2 * ndims * nfabs);
        for (int component = 0; component < nfabs; ++component) {
          const amrex::Box &fabbox = fab.box(component); // valid region
          for (int d = 0; d < ndims; ++d)
            boxes[d + 2 * ndims * component] = fabbox.smallEnd(d);
          for (int d = 0; d < ndims; ++d)
            boxes[d + ndims + 2 * ndims * component] = fabbox.bigEnd(d);
        }
        const int dims[2] = {nfabs, 2 * ndims};
        const string varname =
            dirname + "/" + make_fabarraybasename(leveldata.level);
        ierr = DBWrite(metafile.get(), varname.c_str(), boxes.data(), dims, 2,
                       DB_INT);
        assert(!ierr);
      }
    }

    {
      const string visitname = [&]() {
        ostringstream buf;
        buf << out_dir << "/" << simulation_name << ".visit";
        return buf.str();
      }();
      ofstream visit(visitname, ios::app);
      assert(visit.good());
      visit << make_filename(simulation_name, cctk_iteration) << "\n";
    }

  } // if write metadata

  interval_meta = nullptr;
}

} // namespace CarpetX

#endif
