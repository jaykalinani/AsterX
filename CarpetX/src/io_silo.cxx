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

#include <array>
#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <set>
#include <sstream>
#include <tuple>
#include <vector>

#ifdef HAVE_CAPABILITY_Silo
namespace CarpetX {
using namespace amrex;
using namespace std;

struct mesh_props_t {
  IntVect ngrow;

  auto to_tuple() const { return make_tuple(ngrow); }
  friend bool operator==(const mesh_props_t &p, const mesh_props_t &q) {
    return p.to_tuple() == q.to_tuple();
  }
  friend bool operator<(const mesh_props_t &p, const mesh_props_t &q) {
    return p.to_tuple() < q.to_tuple();
  }
};

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
  }
  ostringstream buf;
  buf << varname;
  if (reflevel >= 0)
    buf << ".rl" << setw(2) << setfill('0') << reflevel //
        << ".c" << setw(8) << setfill('0') << component;
  return DB::legalize_name(buf.str());
}

void OutputSilo(const cGH *restrict const cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int ierr;

  // Set up timers
  static Timer timer("OutputSilo");
  Interval interval(timer);

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

  // Create output file
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
  const vector<bool> group_enabled = [&]() {
    vector<bool> enabled(CCTK_NumGroups(), false);
    const auto callback{
        [](const int index, const char *const optstring, void *const arg) {
          vector<bool> &enabled = *static_cast<vector<bool> *>(arg);
          enabled.at(CCTK_GroupIndexFromVarI(index)) = true;
        }};
    CCTK_TraverseString(out_silo_vars, callback, &enabled, CCTK_GROUP_OR_VAR);
    return enabled;
  }();

  constexpr int ndims = dim;

  // Write data
  {
    DB::ptr<DBfile> file;
    if (write_file) {
      const string filename =
          string(out_dir) + "/" +
          make_filename(simulation_name, cctk_iteration, myproc / ioproc_every);
      file = DB::make(DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL,
                               simulation_name.c_str(), DB_HDF5));
      assert(file);
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
        const MultiFab &mfab = *groupdata.mfab[tl];
        const IndexType &indextype = mfab.ixType();
        const IntVect &ngrow = mfab.nGrowVect();
        const DistributionMapping &dm = mfab.DistributionMap();

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

          const Box &fabbox = mfab.fabbox(component); // exterior

          array<int, ndims> dims;
          for (int d = 0; d < ndims; ++d)
            dims[d] = fabbox.length(d);
          ptrdiff_t zonecount = 1;
          for (int d = 0; d < ndims; ++d)
            zonecount *= dims[d];
          assert(zonecount >= 0 && zonecount <= INT_MAX);

          if (write_file && !have_mesh) {
            const string meshname = make_meshname(leveldata.level, component);

            array<int, ndims> dims_vc;
            for (int d = 0; d < ndims; ++d)
              dims_vc[d] = dims[d] + int(indextype.cellCentered(d));

            const Geometry &geom = ghext->amrcore->Geom(leveldata.level);
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
          const int mpi_tag = 22900; // randomly chosen
          vector<double> buffer;
          const double *data = nullptr;
          if (send_this_fab && write_this_fab) {
            const FArrayBox &fab = mfab[component];
            data = fab.dataPtr();
          } else if (send_this_fab) {
            const FArrayBox &fab = mfab[component];
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

          // Write variable
          if (write_file) {
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

  // Write metadata
  if (write_metafile) {

    const string metafilename =
        string(out_dir) + "/" + make_filename(simulation_name, cctk_iteration);
    const DB::ptr<DBfile> metafile =
        DB::make(DBCreate(metafilename.c_str(), DB_CLOBBER, DB_LOCAL,
                          simulation_name.c_str(), DB_HDF5));
    assert(metafile);

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
      const MultiFab &mfab0 = *groupdata0.mfab[tl];
      const IndexType &indextype = mfab0.ixType();
      const IntVect &ngrow = mfab0.nGrowVect();

      const mesh_props_t mesh_props{ngrow};
      const bool have_mesh = have_meshes.count(mesh_props);

      if (!have_mesh) {
        const string multimeshname = make_meshname();

        vector<string> meshnames;
        for (const auto &leveldata : ghext->leveldata) {
          const auto &groupdata = *leveldata.groupdata.at(gi);
          const MultiFab &mfab = *groupdata.mfab[tl];
          const DistributionMapping &dm = mfab.DistributionMap();
          const int nfabs = dm.size();
          for (int c = 0; c < nfabs; ++c) {
            const int proc = dm[c];
            const string proc_filename = make_filename(
                simulation_name, cctk_iteration, proc / ioproc_every);
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
        typedef array<array<double, ndims>, 2> extent_t;
        vector<extent_t> extents;
        extents.reserve(meshnames.size());
        for (const auto &leveldata : ghext->leveldata) {
          const auto &groupdata = *leveldata.groupdata.at(gi);
          const int tl = 0;
          const MultiFab &mfab = *groupdata.mfab[tl];
          const Geometry &geom = ghext->amrcore->Geom(leveldata.level);
          const double *const x0 = geom.ProbLo();
          const double *const dx = geom.CellSize();
          const int nfabs = mfab.size();
          for (int c = 0; c < nfabs; ++c) {
            const Box &fabbox = mfab.fabbox(c); // exterior
            extent_t extent;
            for (int d = 0; d < ndims; ++d) {
              extent[0][d] = x0[d] + fabbox.smallEnd(d) * dx[d];
              extent[1][d] = x0[d] + fabbox.bigEnd(d) * dx[d];
            }
            extents.push_back(extent);
          }
        }
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
          const MultiFab &mfab = *groupdata.mfab[tl];
          const int nfabs = mfab.size();
          for (int c = 0; c < nfabs; ++c) {
            const Box &fabbox = mfab.fabbox(c); // exterior
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
            const MultiFab &mfab = *groupdata.mfab[tl];
            const DistributionMapping &dm = mfab.DistributionMap();
            const int nfabs = dm.size();
            for (int c = 0; c < nfabs; ++c) {
              const int proc = dm[c];
              const string proc_filename = make_filename(
                  simulation_name, cctk_iteration, proc / ioproc_every);
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
  }   // if write metadata
}

} // namespace CarpetX
#endif
