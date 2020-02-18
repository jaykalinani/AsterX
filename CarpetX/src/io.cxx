#include "driver.hxx"
#include "io.hxx"
#include "reduction.hxx"
#include "schedule.hxx"
#include "timer.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <AMReX.H>
#include <AMReX_Orientation.H>
#include <AMReX_PlotFileUtil.H>

#ifdef HAVE_CAPABILITY_Silo
#include <silo.hxx>
#endif

#include <array>
#include <cctype>
#include <fstream>
#include <functional>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace CarpetX {
using namespace amrex;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

void OutputPlotfile(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("OutputPlotfile");
  Interval interval(timer);

  const int numgroups = CCTK_NumGroups();
  vector<bool> group_enabled(numgroups, false);
  auto enable_group{[](int index, const char *optstring, void *callback) {
    vector<bool> &group_enabled = *static_cast<vector<bool> *>(callback);
    // CCTK_TraverseString expands the given groups into their variables
    // so I condense it down again to groups only
    group_enabled.at(CCTK_GroupIndexFromVarI(index)) = true;
  }};
  CCTK_TraverseString(out_plotfile_groups, enable_group, &group_enabled,
                      CCTK_GROUP_OR_VAR);

  for (int gi = 0; gi < numgroups; ++gi) {
    if (!group_enabled.at(gi))
      continue;
    if (CCTK_GroupTypeI(gi) != CCTK_GF)
      continue;

    auto &restrict groupdata0 = *ghext->leveldata.at(0).groupdata.at(gi);
    if (groupdata0.mfab.size() > 0) {
      const int tl = 0;

      string groupname = unique_C_ptr<char>(CCTK_GroupName(gi)).get();
      groupname = regex_replace(groupname, regex("::"), "-");
      for (auto &c : groupname)
        c = tolower(c);
      const string filename = [&]() {
        ostringstream buf;
        buf << out_dir << "/" << groupname << ".it" << setw(6) << setfill('0')
            << cctk_iteration;
        return buf.str();
      }();

      Vector<string> varnames(groupdata0.numvars);
      for (int vi = 0; vi < groupdata0.numvars; ++vi) {
        ostringstream buf;
        buf << CCTK_VarName(groupdata0.firstvarindex + vi);
        for (int i = 0; i < tl; ++i)
          buf << "_p";
        varnames.at(vi) = buf.str();
      }

      Vector<const MultiFab *> mfabs(ghext->leveldata.size());
      Vector<Geometry> geoms(ghext->leveldata.size());
      Vector<int> iters(ghext->leveldata.size());
      Vector<IntVect> reffacts(ghext->leveldata.size());
      for (const auto &restrict leveldata : ghext->leveldata) {
        mfabs.at(leveldata.level) = &*leveldata.groupdata.at(gi)->mfab.at(tl);
        geoms.at(leveldata.level) = ghext->amrcore->Geom(leveldata.level);
        iters.at(leveldata.level) = cctk_iteration;
        reffacts.at(leveldata.level) = IntVect{2, 2, 2};
      }

      // TODO: Output all groups into a single file
      WriteMultiLevelPlotfile(filename, mfabs.size(), mfabs, varnames, geoms,
                              cctk_time, iters, reffacts);

      const bool is_root = CCTK_MyProc(nullptr) == 0;
      if (is_root) {
        const string visitname = [&]() {
          ostringstream buf;
          buf << out_dir << "/" << groupname << ".visit";
          return buf.str();
        }();
        ofstream visit(visitname, ios::app);
        assert(visit.good());
        // visit << filename << "/Header\n";
        visit << groupname << ".it" << setw(6) << setfill('0') << cctk_iteration
              << "/Header\n";
        visit.close();
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

#ifdef HAVE_CAPABILITY_Silo

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

void OutputSilo(const cGH *restrict const cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int ierr;

  // Set up timers
  static Timer timer("OutputSilo");
  Interval interval(timer);

  const int myproc = CCTK_MyProc(cctkGH);
  const int nprocs = CCTK_nProcs(cctkGH);

  // Notes for later
  out_mode;
  out_proc_every;

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

  const string name = [&]() {
    string name = parfilename;
    const size_t last_slash = name.rfind('/');
    if (last_slash != string::npos && last_slash < name.length())
      name = name.substr(last_slash + 1);
    const size_t last_dot = name.rfind('.');
    if (last_dot != string::npos && last_dot > 0)
      name = name.substr(0, last_dot);
    return name;
  }();

  const string filename = [&]() {
    ostringstream buf;
    buf << out_dir << "/" << name << ".it" << setw(6) << setfill('0')
        << cctk_iteration << ".p" << setw(6) << setfill('0') << myproc
        << ".silo";
    return buf.str();
  }();
  const DB::ptr<DBfile> file = DB::make(
      DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL, name.c_str(), DB_HDF5));
  assert(file);

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

  // Iterate over all groups that should be output
  map<mesh_props_t, string> all_multimeshnames;
  for (int gi = 0; gi < CCTK_NumGroups(); ++gi) {
    if (!group_enabled.at(gi))
      continue;
    if (CCTK_GroupTypeI(gi) != CCTK_GF)
      continue;

    // Write group
    const int v0 = CCTK_FirstVarIndexI(gi);
    const int nv = CCTK_NumVarsInGroupI(gi);
    for (int vi = 0; vi < nv; ++vi) {

      const auto &leveldata0 = *ghext->leveldata.begin();
      // const auto &groupdata0 = *leveldata0.groupdata.at(gi);
      const int tl = 0;
      const MultiFab &mfab0 = *leveldata0.groupdata.at(gi)->mfab[tl];
      const IndexType &indextype = mfab0.ixType();
      const IntVect &ngrow = mfab0.nGrowVect();
      const mesh_props_t mesh_props{ngrow};
      const bool have_multimesh = all_multimeshnames.count(mesh_props);

      constexpr int ndims = dim;

      vector<string> meshnames;
      vector<array<array<double, ndims>, 2> > meshextents;
      vector<int> meshzonecounts;
      vector<string> varnames;
      for (const auto &leveldata : ghext->leveldata) {
        const auto &groupdata = *leveldata.groupdata.at(gi);
        const int tl = 0;
        const MultiFab &mfab = *groupdata.mfab[tl];
        assert(mfab.ixType() == indextype);
        assert(mfab.nGrowVect() == ngrow);

        for (MFIter mfi(mfab); mfi.isValid(); ++mfi) {
          // const Box &box = mfi.validbox();  // interior
          const Box &fabbox = mfi.fabbox(); // exterior
          const FArrayBox &fab = mfab[mfi];
          assert(fab.box() == fabbox);

          const string meshname = [&]() {
            ostringstream buf;
            buf << "box.rl" << setw(2) << setfill('0') << leveldata.level
                << ".c" << setw(8) << setfill('0') << mfi.index();
            return DB::legalize_name(buf.str());
          }();

          if (!have_multimesh) {
            meshnames.push_back(meshname);

            array<int, ndims> dims_vc;
            for (int d = 0; d < ndims; ++d)
              dims_vc[d] = fabbox.length(d) + int(indextype.cellCentered(d));

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

            array<array<double, ndims>, 2> extent;
            for (int d = 0; d < ndims; ++d) {
              extent[0][d] = *coords[d].begin();
              extent[1][d] = *coords[d].rbegin();
            }
            meshextents.push_back(extent);

            int zonecount = 1;
            for (int d = 0; d < ndims; ++d)
              zonecount *= dims_vc[d];
            meshzonecounts.push_back(zonecount);

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

            int column_major = 1;
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
          }

          const string varname = [&]() {
            ostringstream buf;
            buf << CCTK_FullVarName(v0 + vi) << ".rl" << setw(2) << setfill('0')
                << leveldata.level << ".c" << setw(8) << setfill('0')
                << mfi.index();
            return DB::legalize_name(buf.str());
          }();
          varnames.push_back(varname);

          array<int, ndims> dims;
          for (int d = 0; d < ndims; ++d)
            dims[d] = fabbox.length(d);

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

          int column_major = 1;
          ierr = DBAddOption(optlist.get(), DBOPT_MAJORORDER, &column_major);
          assert(!ierr);

          // float time = cctk_time;
          // ierr = DBAddOption(optlist.get(), DBOPT_TIME, &time);
          // assert(!ierr);
          double dtime = cctk_time;
          ierr = DBAddOption(optlist.get(), DBOPT_DTIME, &dtime);
          assert(!ierr);

          // int hide_from_gui = 1;
          // ierr =
          //     DBAddOption(optlist.get(), DBOPT_HIDE_FROM_GUI,
          //     &hide_from_gui);
          // assert(!ierr);

          ierr = DBPutQuadvar1(file.get(), varname.c_str(), meshname.c_str(),
                               fab.dataPtr(vi), dims.data(), ndims, nullptr, 0,
                               DB_DOUBLE, centering, optlist.get());
          assert(!ierr);

        } // box
      }   // level

      const string multimeshname = [&]() {
        ostringstream buf;
        buf << "multimesh";
        return DB::legalize_name(buf.str());
      }();

      if (!have_multimesh) {
        all_multimeshnames[mesh_props] = multimeshname;

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
        assert(sizeof(*meshextents.data()) == extents_size * sizeof(double));
        ierr = DBAddOption(optlist.get(), DBOPT_EXTENTS_SIZE, &extents_size);
        assert(!ierr);
        assert(meshextents.size() == meshname_ptrs.size());
        ierr = DBAddOption(optlist.get(), DBOPT_EXTENTS, meshextents.data());
        assert(!ierr);

        assert(meshzonecounts.size() == meshname_ptrs.size());
        ierr =
            DBAddOption(optlist.get(), DBOPT_ZONECOUNTS, meshzonecounts.data());
        assert(!ierr);

        ierr = DBPutMultimesh(file.get(), multimeshname.c_str(),
                              meshname_ptrs.size(), meshname_ptrs.data(),
                              nullptr, optlist.get());
        assert(!ierr);
      }

      const string multivarname = [&]() {
        ostringstream buf;
        buf << CCTK_FullVarName(v0 + vi);
        return DB::legalize_name(buf.str());
      }();

      vector<const char *> varname_ptrs;
      varname_ptrs.reserve(varnames.size());
      for (const auto &varname : varnames)
        varname_ptrs.push_back(varname.c_str());

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

      ierr =
          DBPutMultivar(file.get(), multivarname.c_str(), varname_ptrs.size(),
                        varname_ptrs.data(), nullptr, optlist.get());
      assert(!ierr);

    } // vi
  }   // gi
}

#endif

////////////////////////////////////////////////////////////////////////////////

void WriteASCII(const cGH *restrict cctkGH, const string &filename, int gi,
                const vector<string> &varnames) {
  ostringstream buf;
  buf << filename << ".tsv";
  ofstream file(buf.str());
  const string sep = "\t";

  // get more precision for floats, could also use
  // https://stackoverflow.com/a/30968371
  file << setprecision(numeric_limits<CCTK_REAL>::digits10 + 1) << scientific;

  // Output header
  file << "# 1:iteration" << sep << "2:time" << sep << "3:level" << sep
       << "4:component" << sep << "5:i" << sep << "6:j" << sep << "7:k" << sep
       << "8:x" << sep << "9:y" << sep << "10:z";
  int col = 11;
  for (const auto &varname : varnames)
    file << sep << col++ << ":" << varname;
  file << "\n";

  for (const auto &leveldata : ghext->leveldata) {
    const auto &groupdata = *leveldata.groupdata.at(gi);
    const int tl = 0;
    const auto &geom = ghext->amrcore->Geom(leveldata.level);
    const auto &mfab = *groupdata.mfab.at(tl);
    for (MFIter mfi(mfab); mfi.isValid(); ++mfi) {
      const Array4<const CCTK_REAL> &vars = mfab.array(mfi);
      const auto &imin = vars.begin;
      const auto &imax = vars.end;
      for (int k = imin.z; k < imax.z; ++k) {
        for (int j = imin.y; j < imax.y; ++j) {
          for (int i = imin.x; i < imax.x; ++i) {
            const array<int, dim> I{i, j, k};
            array<CCTK_REAL, dim> x;
            for (int d = 0; d < dim; ++d)
              x[d] = geom.ProbLo(d) +
                     (I[d] + 0.5 * groupdata.indextype[d]) * geom.CellSize(d);
            file << cctkGH->cctk_iteration << sep << cctkGH->cctk_time << sep
                 << leveldata.level << sep << mfi.index() << sep << I[0] << sep
                 << I[1] << sep << I[2] << sep << x[0] << sep << x[1] << sep
                 << x[2];
            for (int n = 0; n < groupdata.numvars; ++n)
              file << sep << vars(i, j, k, n);
            file << "\n";
          }
        }
      }
    }
  }

  file.close();
}

void OutputASCII(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (!out_tsv)
    return;

  static Timer timer("OutputASCII");
  Interval interval(timer);

  const int numgroups = CCTK_NumGroups();
  for (int gi = 0; gi < numgroups; ++gi) {
    cGroup group;
    int ierr = CCTK_GroupData(gi, &group);
    assert(!ierr);

    if (group.grouptype != CCTK_GF)
      continue;

    auto &restrict groupdata0 = *ghext->leveldata.at(0).groupdata.at(gi);
    if (groupdata0.mfab.size() > 0) {
      const int tl = 0;

      string groupname = unique_C_ptr<char>(CCTK_GroupName(gi)).get();
      groupname = regex_replace(groupname, regex("::"), "-");
      for (auto &c : groupname)
        c = tolower(c);
      ostringstream buf;
      buf << out_dir << "/" << groupname;
      buf << ".it" << setw(6) << setfill('0') << cctk_iteration;
      buf << ".p" << setw(4) << setfill('0') << CCTK_MyProc(nullptr);
      const string filename = buf.str();

      Vector<string> varnames(groupdata0.numvars);
      for (int vi = 0; vi < groupdata0.numvars; ++vi) {
        ostringstream buf;
        buf << CCTK_VarName(groupdata0.firstvarindex + vi);
        for (int i = 0; i < tl; ++i)
          buf << "_p";
        varnames.at(vi) = buf.str();
      }

      WriteASCII(cctkGH, filename, gi, varnames);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void OutputNorms(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("OutputNorms");
  Interval interval(timer);

  const bool is_root = CCTK_MyProc(nullptr) == 0;
  const string sep = "\t";

  ofstream file;
  if (is_root) {
    ostringstream buf;
    buf << out_dir << "/norms.it" << setw(6) << setfill('0') << cctk_iteration
        << ".tsv";
    const string filename = buf.str();
    file = ofstream(filename);

    // get more precision for floats, could also use
    // https://stackoverflow.com/a/30968371
    file << setprecision(numeric_limits<CCTK_REAL>::digits10 + 1) << scientific;

    int col = 0;
    file << "# " << ++col << ":iteration";
    file << sep << ++col << ":time";
    file << sep << ++col << ":varname";
    file << sep << ++col << ":min";
    file << sep << ++col << ":max";
    file << sep << ++col << ":sum";
    file << sep << ++col << ":avg";
    file << sep << ++col << ":stddev";
    file << sep << ++col << ":volume";
    file << sep << ++col << ":L1norm";
    file << sep << ++col << ":L2norm";
    file << sep << ++col << ":maxabs";
    for (int d = 0; d < dim; ++d)
      file << sep << ++col << ":minloc[" << d << "]";
    for (int d = 0; d < dim; ++d)
      file << sep << ++col << ":maxloc[" << d << "]";
    file << "\n";
  }

  const int numgroups = CCTK_NumGroups();
  for (int gi = 0; gi < numgroups; ++gi) {
    if (CCTK_GroupTypeI(gi) != CCTK_GF)
      continue;

    const int level = 0;
    const GHExt::LevelData &restrict leveldata = ghext->leveldata.at(level);
    const GHExt::LevelData::GroupData &restrict groupdata =
        *leveldata.groupdata.at(gi);

    const int tl = 0;
    for (int vi = 0; vi < groupdata.numvars; ++vi) {

      const reduction<CCTK_REAL, dim> red = reduce(gi, vi, tl);

      if (is_root) {
        file << cctk_iteration << sep << cctk_time << sep
             << CCTK_FullVarName(groupdata.firstvarindex + vi) << sep << red.min
             << sep << red.max << sep << red.sum << sep << red.avg() << sep
             << red.sdv() << sep << red.norm0() << sep << red.norm1() << sep
             << red.norm2() << sep << red.norm_inf();
        for (int d = 0; d < dim; ++d)
          file << sep << red.minloc[d];
        for (int d = 0; d < dim; ++d)
          file << sep << red.maxloc[d];
        file << "\n";
      }
    }
  }

  if (is_root)
    file.close();
}

////////////////////////////////////////////////////////////////////////////////

void OutputScalars(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("OutputNorms");
  Interval interval(timer);

  const bool is_root = CCTK_MyProc(nullptr) == 0;
  if (!is_root)
    return;
  const string sep = "\t";

  const int numgroups = CCTK_NumGroups();
  for (int gi = 0; gi < numgroups; ++gi) {
    cGroup group;
    int ierr = CCTK_GroupData(gi, &group);
    assert(!ierr);

    if (group.grouptype != CCTK_SCALAR)
      continue;

    string groupname = unique_C_ptr<char>(CCTK_GroupName(gi)).get();
    groupname = regex_replace(groupname, regex("::"), "-");
    for (auto &c : groupname)
      c = tolower(c);
    ostringstream buf;
    buf << out_dir << "/" << groupname << ".tsv";
    const string filename = buf.str();

    ofstream file(filename.c_str(), std::ofstream::out | std::ofstream::app);

    // TODO: write proper header (once)
#if 0
    // Output header
    file << "#"
         << sep
         << "iteration"
         << sep
         << "time";
    for (int vi = 0; vi < groupdata.numvars; ++vi) {
      file << sep << CCTK_VarName(group.firstvarindex + vi);
    file << "\n";
#endif

    // Output data
    const GHExt::GlobalData &restrict globaldata = ghext->globaldata;
    const GHExt::GlobalData::ScalarGroupData &restrict scalargroupdata =
        *globaldata.scalargroupdata.at(gi);
    const int tl = 0;
    file << cctkGH->cctk_iteration << sep << cctkGH->cctk_time;
    for (int vi = 0; vi < scalargroupdata.numvars; ++vi) {
      file << sep << scalargroupdata.data.at(tl).at(vi);
    }
    file << "\n";
  }
}

////////////////////////////////////////////////////////////////////////////////

int OutputGH(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("OutputGH");
  Interval interval(timer);

  const bool is_root = CCTK_MyProc(nullptr) == 0;
  if (is_root)
    cout << "OutputGH: iteration " << cctk_iteration << ", time " << cctk_time
         << ", run time " << CCTK_RunTime() << " s\n";

  if (out_every > 0 && cctk_iteration % out_every == 0) {
    OutputPlotfile(cctkGH);
#ifdef HAVE_CAPABILITY_Silo
    // TODO: Stop at paramcheck time when Silo output parameters are
    // set, but Silo is not available
    OutputSilo(cctkGH);
#endif
    OutputASCII(cctkGH);
    OutputNorms(cctkGH);
    OutputScalars(cctkGH);
  }

  // TODO: This should be the number of variables output
  return 0;
}

} // namespace CarpetX
