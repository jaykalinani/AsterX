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

#include <cctype>
#include <fstream>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <string>
#include <type_traits>
#include <utility>

namespace AMReX {
using namespace amrex;
using namespace std;

void OutputPlotfile(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("OutputPlotfile");
  Interval interval(timer);

  const int numgroups = CCTK_NumGroups();
  for (int gi = 0; gi < numgroups; ++gi) {
    cGroup group;
    int ierr = CCTK_GroupData(gi, &group);
    assert(!ierr);

    if (group.grouptype != CCTK_GF)
      continue;

    auto &restrict groupdata0 = ghext->leveldata.at(0).groupdata.at(gi);
    if (groupdata0.mfab.size() > 0) {
      const int tl = 0;

      string groupname = unique_ptr<char>(CCTK_GroupName(gi)).get();
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
        mfabs.at(leveldata.level) = &*leveldata.groupdata.at(gi).mfab.at(tl);
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

#warning "TODO: Get these from schedule.cxx"
struct TileBox {
  array<int, dim> tile_min;
  array<int, dim> tile_max;
};

void enter_level_mode(cGH *restrict cctkGH,
                      const GHExt::LevelData &restrict leveldata);
void leave_level_mode(cGH *restrict cctkGH,
                      const GHExt::LevelData &restrict leveldata);
void enter_local_mode(cGH *restrict cctkGH, TileBox &restrict tilebox,
                      const GHExt::LevelData &restrict leveldata,
                      const MFIter &mfi);
void leave_local_mode(cGH *restrict cctkGH, TileBox &restrict tilebox,
                      const GHExt::LevelData &restrict leveldata,
                      const MFIter &mfi);

void WriteASCII(const cGH *restrict cctkGH, const string &filename, int gi,
                const vector<string> &varnames) {
  ostringstream buf;
  buf << filename << ".tsv";
  ofstream file(buf.str());

  // Output header
  file << "# 1:iteration"
       << "\t"
       << "2:time"
       << "\t"
       << "3:level"
       << "\t"
       << "4:component"
       << "\t"
       << "5:i"
       << "\t"
       << "6:j"
       << "\t"
       << "7:k"
       << "\t"
       << "8:x"
       << "\t"
       << "9:y"
       << "\t"
       << "10:z";
  int col = 11;
  for (const auto &varname : varnames)
    file << "\t" << col++ << ":" << varname;
  file << "\n";

  // get more precision for floats, could also use
  // https://stackoverflow.com/a/30968371
  file << setprecision(numeric_limits<CCTK_REAL>::digits10 + 1) << scientific;

  for (const auto &leveldata : ghext->leveldata) {
    enter_level_mode(const_cast<cGH *>(cctkGH), leveldata);
    const auto &groupdata = leveldata.groupdata.at(gi);
    const int tl = 0;
    for (MFIter mfi(*leveldata.mfab0); mfi.isValid(); ++mfi) {
      TileBox tilebox;
      enter_local_mode(const_cast<cGH *>(cctkGH), tilebox, leveldata, mfi);
      GridPtrDesc1 grid(leveldata, groupdata, mfi);
      const Array4<const CCTK_REAL> &vars = groupdata.mfab.at(tl)->array(mfi);
      vector<GF3D1<const CCTK_REAL> > ptrs_;
      ptrs_.reserve(groupdata.numvars);
      for (int vi = 0; vi < groupdata.numvars; ++vi)
        ptrs_.push_back(grid.gf3d(vars, vi));
      // write_arrays(file, cctkGH, leveldata.level, mfi.index(), ptrs, grid);
      Loop::loop_idx(cctkGH, where_t::everywhere, groupdata.indextype,
                     [&](const Loop::PointDesc &p) {
                       file << cctkGH->cctk_iteration << "\t"
                            << cctkGH->cctk_time << "\t" << leveldata.level
                            << "\t"
                            << mfi.index()
                            // << "\t" << isghost
                            << "\t" << (grid.lbnd[0] + p.i) << "\t"
                            << (grid.lbnd[1] + p.j) << "\t"
                            << (grid.lbnd[2] + p.k) << "\t" << p.x << "\t"
                            << p.y << "\t" << p.z;
                       for (const auto &ptr_ : ptrs_)
                         file << "\t" << ptr_(p.I);
                       file << "\n";
                     });
      leave_local_mode(const_cast<cGH *>(cctkGH), tilebox, leveldata, mfi);
    }
    leave_level_mode(const_cast<cGH *>(cctkGH), leveldata);
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

    auto &restrict groupdata0 = ghext->leveldata.at(0).groupdata.at(gi);
    if (groupdata0.mfab.size() > 0) {
      const int tl = 0;

      string groupname = unique_ptr<char>(CCTK_GroupName(gi)).get();
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

void OutputNorms(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("OutputNorms");
  Interval interval(timer);

  const bool is_root = CCTK_MyProc(nullptr) == 0;

  ofstream file;
  if (is_root) {
    ostringstream buf;
    buf << out_dir << "/norms.it" << setw(6) << setfill('0') << cctk_iteration
        << ".tsv";
    const string filename = buf.str();
    file = ofstream(filename);
    file << "# 1:iteration"
         << "\t"
         << "2:time"
         << "\t"
         << "3:varname"
         << "\t"
         << "4:min"
         << "\t"
         << "5:max"
         << "\t"
         << "6:sum"
         << "\t"
         << "7:avg"
         << "\t"
         << "8:stddev"
         << "\t"
         << "9:volume"
         << "\t"
         << "10:maxabs"
         << "\t"
         << "11:L1norm"
         << "\t"
         << "12:L2norm"
         << "\n";
  }

  const int numgroups = CCTK_NumGroups();
  for (int gi = 0; gi < numgroups; ++gi) {
    if (CCTK_GroupTypeI(gi) != CCTK_GF)
      continue;

    const int level = 0;
    const GHExt::LevelData &restrict leveldata = ghext->leveldata.at(level);
    const GHExt::LevelData::GroupData &restrict groupdata =
        leveldata.groupdata.at(gi);

    const int tl = 0;
    for (int vi = 0; vi < groupdata.numvars; ++vi) {

#if 1
      reduction<CCTK_REAL> red = reduce(gi, vi, tl);
#else
      reduction<CCTK_REAL> red;
      for (auto &restrict leveldata : ghext->leveldata) {
        auto &restrict groupdata = leveldata.groupdata.at(gi);
        MultiFab &mfab = *groupdata.mfab.at(tl);
        reduction<CCTK_REAL> red1;
        red1.min = mfab.min(vi);
        red1.max = mfab.max(vi);
        red1.sum = mfab.sum(vi);
        // red1.sum2 = mfab.sum2(vi);
        // red1.vol = mfab.vol(vi);
        red1.maxabs = mfab.norminf(vi);
        red1.sumabs = mfab.norm1(vi, mfab.fb_period);
        red1.sum2abs = pow(mfab.norm2(vi), 2);
        red += red1;
      }
#endif

      if (is_root)
        file << cctk_iteration << "\t" << cctk_time << "\t"
             << CCTK_FullVarName(groupdata.firstvarindex + vi) << "\t"
             << red.min << "\t" << red.max << "\t" << red.sum << "\t"
             << red.avg() << "\t" << red.sdv() << "\t" << red.norm0() << "\t"
             << red.norm1() << "\t" << red.norm2() << "\t" << red.norm_inf()
             << "\n";
    }
  }

  if (is_root)
    file.close();
}

void OutputScalars(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("OutputNorms");
  Interval interval(timer);

  const bool is_root = CCTK_MyProc(nullptr) == 0;
  if (!is_root)
    return;

  const int numgroups = CCTK_NumGroups();
  for (int gi = 0; gi < numgroups; ++gi) {
    cGroup group;
    int ierr = CCTK_GroupData(gi, &group);
    assert(!ierr);

    if (group.grouptype != CCTK_SCALAR)
      continue;

    string groupname = unique_ptr<char>(CCTK_GroupName(gi)).get();
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
         << "\t"
         << "iteration"
         << "\t"
         << "time";
    for (int vi = 0; vi < groupdata.numvars; ++vi) {
      file << "\t" << CCTK_VarName(group.firstvarindex + vi);
    file << "\n";
#endif

    // Output data
    const GHExt::GlobalData &restrict globaldata = ghext->globaldata;
    const GHExt::GlobalData::ScalarGroupData &restrict scalargroupdata =
        globaldata.scalargroupdata.at(gi);
    const int tl = 0;
    file << cctkGH->cctk_iteration << "\t" << cctkGH->cctk_time;
    for (int vi = 0; vi < scalargroupdata.numvars; ++vi) {
      file << "\t" << scalargroupdata.data.at(tl).at(vi);
    }
    file << "\n";
  }
}

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
    OutputASCII(cctkGH);
    OutputNorms(cctkGH);
    OutputScalars(cctkGH);
  }

  // TODO: This should be the number of variables output
  return 0;
}

} // namespace AMReX
