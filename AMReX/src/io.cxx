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
extern vector<TileBox> thread_local_tilebox;
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
  file << "#"
       << "\t"
       << "time"
       << "\t"
       << "level"
       << "\t"
       << "component"
       << "\t"
       << "i"
       << "\t"
       << "j"
       << "\t"
       << "k"
       << "\t"
       << "x"
       << "\t"
       << "y"
       << "\t"
       << "z";
  for (const auto &varname : varnames)
    file << "\t" << varname;
  file << "\n";

  for (const auto &leveldata : ghext->leveldata) {
    enter_level_mode(const_cast<cGH *>(cctkGH), leveldata);
    const auto &groupdata = leveldata.groupdata.at(gi);
    const int tl = 0;
    for (MFIter mfi(*leveldata.mfab0); mfi.isValid(); ++mfi) {
      TileBox &tilebox = thread_local_tilebox.at(omp_get_thread_num());
      enter_local_mode(const_cast<cGH *>(cctkGH), tilebox, leveldata, mfi);
      GridPtrDesc grid(leveldata, mfi);
      const Array4<CCTK_REAL> &vars = groupdata.mfab.at(tl)->array(mfi);
      vector<const CCTK_REAL *> ptrs(groupdata.numvars);
      for (int vi = 0; vi < groupdata.numvars; ++vi)
        ptrs.at(vi) = grid.ptr(vars, vi);
      // write_arrays(file, cctkGH, leveldata.level, mfi.index(), ptrs, grid);
      Loop::loop_all(
          cctkGH, groupdata.indextype, [&](const Loop::PointDesc &p) {
            file << cctkGH->cctk_iteration << "\t" << cctkGH->cctk_time << "\t"
                 << leveldata.level << "\t"
                 << mfi.index()
                 // << "\t" << isghost
                 << "\t" << (grid.lbnd[0] + p.i) << "\t" << (grid.lbnd[1] + p.j)
                 << "\t" << (grid.lbnd[2] + p.k) << "\t" << p.x << "\t" << p.y
                 << "\t" << p.z;
            for (const auto &ptr : ptrs)
              file << "\t" << ptr[p.idx];
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

  const int numgroups = CCTK_NumGroups();
  for (int gi = 0; gi < numgroups; ++gi) {
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
        cout << "  "
             << unique_ptr<char>(CCTK_FullName(groupdata.firstvarindex + vi))
                    .get()
             << ": maxabs=" << red.maxabs << " sum=" << red.sum
             << " vol=" << red.vol << "\n";
    }
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
  }

  // TODO: This should be the number of variables output
  return 0;
}

} // namespace AMReX
