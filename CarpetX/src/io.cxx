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

namespace CarpetX {
using namespace amrex;
using namespace std;

void OutputPlotfile(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("OutputPlotfile");
  Interval interval(timer);

  const int numgroups = CCTK_NumGroups();
  vector<bool> group_enabled(numgroups, false);
  auto enable_group{[](int index, const char *optstring, void *callback) {
    vector<bool> &group_enabled = *static_cast<vector<bool> *>(callback);
    group_enabled.at(index) = true;
  }};
  CCTK_TraverseString(out_plotfile_groups, enable_group, &group_enabled,
                      CCTK_GROUP);

  for (int gi = 0; gi < numgroups; ++gi) {
    if (!group_enabled.at(gi))
      continue;
    if (CCTK_GroupTypeI(gi) != CCTK_GF)
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
    const auto &groupdata = leveldata.groupdata.at(gi);
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

    file << "# 1:iteration" << sep << "2:time" << sep << "3:varname" << sep
         << "4:min" << sep << "5:max" << sep << "6:sum" << sep << "7:avg" << sep
         << "8:stddev" << sep << "9:volume" << sep << "10:L1norm" << sep
         << "11:L2norm" << sep << "12:maxabs"
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
        file << cctk_iteration << sep << cctk_time << sep
             << CCTK_FullVarName(groupdata.firstvarindex + vi) << sep << red.min
             << sep << red.max << sep << red.sum << sep << red.avg() << sep
             << red.sdv() << sep << red.norm0() << sep << red.norm1() << sep
             << red.norm2() << sep << red.norm_inf() << "\n";
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
  const string sep = "\t";

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
        globaldata.scalargroupdata.at(gi);
    const int tl = 0;
    file << cctkGH->cctk_iteration << sep << cctkGH->cctk_time;
    for (int vi = 0; vi < scalargroupdata.numvars; ++vi) {
      file << sep << scalargroupdata.data.at(tl).at(vi);
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

} // namespace CarpetX
