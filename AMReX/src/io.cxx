#include <driver.hxx>
#include <io.hxx>
#include <reduction.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>

#include <cctype>
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

int OutputGH(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (out_every <= 0 || cctk_iteration % out_every != 0)
    return 0;

  int count_vars = 0;
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
      buf << "wavetoy/" << groupname;
      buf << ".it" << setw(6) << setfill('0') << cctk_iteration;
      string filename = buf.str();

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

      count_vars += ghext->leveldata.at(0).groupdata.at(gi).numvars;
    }
  }

  for (int gi = 0; gi < numgroups; ++gi) {
    const int level = 0;
    const GHExt::LevelData &restrict leveldata = ghext->leveldata.at(level);
    const GHExt::LevelData::GroupData &restrict groupdata =
        leveldata.groupdata.at(gi);
    const int numvars = groupdata.numvars;
    for (int vi = 0; vi < numvars; ++vi) {
      const int tl = 0;
      reduction<CCTK_REAL> red = reduce(gi, vi, tl);
      CCTK_VINFO("maxabs(%s)=%g", CCTK_VarName(groupdata.firstvarindex + vi),
                 red.maxabs);
      // CCTK_REAL maxabs = 0.0;
      // for (auto &restrict leveldata : ghext->leveldata) {
      //   auto &restrict groupdata = leveldata.groupdata.at(gi);
      //   MultiFab &mfab = *groupdata.mfab.at(tl);
      //   maxabs = fmax(maxabs, mfab.norminf(vi));
      // }
      // CCTK_VINFO("            %g", maxabs);
    }
  }

  return count_vars;
}

} // namespace AMReX
