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
    for (const auto &restrict leveldata : ghext->leveldata) {
      auto &restrict groupdata = leveldata.groupdata.at(gi);
      if (groupdata.mfab.size() > 0) {
        const int tl = 0;

        string groupname = unique_ptr<char>(CCTK_GroupName(gi)).get();
        groupname = regex_replace(groupname, regex("::"), "-");
        for (auto &c : groupname)
          c = tolower(c);
        ostringstream buf;
        buf << "wavetoy/" << groupname;
        buf << ".rl" << setw(2) << setfill('0') << leveldata.level;
        buf << ".it" << setw(6) << setfill('0') << cctk_iteration;
        string filename = buf.str();

        Vector<string> varnames(groupdata.numvars);
        for (int vi = 0; vi < groupdata.numvars; ++vi) {
          ostringstream buf;
          buf << CCTK_VarName(groupdata.firstvarindex + vi);
          for (int i = 0; i < tl; ++i)
            buf << "_p";
          varnames.at(vi) = buf.str();
        }

        // TODO: Output all levels into a single file
        // TODO: Output all groups into a single file
        WriteSingleLevelPlotfile(filename, *groupdata.mfab.at(tl), varnames,
                                 ghext->amrmesh->Geom(leveldata.level),
                                 cctk_time, cctk_iteration);
      }
    }
    count_vars += ghext->leveldata.at(0).groupdata.at(gi).numvars;
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
    }
  }

  return count_vars;
}

} // namespace AMReX
