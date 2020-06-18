#include "io_tsv.hxx"

#include "driver.hxx"
#include "timer.hxx"

#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <fstream>
#include <iomanip>
#include <limits>
#include <regex>
#include <string>
#include <vector>

namespace CarpetX {
using namespace std;

void WriteTSV(const cGH *restrict cctkGH, const string &filename, int gi,
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
    for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
      const amrex::Array4<const CCTK_REAL> &vars = mfab.array(mfi);
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

void OutputTSV(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (!out_tsv)
    return;

  static Timer timer("OutputTSV");
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

      amrex::Vector<string> varnames(groupdata0.numvars);
      for (int vi = 0; vi < groupdata0.numvars; ++vi) {
        ostringstream buf;
        buf << CCTK_VarName(groupdata0.firstvarindex + vi);
        for (int i = 0; i < tl; ++i)
          buf << "_p";
        varnames.at(vi) = buf.str();
      }

      WriteTSV(cctkGH, filename, gi, varnames);
    }
  }
}

} // namespace CarpetX
