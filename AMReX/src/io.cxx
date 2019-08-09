#include "driver.hxx"
#include "io.hxx"
#include "reduction.hxx"
#include "schedule.hxx"

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
    }
  }
}

template <typename T>
void write_arrays(ostream &os, const cGH *restrict cctkGH, int level,
                  int component, const vector<const T *> ptrs,
                  const GridDesc &grid) {
  DECLARE_CCTK_ARGUMENTS;

  const int levfac = 1 << level;
  CCTK_REAL x0[dim], dx[dim];
  for (int d = 0; d < dim; ++d) {
    dx[d] = cctk_delta_space[d] / levfac;
    x0[d] = cctk_origin_space[d] + 0.5 * dx[d];
  }

  int imin[dim], imax[dim];
  for (int d = 0; d < dim; ++d) {
    imin[d] = grid.bbox[2 * d] ? 0 : grid.nghostzones[d];
    imax[d] = grid.lsh[d] - (grid.bbox[2 * d + 1] ? 0 : grid.nghostzones[d]);
  }

  constexpr int di = 1;
  const int dj = di * grid.ash[0];
  const int dk = dj * grid.ash[1];

  for (int k = 0; k < grid.lsh[2]; ++k) {
    for (int j = 0; j < grid.lsh[1]; ++j) {
      for (int i = 0; i < grid.lsh[0]; ++i) {
        const int idx = di * i + dj * j + dk * k;
        const int isghost = !(i >= imin[0] && i < imax[0] && j >= imin[1] &&
                              j < imax[1] && k >= imin[2] && k < imax[2]);
        T x = x0[0] + (grid.lbnd[0] + i) * dx[0];
        T y = x0[1] + (grid.lbnd[1] + j) * dx[1];
        T z = x0[2] + (grid.lbnd[2] + k) * dx[2];
        os << cctk_iteration << "\t" << cctk_time << "\t" << level << "\t"
           << component << "\t" << isghost << "\t" << (grid.lbnd[0] + i) << "\t"
           << (grid.lbnd[1] + j) << "\t" << (grid.lbnd[2] + k) << "\t" << x
           << "\t" << y << "\t" << z;
        for (const auto &ptr : ptrs)
          os << "\t" << ptr[idx];
        os << "\n";
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
       // << "\t"
       // << "isghost"
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
      string filename = buf.str();

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

  const bool is_root = CCTK_MyProc(nullptr) == 0;

  const int numgroups = CCTK_NumGroups();
  for (int gi = 0; gi < numgroups; ++gi) {
    const int level = 0;
    const GHExt::LevelData &restrict leveldata = ghext->leveldata.at(level);
    const GHExt::LevelData::GroupData &restrict groupdata =
        leveldata.groupdata.at(gi);
    const int numvars = groupdata.numvars;
    for (int vi = 0; vi < numvars; ++vi) {
      const int tl = 0;
      reduction<CCTK_REAL> red = reduce(gi, vi, tl);
      if (is_root)
        cout << "  "
             << unique_ptr<char>(CCTK_FullName(groupdata.firstvarindex + vi))
                    .get()
             << ": maxabs=" << red.maxabs << "\n";
      // CCTK_REAL maxabs = 0.0;
      // for (auto &restrict leveldata : ghext->leveldata) {
      //   auto &restrict groupdata = leveldata.groupdata.at(gi);
      //   MultiFab &mfab = *groupdata.mfab.at(tl);
      //   maxabs = fmax(maxabs, mfab.norminf(vi));
      // }
      // CCTK_VINFO("            %g", maxabs);
    }
  }
}

int OutputGH(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const bool is_root = CCTK_MyProc(nullptr) == 0;
  if (is_root)
    cout << "OutputGH: iteration " << cctk_iteration << ", time " << cctk_time
         << ", run time " << CCTK_RunTime() << " s\n";

  if (out_every <= 0 || cctk_iteration % out_every != 0)
    return 0;

  OutputPlotfile(cctkGH);
  OutputASCII(cctkGH);
  OutputNorms(cctkGH);

  // TODO: This should be the number of variables output
  return 0;
}

} // namespace AMReX
