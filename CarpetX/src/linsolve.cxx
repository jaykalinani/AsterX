#include "driver.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <AMReX.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>

#include <vector>

namespace CarpetX {
using namespace amrex;
using namespace std;

extern "C" void CarpetX_SolvePoisson() {
  // Create operator

  Vector<Geometry> geoms(ghext->leveldata.size());
  Vector<BoxArray> grids(ghext->leveldata.size());
  Vector<DistributionMapping> dmaps(ghext->leveldata.size());
  for (int level = 0; level < int(ghext->leveldata.size()); ++level) {
    geoms.at(level) = ghext->amrcore->Geom(level);
    grids.at(level) = ghext->amrcore->boxArray(level);
    dmaps.at(level) = ghext->amrcore->DistributionMap(level);
  }

  MLNodeLaplacian mlnodelaplacian(geoms, grids, dmaps);

  mlnodelaplacian.setDomainBC(
      {LinOpBCType::Periodic, LinOpBCType::Periodic, LinOpBCType::Periodic},
      {LinOpBCType::Periodic, LinOpBCType::Periodic, LinOpBCType::Periodic});

  vector<MultiFab> sigmas(ghext->leveldata.size());
  for (int level = 0; level < int(ghext->leveldata.size()); ++level) {
    auto &sigma = sigmas.at(level);
    sigma.define(ghext->amrcore->boxArray(level),
                 ghext->amrcore->DistributionMap(level), 1, 0);
    sigma.setVal(1.0);
    mlnodelaplacian.setSigma(level, sigma);
  }

  mlnodelaplacian.setVerbose(10);

  // Create solver

  MLMG mlmg(mlnodelaplacian);

  const int gi_rhs = CCTK_GroupIndex("Maxwell::dive");
  assert(gi_rhs >= 0);
  const int gi_sol = CCTK_GroupIndex("Maxwell::phi1");
  assert(gi_sol >= 0);
  const int gi_res = CCTK_GroupIndex("Maxwell::dive1");
  assert(gi_res >= 0);
  const int tl = 0;
  const int vi = 0;

  Vector<MultiFab *> ress(ghext->leveldata.size());
  Vector<MultiFab *> sols(ghext->leveldata.size());
  Vector<const MultiFab *> rhss(ghext->leveldata.size());
  for (int level = 0; level < int(ghext->leveldata.size()); ++level) {
    const auto &restrict leveldata = ghext->leveldata.at(level);
    const auto &restrict groupdata_rhs = *leveldata.groupdata.at(gi_rhs);
    rhss.at(level) = groupdata_rhs.mfab.at(tl).get();
    const auto &restrict groupdata_sol = *leveldata.groupdata.at(gi_sol);
    sols.at(level) = groupdata_sol.mfab.at(tl).get();
    const auto &restrict groupdata_res = *leveldata.groupdata.at(gi_res);
    ress.at(level) = groupdata_res.mfab.at(tl).get();
  }

  mlmg.setVerbose(10);
  mlmg.setBottomVerbose(10);

  // Initial guess

  for (int level = 0; level < int(ghext->leveldata.size()); ++level) {
    sols.at(level)->setVal(0.0, vi);
  }

  // Solve

  mlmg.compResidual(ress, sols, rhss);
  for (int level = 0; level < int(ghext->leveldata.size()); ++level)
    CCTK_VINFO("norm_inf rhs[%d]: %g", level,
               double(rhss.at(level)->norminf(0, 0, false, true)));
  for (int level = 0; level < int(ghext->leveldata.size()); ++level)
    CCTK_VINFO("norm_inf sol[%d]: %g", level,
               double(sols.at(level)->norminf(0, 0, false, true)));
  for (int level = 0; level < int(ghext->leveldata.size()); ++level)
    CCTK_VINFO("norm_inf res[%d]: %g", level,
               double(ress.at(level)->norminf(0, 0, false, true)));

  const CCTK_REAL rtol = 0.0;
  const CCTK_REAL atol = 1.0e-12;
  mlmg.solve(sols, rhss, rtol, atol);
}

} // namespace CarpetX
