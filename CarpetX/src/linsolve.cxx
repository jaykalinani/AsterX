#include "driver.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <AMReX.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>

#include <vector>

namespace CarpetX {
using namespace std;

extern "C" void CarpetX_SolvePoisson() {
  // Create operator

  amrex::Vector<amrex::Geometry> geoms(ghext->leveldata.size());
  amrex::Vector<amrex::BoxArray> grids(ghext->leveldata.size());
  amrex::Vector<amrex::DistributionMapping> dmaps(ghext->leveldata.size());
  for (int level = 0; level < int(ghext->leveldata.size()); ++level) {
    geoms.at(level) = ghext->amrcore->Geom(level);
    grids.at(level) = ghext->amrcore->boxArray(level);
    dmaps.at(level) = ghext->amrcore->DistributionMap(level);
  }

  amrex::MLNodeLaplacian mlnodelaplacian(geoms, grids, dmaps);

  // TODO: Look at Geometry for this
  mlnodelaplacian.setDomainBC(
      {amrex::LinOpBCType::Dirichlet, amrex::LinOpBCType::Dirichlet,
       amrex::LinOpBCType::Dirichlet},
      {amrex::LinOpBCType::Dirichlet, amrex::LinOpBCType::Dirichlet,
       amrex::LinOpBCType::Dirichlet});

  vector<amrex::MultiFab> sigmas(ghext->leveldata.size());
  for (int level = 0; level < int(ghext->leveldata.size()); ++level) {
    auto &sigma = sigmas.at(level);
    sigma.define(ghext->amrcore->boxArray(level),
                 ghext->amrcore->DistributionMap(level), 1, 0);
    sigma.setVal(1.0);
    mlnodelaplacian.setSigma(level, sigma);
  }

  mlnodelaplacian.setVerbose(10);

  // Create solver

  amrex::MLMG mlmg(mlnodelaplacian);

  const int gi_rhs = CCTK_GroupIndex("Poisson::rhs");
  assert(gi_rhs >= 0);
  const int gi_sol = CCTK_GroupIndex("Poisson::phi");
  assert(gi_sol >= 0);
  const int gi_res = CCTK_GroupIndex("Poisson::res");
  assert(gi_res >= 0);
  const int tl = 0;
  const int vi = 0;

  amrex::Vector<amrex::MultiFab *> ress(ghext->leveldata.size());
  amrex::Vector<amrex::MultiFab *> sols(ghext->leveldata.size());
  amrex::Vector<const amrex::MultiFab *> rhss(ghext->leveldata.size());
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

  // Solve

  mlmg.compResidual(ress, sols, rhss);
#pragma omp critical
  {
    CCTK_VINFO("Before solving:");
    for (int i = -1; i < 10; ++i) {
      const int j = 4, k = 4;
      CCTK_VINFO("phi[%d,%d,%d]=%g", i, j, k,
                 double(sols.at(0)->array(0)(i, j, k, vi)));
    }
    for (int level = 0; level < int(ghext->leveldata.size()); ++level)
      CCTK_VINFO("norm_inf rhs[%d]: %g", level,
                 double(rhss.at(level)->norminf(vi, 0, false, true)));
    for (int level = 0; level < int(ghext->leveldata.size()); ++level)
      CCTK_VINFO("norm_inf sol[%d]: %g", level,
                 double(sols.at(level)->norminf(vi, 0, false, true)));
    for (int level = 0; level < int(ghext->leveldata.size()); ++level)
      CCTK_VINFO("norm_inf res[%d]: %g", level,
                 double(ress.at(level)->norminf(vi, 0, false, true)));
  }

  const CCTK_REAL rtol = 0.0;
  const CCTK_REAL atol = 1.0e-12;
  const CCTK_REAL maxerr = mlmg.solve(sols, rhss, rtol, atol);
#pragma omp critical
  CCTK_VINFO("Solution error (norm_inf): %g", double(maxerr));

  mlmg.compResidual(ress, sols, rhss);
#pragma omp critical
  {
    CCTK_VINFO("After solving:");
    for (int i = -1; i < 10; ++i) {
      const int j = 4, k = 4;
      CCTK_VINFO("phi[%d,%d,%d]=%g", i, j, k,
                 double(sols.at(0)->array(0)(i, j, k, vi)));
    }
    for (int level = 0; level < int(ghext->leveldata.size()); ++level)
      CCTK_VINFO("norm_inf rhs[%d]: %g", level,
                 double(rhss.at(level)->norminf(vi, 0, false, true)));
    for (int level = 0; level < int(ghext->leveldata.size()); ++level)
      CCTK_VINFO("norm_inf sol[%d]: %g", level,
                 double(sols.at(level)->norminf(vi, 0, false, true)));
    for (int level = 0; level < int(ghext->leveldata.size()); ++level)
      CCTK_VINFO("norm_inf res[%d]: %g", level,
                 double(ress.at(level)->norminf(vi, 0, false, true)));
  }
}

} // namespace CarpetX
