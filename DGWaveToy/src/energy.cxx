#include <dg.hxx>
#include <loop.hxx>
#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>

namespace DGWaveToy {
using namespace DG;
using namespace Arith;
using namespace std;

using Loop::dim;

extern "C" void DGWaveToy_Energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_DGWaveToy_Energy;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> ft_rhs_(cctkGH, ft_rhs);

  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> fx_(cctkGH, fx);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> fy_(cctkGH, fy);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> fz_(cctkGH, fz);

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> eps_(cctkGH, eps);

  const Loop::GridDescBase grid(cctkGH);

  switch (dg_npoints) {

  case 2: {
    constexpr int N = 2;
    grid.loop_int<1, 1, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
      const vect<CCTK_REAL, dim> dfx = get_derivs<N>(grid, fx_, p);
      const vect<CCTK_REAL, dim> dfy = get_derivs<N>(grid, fy_, p);
      const vect<CCTK_REAL, dim> dfz = get_derivs<N>(grid, fz_, p);

      eps_(p.I) = (pow(ft_rhs_(p.I), 2) + pow(dfx[0], 2) + pow(dfy[1], 2) +
                   pow(dfz[2], 2));
    });
    break;
  }

  case 4: {
    constexpr int N = 4;
    grid.loop_int<1, 1, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
      const vect<CCTK_REAL, dim> dfx = get_derivs<N>(grid, fx_, p);
      const vect<CCTK_REAL, dim> dfy = get_derivs<N>(grid, fy_, p);
      const vect<CCTK_REAL, dim> dfz = get_derivs<N>(grid, fz_, p);

      eps_(p.I) = (pow(ft_rhs_(p.I), 2) + pow(dfx[0], 2) + pow(dfy[1], 2) +
                   pow(dfz[2], 2));
    });
    break;
  }

  case 8: {
    constexpr int N = 8;
    grid.loop_int<1, 1, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
      const vect<CCTK_REAL, dim> dfx = get_derivs<N>(grid, fx_, p);
      const vect<CCTK_REAL, dim> dfy = get_derivs<N>(grid, fy_, p);
      const vect<CCTK_REAL, dim> dfz = get_derivs<N>(grid, fz_, p);

      eps_(p.I) = (pow(ft_rhs_(p.I), 2) + pow(dfx[0], 2) + pow(dfy[1], 2) +
                   pow(dfz[2], 2));
    });
    break;
  }

  default:
    assert(0);
  }
}

extern "C" void DGWaveToy_EnergyBoundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_DGWaveToy_EnergyBoundaries;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> eps_(cctkGH, eps);

  const Loop::GridDescBase grid(cctkGH);

  grid.loop_bnd<1, 1, 1>(grid.nghostzones,
                         [&](const Loop::PointDesc &p) { eps_(p.I) = 0; });
}

} // namespace DGWaveToy
