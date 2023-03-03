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
using namespace std;

using Loop::dim;

extern "C" void DGWaveToy_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_DGWaveToy_RHS;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> u_(cctkGH, u);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> ft_(cctkGH, ft);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> fx_(cctkGH, fx);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> fy_(cctkGH, fy);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> fz_(cctkGH, fz);

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> u_rhs_(cctkGH, u_rhs);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> ft_rhs_(cctkGH, ft_rhs);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> fx_rhs_(cctkGH, fx_rhs);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> fy_rhs_(cctkGH, fy_rhs);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> fz_rhs_(cctkGH, fz_rhs);

  const Loop::GridDescBase grid(cctkGH);

  grid.loop_all<1, 1, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
    u_rhs_(p.I) = ft_(p.I);
  });

  switch (dg_npoints) {

  case 2: {
    constexpr int N = 2;
    grid.loop_int<1, 1, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
      const vect<CCTK_REAL, dim> dft = get_derivs<N>(grid, ft_, p);
      const vect<CCTK_REAL, dim> dfx = get_derivs<N>(grid, fx_, p);
      const vect<CCTK_REAL, dim> dfy = get_derivs<N>(grid, fy_, p);
      const vect<CCTK_REAL, dim> dfz = get_derivs<N>(grid, fz_, p);

      ft_rhs_(p.I) = dfx[0] + dfy[1] + dfz[2];
      fx_rhs_(p.I) = dft[0];
      fy_rhs_(p.I) = dft[1];
      fz_rhs_(p.I) = dft[2];
    });
    break;
  }

  case 4: {
    constexpr int N = 4;
    grid.loop_int<1, 1, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
      const vect<CCTK_REAL, dim> dft = get_derivs<N>(grid, ft_, p);
      const vect<CCTK_REAL, dim> dfx = get_derivs<N>(grid, fx_, p);
      const vect<CCTK_REAL, dim> dfy = get_derivs<N>(grid, fy_, p);
      const vect<CCTK_REAL, dim> dfz = get_derivs<N>(grid, fz_, p);

      ft_rhs_(p.I) = dfx[0] + dfy[1] + dfz[2];
      fx_rhs_(p.I) = dft[0];
      fy_rhs_(p.I) = dft[1];
      fz_rhs_(p.I) = dft[2];
    });
    break;
  }

  case 8: {
    constexpr int N = 8;
    grid.loop_int<1, 1, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
      const vect<CCTK_REAL, dim> dft = get_derivs<N>(grid, ft_, p);
      const vect<CCTK_REAL, dim> dfx = get_derivs<N>(grid, fx_, p);
      const vect<CCTK_REAL, dim> dfy = get_derivs<N>(grid, fy_, p);
      const vect<CCTK_REAL, dim> dfz = get_derivs<N>(grid, fz_, p);

      ft_rhs_(p.I) = dfx[0] + dfy[1] + dfz[2];
      fx_rhs_(p.I) = dft[0];
      fy_rhs_(p.I) = dft[1];
      fz_rhs_(p.I) = dft[2];
    });
    break;
  }

  default:
    assert(0);
  }
}

extern "C" void DGWaveToy_RHSBoundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_DGWaveToy_RHSBoundaries;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> ft_rhs_(cctkGH, ft_rhs);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> fx_rhs_(cctkGH, fx_rhs);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> fy_rhs_(cctkGH, fy_rhs);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> fz_rhs_(cctkGH, fz_rhs);

  const Loop::GridDescBase grid(cctkGH);

  grid.loop_bnd<1, 1, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
    ft_rhs_(p.I) = 0;
    fx_rhs_(p.I) = 0;
    fy_rhs_(p.I) = 0;
    fz_rhs_(p.I) = 0;
  });
}

} // namespace DGWaveToy
