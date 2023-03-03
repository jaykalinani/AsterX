#include "defs.hxx"

#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>
#include <cmath>

namespace Hydro {
using namespace std;

////////////////////////////////////////////////////////////////////////////////

extern "C" void Hydro_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Hydro_Initialize;
  DECLARE_CCTK_PARAMETERS;

  const Loop::vect<CCTK_REAL, dim> x0{
      {CCTK_ORIGIN_SPACE(0), CCTK_ORIGIN_SPACE(1), CCTK_ORIGIN_SPACE(2)}};
  const Loop::vect<CCTK_REAL, dim> dx{
      {CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1), CCTK_DELTA_SPACE(2)}};

  const Loop::vect<int, dim> lbnd{{cctk_lbnd[0], cctk_lbnd[1], cctk_lbnd[2]}};
  const Loop::vect<int, dim> lsh{{cctk_lsh[0], cctk_lsh[1], cctk_lsh[2]}};
  const Loop::vect<int, dim> nghostzones{
      {cctk_nghostzones[0], cctk_nghostzones[1], cctk_nghostzones[2]}};

  const Loop::vect<int, dim> tmin{cctk_tile_min[0], cctk_tile_min[1], cctk_tile_min[2]};
  const Loop::vect<int, dim> tmax{cctk_tile_max[0], cctk_tile_max[1], cctk_tile_max[2]};

  const Loop::vect<int, dim> imin = max(tmin, nghostzones);
  const Loop::vect<int, dim> imax = min(tmax, lsh - nghostzones);

  const Loop::vect<int, dim> ash{{cctk_ash[0], cctk_ash[1], cctk_ash[2]}};
  constexpr ptrdiff_t di = 1;
  const ptrdiff_t dj = di * cctk_ash[0];
  const ptrdiff_t dk = dj * cctk_ash[1];

  if (CCTK_EQUALS(setup, "equilibrium")) {

    for (int k = imin[2]; k < imax[2]; ++k) {
      for (int j = imin[1]; j < imax[1]; ++j) {
        for (int i = imin[0]; i < imax[0]; i += vsize) {
          CCTK_BOOLVEC mask = mask_for_loop_tail<CCTK_BOOLVEC>(i, imax[0]);
          ptrdiff_t ind = i + dj * j + dk * k;
          CCTK_REALVEC rho1 = CCTK_REAL(1.0);
          CCTK_REALVEC velx1 = CCTK_REAL(0.0);
          CCTK_REALVEC vely1 = CCTK_REAL(0.0);
          CCTK_REALVEC velz1 = CCTK_REAL(0.0);
          CCTK_REALVEC eint1 = CCTK_REAL(1.0);
          mask_storeu(mask, &rho[ind], rho1);
          mask_storeu(mask, &velx[ind], velx1);
          mask_storeu(mask, &vely[ind], vely1);
          mask_storeu(mask, &velz[ind], velz1);
          mask_storeu(mask, &eint[ind], eint1);
        }
      }
    }

  } else if (CCTK_EQUALS(setup, "sound wave")) {

    for (int k = imin[2]; k < imax[2]; ++k) {
      for (int j = imin[1]; j < imax[1]; ++j) {
        for (int i = imin[0]; i < imax[0]; i += vsize) {
          CCTK_BOOLVEC mask = mask_for_loop_tail<CCTK_BOOLVEC>(i, imax[0]);
          ptrdiff_t ind = i + dj * j + dk * k;
          CCTK_REALVEC x1 =
              x0[0] + (lbnd[0] + i) * dx[0] + iota<CCTK_REALVEC>() * dx[0];
          CCTK_REALVEC rho1 = CCTK_REAL(1.0);
          CCTK_REALVEC velx1 =
              CCTK_REAL(0.0) + CCTK_REAL(amplitude) * sin(CCTK_REAL(M_PI) * x1);
          CCTK_REALVEC vely1 = CCTK_REAL(0.0);
          CCTK_REALVEC velz1 = CCTK_REAL(0.0);
          CCTK_REALVEC eint1 = CCTK_REAL(1.0);
          mask_storeu(mask, &rho[ind], rho1);
          mask_storeu(mask, &velx[ind], velx1);
          mask_storeu(mask, &vely[ind], vely1);
          mask_storeu(mask, &velz[ind], velz1);
          mask_storeu(mask, &eint[ind], eint1);
        }
      }
    }

  } else if (CCTK_EQUALS(setup, "shock tube")) {

    for (int k = imin[2]; k < imax[2]; ++k) {
      for (int j = imin[1]; j < imax[1]; ++j) {
        for (int i = imin[0]; i < imax[0]; i += vsize) {
          CCTK_BOOLVEC mask = mask_for_loop_tail<CCTK_BOOLVEC>(i, imax[0]);
          ptrdiff_t ind = i + dj * j + dk * k;
          CCTK_REALVEC x1 =
              x0[0] + (lbnd[0] + i) * dx[0] + iota<CCTK_REALVEC>() * dx[0];
          // left state
          CCTK_REALVEC rho1l = CCTK_REAL(2.0);
          CCTK_REALVEC velx1l = CCTK_REAL(0.0);
          CCTK_REALVEC vely1l = CCTK_REAL(0.0);
          CCTK_REALVEC velz1l = CCTK_REAL(0.0);
          CCTK_REALVEC eint1l = CCTK_REAL(2.0);
          // right state
          CCTK_REALVEC rho1r = CCTK_REAL(1.0);
          CCTK_REALVEC velx1r = CCTK_REAL(0.0);
          CCTK_REALVEC vely1r = CCTK_REAL(0.0);
          CCTK_REALVEC velz1r = CCTK_REAL(0.0);
          CCTK_REALVEC eint1r = CCTK_REAL(1.0);
          // choose
          CCTK_REALVEC rho1 = if_else(x1 <= CCTK_REAL(0.0), rho1l, rho1r);
          CCTK_REALVEC velx1 = if_else(x1 <= CCTK_REAL(0.0), velx1l, velx1r);
          CCTK_REALVEC vely1 = if_else(x1 <= CCTK_REAL(0.0), vely1l, vely1r);
          CCTK_REALVEC velz1 = if_else(x1 <= CCTK_REAL(0.0), velz1l, velz1r);
          CCTK_REALVEC eint1 = if_else(x1 <= CCTK_REAL(0.0), eint1l, eint1r);
          mask_storeu(mask, &rho[ind], rho1);
          mask_storeu(mask, &velx[ind], velx1);
          mask_storeu(mask, &vely[ind], vely1);
          mask_storeu(mask, &velz[ind], velz1);
          mask_storeu(mask, &eint[ind], eint1);
        }
      }
    }

  } else if (CCTK_EQUALS(setup, "spherical shock")) {

    for (int k = imin[2]; k < imax[2]; ++k) {
      for (int j = imin[1]; j < imax[1]; ++j) {
        for (int i = imin[0]; i < imax[0]; i += vsize) {
          CCTK_BOOLVEC mask = mask_for_loop_tail<CCTK_BOOLVEC>(i, imax[0]);
          ptrdiff_t ind = i + dj * j + dk * k;
          CCTK_REALVEC x1 =
              x0[0] + (lbnd[0] + i) * dx[0] + iota<CCTK_REALVEC>() * dx[0];
          CCTK_REALVEC y1 = x0[1] + (lbnd[1] + j) * dx[1];
          CCTK_REALVEC z1 = x0[2] + (lbnd[2] + k) * dx[2];
          CCTK_REALVEC r2 = pow2(x1) + (pow2(y1) + pow2(z1));
          CCTK_REAL sr2 = pow2(shock_radius);
          // inner state
          CCTK_REALVEC rho1i = CCTK_REAL(2.0);
          CCTK_REALVEC velx1i = CCTK_REAL(0.0);
          CCTK_REALVEC vely1i = CCTK_REAL(0.0);
          CCTK_REALVEC velz1i = CCTK_REAL(0.0);
          CCTK_REALVEC eint1i = CCTK_REAL(2.0);
          // outer state
          CCTK_REALVEC rho1o = CCTK_REAL(1.0);
          CCTK_REALVEC velx1o = CCTK_REAL(0.0);
          CCTK_REALVEC vely1o = CCTK_REAL(0.0);
          CCTK_REALVEC velz1o = CCTK_REAL(0.0);
          CCTK_REALVEC eint1o = CCTK_REAL(1.0);
          // choose
          CCTK_REALVEC rho1 = if_else(r2 <= sr2, rho1i, rho1o);
          CCTK_REALVEC velx1 = if_else(r2 <= sr2, velx1i, velx1o);
          CCTK_REALVEC vely1 = if_else(r2 <= sr2, vely1i, vely1o);
          CCTK_REALVEC velz1 = if_else(r2 <= sr2, velz1i, velz1o);
          CCTK_REALVEC eint1 = if_else(r2 <= sr2, eint1i, eint1o);
          mask_storeu(mask, &rho[ind], rho1);
          mask_storeu(mask, &velx[ind], velx1);
          mask_storeu(mask, &vely[ind], vely1);
          mask_storeu(mask, &velz[ind], velz1);
          mask_storeu(mask, &eint[ind], eint1);
        }
      }
    }

  } else {
    assert(0);
  }

  // Set an ignorable pressure
  for (int k = imin[2]; k < imax[2]; ++k) {
    for (int j = imin[1]; j < imax[1]; ++j) {
      for (int i = imin[0]; i < imax[0]; i += vsize) {
        CCTK_BOOLVEC mask = mask_for_loop_tail<CCTK_BOOLVEC>(i, imax[0]);
        ptrdiff_t ind = i + dj * j + dk * k;
        CCTK_REALVEC press1 = CCTK_REAL(-1.0);
        mask_storeu(mask, &press[ind], press1);
      }
    }
  }
}

} // namespace Hydro
