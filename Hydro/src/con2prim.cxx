#include "defs.hxx"

#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <array>
#include <cmath>

namespace Hydro {
using namespace std;

////////////////////////////////////////////////////////////////////////////////

extern "C" void Hydro_Con2prim(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Hydro_Con2prim;
  DECLARE_CCTK_PARAMETERS;

  const Loop::vect<int, dim> tmin{cctk_tile_min[0], cctk_tile_min[1], cctk_tile_min[2]};
  const Loop::vect<int, dim> tmax{cctk_tile_max[0], cctk_tile_max[1], cctk_tile_max[2]};

  const Loop::vect<int, dim> imin = tmin;
  const Loop::vect<int, dim> imax = tmax;

  const Loop::vect<int, dim> ash{{cctk_ash[0], cctk_ash[1], cctk_ash[2]}};
  constexpr ptrdiff_t di = 1;
  const ptrdiff_t dj = di * cctk_ash[0];
  const ptrdiff_t dk = dj * cctk_ash[1];

  for (int k = imin[2]; k < imax[2]; ++k) {
    for (int j = imin[1]; j < imax[1]; ++j) {
      for (int i = imin[0]; i < imax[0]; i += vsize) {
        CCTK_BOOLVEC mask = mask_for_loop_tail<CCTK_BOOLVEC>(i, imax[0]);
        ptrdiff_t ind = i + dj * j + dk * k;

        CCTK_REALVEC dens1 = masko_loadu(mask, &dens[ind], 1);
        CCTK_REALVEC momx1 = maskz_loadu(mask, &momx[ind]);
        CCTK_REALVEC momy1 = maskz_loadu(mask, &momy[ind]);
        CCTK_REALVEC momz1 = maskz_loadu(mask, &momz[ind]);
        CCTK_REALVEC etot1 = masko_loadu(mask, &etot[ind], 1);

        CCTK_REALVEC rho1 = dens1;

        CCTK_REALVEC rho1_inv = 1 / rho1;

        CCTK_REALVEC ekin1 = CCTK_REAL(0.5) * rho1_inv *
                             sqrt(pow2(momx1) + pow2(momy1) + pow2(momz1));
        CCTK_REALVEC eint1 = etot1 - ekin1;

        // Equation of state: p = (gamma - 1) e
        CCTK_REALVEC press1 = (gamma - 1) * eint1;

        // vel^j = delta^j_i mom_i / rho
        CCTK_REALVEC velx1 = rho1_inv * momx1;
        CCTK_REALVEC vely1 = rho1_inv * momy1;
        CCTK_REALVEC velz1 = rho1_inv * momz1;

        mask_storeu(mask, &rho[ind], rho1);
        mask_storeu(mask, &velx[ind], velx1);
        mask_storeu(mask, &vely[ind], vely1);
        mask_storeu(mask, &velz[ind], velz1);
        mask_storeu(mask, &eint[ind], eint1);
        mask_storeu(mask, &press[ind], press1);
      }
    }
  }
}

} // namespace Hydro
