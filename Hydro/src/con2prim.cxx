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

  array<CCTK_INT, dim> tmin_, tmax_;
  GetTileExtent(cctkGH, tmin_.data(), tmax_.data());
  const Loop::vect<int, dim> tmin(tmin_);
  const Loop::vect<int, dim> tmax(tmax_);

  const Loop::vect<int, dim> imin = tmin;
  const Loop::vect<int, dim> imax = tmax;

  const Loop::vect<int, dim> ash{{cctk_ash[0], cctk_ash[1], cctk_ash[2]}};
  constexpr ptrdiff_t di = 1;
  const ptrdiff_t dj = di * cctk_ash[0];
  const ptrdiff_t dk = dj * cctk_ash[1];

  for (int k = imin[2]; k < imax[2]; ++k) {
    for (int j = imin[1]; j < imax[1]; ++j) {
      for (int i = imin[0]; i < imax[0]; i += VS) {
        ptrdiff_t ind = i + dj * j + dk * k;

        CCTK_REALVEC rho1 = vloadu(rho[ind]);
        CCTK_REALVEC momx1 = vloadu(momx[ind]);
        CCTK_REALVEC momy1 = vloadu(momy[ind]);
        CCTK_REALVEC momz1 = vloadu(momz[ind]);
        CCTK_REALVEC etot1 = vloadu(etot[ind]);

        CCTK_REALVEC rho_inv = CCTK_REAL(1.0) / rho1;

        CCTK_REALVEC ekin = CCTK_REAL(0.5) * rho_inv *
                            sqrt(pow2(momx1) + pow2(momy1) + pow2(momz1));
        CCTK_REALVEC eint1 = etot1 - ekin;

        // Equation of state: p = (gamma - 1) e
        CCTK_REALVEC press1 = CCTK_REAL(gamma - 1) * eint1;

        // vel^j = delta^j_i mom_i / rho
        CCTK_REALVEC velx1 = rho_inv * momx1;
        CCTK_REALVEC vely1 = rho_inv * momy1;
        CCTK_REALVEC velz1 = rho_inv * momz1;

        press1.storeu_partial(press[ind], i, imin[0], imax[0]);
        velx1.storeu_partial(velx[ind], i, imin[0], imax[0]);
        vely1.storeu_partial(vely[ind], i, imin[0], imax[0]);
        velz1.storeu_partial(velz[ind], i, imin[0], imax[0]);
        eint1.storeu_partial(eint[ind], i, imin[0], imax[0]);
      }
    }
  }
}

} // namespace Hydro
