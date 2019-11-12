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

extern "C" void Hydro_Prim2Con(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Hydro_Prim2Con;
  DECLARE_CCTK_PARAMETERS;

  const Loop::vect<int, dim> lsh{{cctk_lsh[0], cctk_lsh[1], cctk_lsh[2]}};

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
      for (int i = imin[0]; i < imax[0]; i += vsize) {
        ptrdiff_t ind = i + dj * j + dk * k;

        CCTK_REALVEC rho1 = vloadu(rho[ind]);
        CCTK_REALVEC velx1 = vloadu(velx[ind]);
        CCTK_REALVEC vely1 = vloadu(vely[ind]);
        CCTK_REALVEC velz1 = vloadu(velz[ind]);
        CCTK_REALVEC eint1 = vloadu(eint[ind]);

        CCTK_REALVEC ekin1 =
            CCTK_REAL(0.5) * rho1 * (pow2(velx1) + pow2(vely1) + pow2(velz1));

        CCTK_REALVEC dens1 = rho1;
        CCTK_REALVEC momx1 = rho1 * velx1;
        CCTK_REALVEC momy1 = rho1 * vely1;
        CCTK_REALVEC momz1 = rho1 * velz1;
        CCTK_REALVEC etot1 = ekin1 + eint1;

        dens1.storeu_partial(dens[ind], i, imin[0], imax[0]);
        momx1.storeu_partial(momx[ind], i, imin[0], imax[0]);
        momy1.storeu_partial(momy[ind], i, imin[0], imax[0]);
        momz1.storeu_partial(momz[ind], i, imin[0], imax[0]);
        etot1.storeu_partial(etot[ind], i, imin[0], imax[0]);
      }
    }
  }
}

} // namespace Hydro
