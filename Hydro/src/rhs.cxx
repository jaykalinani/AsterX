#include "defs.hxx"

#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

namespace Hydro {
using namespace std;

////////////////////////////////////////////////////////////////////////////////

extern "C" void Hydro_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Hydro_RHS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);
  const CCTK_REAL dV1 = 1 / (dx * dy * dz);

  const Loop::vect<int, dim> nx{{1, 0, 0}};
  const Loop::vect<int, dim> ny{{0, 1, 0}};
  const Loop::vect<int, dim> nz{{0, 0, 1}};

  const Loop::vect<int, dim> lsh{{cctk_lsh[0], cctk_lsh[1], cctk_lsh[2]}};
  const Loop::vect<int, dim> nghostzones{
      {cctk_nghostzones[0], cctk_nghostzones[1], cctk_nghostzones[2]}};

  array<CCTK_INT, dim> tmin_, tmax_;
  GetTileExtent(cctkGH, tmin_.data(), tmax_.data());
  const Loop::vect<int, dim> tmin(tmin_);
  const Loop::vect<int, dim> tmax(tmax_);

  const Loop::vect<int, dim> imin = max(tmin, nghostzones);
  const Loop::vect<int, dim> imax = min(tmax, lsh - nghostzones);

  const Loop::vect<int, dim> ash{{cctk_ash[0], cctk_ash[1], cctk_ash[2]}};
  constexpr ptrdiff_t di = 1;
  const ptrdiff_t dj = di * ash[0];
  const ptrdiff_t dk = dj * ash[1];

  // x
  // fluxes: face-centred without ghosts
  const Loop::vect<int, dim> ashx = ash + nx - 2 * nghostzones;
  constexpr ptrdiff_t dix = 1;
  const ptrdiff_t djx = dix * ashx[0];
  const ptrdiff_t dkx = djx * ashx[1];
  const ptrdiff_t offx =
      nghostzones[0] * dix + nghostzones[1] * djx + nghostzones[2] * dkx;

  // y
  // fluxes: face-centred without ghosts
  const Loop::vect<int, dim> ashy = ash + ny - 2 * nghostzones;
  constexpr ptrdiff_t diy = 1;
  const ptrdiff_t djy = diy * ashy[0];
  const ptrdiff_t dky = djy * ashy[1];
  const ptrdiff_t offy =
      nghostzones[0] * diy + nghostzones[1] * djy + nghostzones[2] * dky;

  // z
  // fluxes: face-centred without ghosts
  const Loop::vect<int, dim> ashz = ash + nz - 2 * nghostzones;
  constexpr ptrdiff_t diz = 1;
  const ptrdiff_t djz = diz * ashz[0];
  const ptrdiff_t dkz = djz * ashz[1];
  const ptrdiff_t offz =
      nghostzones[0] * diz + nghostzones[1] * djz + nghostzones[2] * dkz;

  const auto calcrhs =
      [&](const CCTK_REAL *restrict const fxvar,
          const CCTK_REAL *restrict const fyvar,
          const CCTK_REAL *restrict const fzvar,
          CCTK_REAL *restrict const dtvar) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        for (int k = imin[2]; k < imax[2]; ++k) {
          for (int j = imin[1]; j < imax[1]; ++j) {
            for (int i = imin[0]; i < imax[0]; i += vsize) {
              ptrdiff_t ind = i + dj * j + dk * k;
              ptrdiff_t indx = i + djx * j + dkx * k - offx;
              ptrdiff_t indy = i + djy * j + dky * k - offy;
              ptrdiff_t indz = i + djz * j + dkz * k - offz;

              CCTK_REALVEC dtvar1 =
                  dV1 * ((vloadu(fxvar[indx + dix]) - vloadu(fxvar[indx])) +
                         (vloadu(fyvar[indy + djy]) - vloadu(fyvar[indy])) +
                         (vloadu(fzvar[indz + dkz]) - vloadu(fzvar[indz])));

              dtvar1.storeu_partial(dtvar[ind], i, imin[0], imax[0]);
            }
          }
        }
      };

  calcrhs(fxdens, fydens, fzdens, dtdens);
  calcrhs(fxmomx, fymomx, fzmomx, dtmomx);
  calcrhs(fxmomy, fymomy, fzmomy, dtmomy);
  calcrhs(fxmomz, fymomz, fzmomz, dtmomz);
  calcrhs(fxetot, fyetot, fzetot, dtetot);
}

} // namespace Hydro
