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

namespace {

// LLF (local Lax-Friedrichs) Riemann solver
template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE T llf(const array<T, 2> &var,
                                          const array<T, 2> &flux) {
  return T(0.5) * ((flux[0] + flux[1]) - (var[1] - var[0]));
}

} // namespace

////////////////////////////////////////////////////////////////////////////////

extern "C" void Hydro_Fluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Hydro_Fluxes;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);
  const CCTK_REAL dAx = dy * dz;
  const CCTK_REAL dAy = dx * dz;
  const CCTK_REAL dAz = dx * dy;

  const Loop::vect<int, dim> zero{{0, 0, 0}};
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

  const Loop::vect<int, dim> ash{{cctk_ash[0], cctk_ash[1], cctk_ash[2]}};
  constexpr ptrdiff_t di = 1;
  const ptrdiff_t dj = di * ash[0];
  const ptrdiff_t dk = dj * ash[1];

  const auto calcflux = [&](const int dir, const int ddir,
                            const CCTK_REAL dAdir,
                            const Loop::vect<int, dim> &ndir,
                            const CCTK_REAL *restrict const veldir,
                            CCTK_REAL *restrict const fdirdens,
                            CCTK_REAL *restrict const fdirmomx,
                            CCTK_REAL *restrict const fdirmomy,
                            CCTK_REAL *restrict const fdirmomz,
                            CCTK_REAL *restrict const
                                fdiretot) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const Loop::vect<int, dim> imindir = max(tmin, nghostzones);
    const Loop::vect<int, dim> imaxdir =
        min(tmax + (tmax >= lsh).ifelse(ndir, zero), lsh + ndir - nghostzones);

    // fluxes: face-centred without ghosts
    const Loop::vect<int, dim> ashdir = ash + ndir - 2 * nghostzones;
    constexpr ptrdiff_t didir = 1;
    const ptrdiff_t djdir = didir * ashdir[0];
    const ptrdiff_t dkdir = djdir * ashdir[1];
    const ptrdiff_t offdir = nghostzones[0] * didir + nghostzones[1] * djdir +
                             nghostzones[2] * dkdir;

    for (int k = imindir[2]; k < imaxdir[2]; ++k) {
      for (int j = imindir[1]; j < imaxdir[1]; ++j) {
        for (int i = imindir[0]; i < imaxdir[0]; i += vsize) {
          ptrdiff_t ind = i + dj * j + dk * k;
          ptrdiff_t inddir = i + djdir * j + dkdir * k - offdir;

          // Read conserved and primitive variables
          array<CCTK_REALVEC, 2> dens2, momx2, momy2, momz2, etot2, press2,
              veldir2;
          for (int n = 0; n < 2; ++n) {
            dens2[n] = vloadu(dens[ind + (n - 1) * ddir]);
            momx2[n] = vloadu(momx[ind + (n - 1) * ddir]);
            momy2[n] = vloadu(momy[ind + (n - 1) * ddir]);
            momz2[n] = vloadu(momz[ind + (n - 1) * ddir]);
            etot2[n] = vloadu(etot[ind + (n - 1) * ddir]);
            press2[n] = vloadu(press[ind + (n - 1) * ddir]);
            veldir2[n] = vloadu(veldir[ind + (n - 1) * ddir]);
          }

#warning "TODO: This calculates the fluxes twice; change this"
          // Calculate centred fluxes
          array<CCTK_REALVEC, 2> fdirdens2, fdirmomx2, fdirmomy2, fdirmomz2,
              fdiretot2;
          for (int n = 0; n < 2; ++n) {
            fdirdens2[n] = dens2[n] * veldir2[n];
            fdirmomx2[n] = momx2[n] * veldir2[n] + (dir == 0 ? press2[n] : 0);
            fdirmomy2[n] = momy2[n] * veldir2[n] + (dir == 1 ? press2[n] : 0);
            fdirmomz2[n] = momz2[n] * veldir2[n] + (dir == 2 ? press2[n] : 0);
            fdiretot2[n] = (etot2[n] + press2[n]) * veldir2[n];
          }

          // Calculate face fluxes
          CCTK_REALVEC fdirdens1 = dAdir * llf(dens2, fdirdens2);
          CCTK_REALVEC fdirmomx1 = dAdir * llf(momx2, fdirmomx2);
          CCTK_REALVEC fdirmomy1 = dAdir * llf(momy2, fdirmomy2);
          CCTK_REALVEC fdirmomz1 = dAdir * llf(momz2, fdirmomz2);
          CCTK_REALVEC fdiretot1 = dAdir * llf(etot2, fdiretot2);

          // Store fluxes
          fdirdens1.storeu_partial(fdirdens[inddir], i, imindir[0], imaxdir[0]);
          fdirmomx1.storeu_partial(fdirmomx[inddir], i, imindir[0], imaxdir[0]);
          fdirmomy1.storeu_partial(fdirmomy[inddir], i, imindir[0], imaxdir[0]);
          fdirmomz1.storeu_partial(fdirmomz[inddir], i, imindir[0], imaxdir[0]);
          fdiretot1.storeu_partial(fdiretot[inddir], i, imindir[0], imaxdir[0]);
        }
      }
    }
  };

  calcflux(0, di, dAx, nx, velx, fxdens, fxmomx, fxmomy, fxmomz, fxetot);
  calcflux(1, dj, dAy, ny, vely, fydens, fymomx, fymomy, fymomz, fyetot);
  calcflux(2, dk, dAz, nz, velz, fzdens, fzmomx, fzmomy, fzmomz, fzetot);
}

} // namespace Hydro
