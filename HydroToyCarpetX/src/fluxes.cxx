#include <loop.hxx>

#include <vectors.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>

namespace HydroToyCarpetX {
using namespace std;

constexpr int dim = 3;

extern "C" void HydroToyCarpetX_Fluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyCarpetX_Fluxes;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);
  const CCTK_REAL dAx = dt * dy * dz;
  const CCTK_REAL dAy = dt * dx * dz;
  const CCTK_REAL dAz = dt * dx * dy;

  const Loop::vect<int, dim> zero{{0, 0, 0}};
  const Loop::vect<int, dim> di{{1, 0, 0}};
  const Loop::vect<int, dim> dj{{0, 1, 0}};
  const Loop::vect<int, dim> dk{{0, 0, 1}};

  const Loop::vect<int, dim> lsh{{cctk_lsh[0], cctk_lsh[1], cctk_lsh[2]}};
  const Loop::vect<int, dim> nghostzones{
      {cctk_nghostzones[0], cctk_nghostzones[1], cctk_nghostzones[2]}};

  array<CCTK_INT, dim> tmin_, tmax_;
  GetTileExtent(cctkGH, tmin_.data(), tmax_.data());
  const Loop::vect<int, dim> tmin(tmin_);
  const Loop::vect<int, dim> tmax(tmax_);

  // const Loop::vect<int, dim> imin = max(tmin, nghostzones);
  // const Loop::vect<int, dim> imax = min(tmax, lsh - nghostzones);

  const Loop::vect<int, dim> imin_fx = max(tmin, nghostzones);
  const Loop::vect<int, dim> imax_fx =
      min(tmax + (tmax >= lsh).ifelse(di, zero), lsh + di - nghostzones);
  const Loop::vect<int, dim> imin_fy = max(tmin, nghostzones);
  const Loop::vect<int, dim> imax_fy =
      min(tmax + (tmax >= lsh).ifelse(di, zero), lsh + dj - nghostzones);
  const Loop::vect<int, dim> imin_fz = max(tmin, nghostzones);
  const Loop::vect<int, dim> imax_fz =
      min(tmax + (tmax >= lsh).ifelse(di, zero), lsh + dk - nghostzones);

  // regular, cell centred variables
  const Loop::vect<int, dim> ash{{cctk_ash[0], cctk_ash[1], cctk_ash[2]}};
  // variables without ghosts
  const Loop::vect<int, dim> ash_ng = ash - 2 * nghostzones;

  // fluxes: face-centred without ghosts
  const Loop::vect<int, dim> ash_fx = ash_ng + di;
  const Loop::vect<int, dim> ash_fy = ash_ng + dj;
  const Loop::vect<int, dim> ash_fz = ash_ng + dk;

  const auto mkstr = [](const Loop::vect<int, dim> &ash) {
    Loop::vect<ptrdiff_t, dim> str;
    ptrdiff_t s = 1;
    for (int d = 0; d < dim; ++d) {
      str[d] = s;
      s *= ash[d];
    }
    return str;
  };

  const Loop::vect<ptrdiff_t, dim> str = mkstr(ash);
  const Loop::vect<ptrdiff_t, dim> str_fx = mkstr(ash_fx);
  const Loop::vect<ptrdiff_t, dim> str_fy = mkstr(ash_fy);
  const Loop::vect<ptrdiff_t, dim> str_fz = mkstr(ash_fz);

  constexpr ptrdiff_t off = 0;
  const ptrdiff_t off_fx = -str_fx[0] - str_fx[1] - str_fx[2];
  const ptrdiff_t off_fy = -str_fy[0] - str_fy[1] - str_fy[2];
  const ptrdiff_t off_fz = -str_fz[0] - str_fz[1] - str_fz[2];

  typedef vectype<CCTK_REAL> CCTK_VEC_REAL;
  constexpr int VS = vecprops<CCTK_REAL>::size();
  // assert((imax_fx[0] - imin_fx[0]) % VS == 0);
  // assert((imax_fy[0] - imin_fy[0]) % VS == 0);
  // assert((imax_fz[0] - imin_fz[0]) % VS == 0);

  const auto vec_dAx = CCTK_VEC_REAL::set1(dAx);
  const auto vec_dAy = CCTK_VEC_REAL::set1(dAy);
  const auto vec_dAz = CCTK_VEC_REAL::set1(dAz);

  // LLF (local Lax-Friedrichs) Riemann solver
  const auto llf = [](auto lambda_m, auto lambda_p, auto var_m, auto var_p,
                      auto flux_m, auto flux_p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    return CCTK_VEC_REAL::set1(0.5) *
           ((flux_m + flux_p) -
            fmax(fabs(lambda_m), fabs(lambda_p)) * (var_p - var_m));
  };

  // Fluxes
  // frho^i = rho vel^i
  // fmom^i_j = mom_j vel^i + delta^i_j press
  // fetot^i = (etot + press) vel^i

  // For refluxing, AMReX expects fluxes to be integrated over the
  // faces both in space and time.

  const auto calcflux_x = [&](auto &fx_, const auto &u_, const auto &flux_x) {
    for (int k = imin_fx[2]; k < imax_fx[2]; ++k) {
      for (int j = imin_fx[1]; j < imax_fx[1]; ++j) {
        ptrdiff_t idx_i0 = off + str[1] * j + str[2] * k;
        ptrdiff_t idx_fx_i0 = off_fx + str_fx[1] * j + str_fx[2] * k;
        for (int i = imin_fx[0]; i < imax_fx[0]; i += VS) {
          ptrdiff_t idx = idx_i0 + i;
          ptrdiff_t idx_fx = idx_fx_i0 + i;

          CCTK_VEC_REAL lambda_m(1.0);
          CCTK_VEC_REAL lambda_p(-1.0);
          auto var_m = CCTK_VEC_REAL::loadu(u_[idx - str[0]]);
          auto var_p = CCTK_VEC_REAL::loadu(u_[idx]);
          auto flux_m = flux_x(idx - str[0]);
          auto flux_p = flux_x(idx);
          auto fx =
              vec_dAx * llf(lambda_m, lambda_p, var_m, var_p, flux_m, flux_p);
          fx.storeu_partial(fx_[idx_fx], i, i, imax_fx[0]);
        }
      }
    }
  };

  calcflux_x(fxrho, rho, [&](ptrdiff_t idx) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    return CCTK_VEC_REAL::loadu(rho[idx]) * CCTK_VEC_REAL::loadu(velx[idx]);
  });
  calcflux_x(fxmomx, momx, [&](ptrdiff_t idx) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    return CCTK_VEC_REAL::loadu(momx[idx]) * CCTK_VEC_REAL::loadu(velx[idx]) +
           CCTK_VEC_REAL::loadu(press[idx]);
  });
  calcflux_x(fxmomy, momy, [&](ptrdiff_t idx) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    return CCTK_VEC_REAL::loadu(momy[idx]) * CCTK_VEC_REAL::loadu(velx[idx]);
  });
  calcflux_x(fxmomz, momz, [&](ptrdiff_t idx) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    return CCTK_VEC_REAL::loadu(momz[idx]) * CCTK_VEC_REAL::loadu(velx[idx]);
  });
  calcflux_x(fxetot, etot, [&](ptrdiff_t idx) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    return (CCTK_VEC_REAL::loadu(etot[idx]) +
            CCTK_VEC_REAL::loadu(press[idx])) *
           CCTK_VEC_REAL::loadu(velx[idx]);
  });

  const auto calcflux_y = [&](auto &fy_, const auto &u_, const auto &flux_y) {
    for (int k = imin_fy[2]; k < imax_fy[2]; ++k) {
      for (int j = imin_fy[1]; j < imax_fy[1]; ++j) {
        ptrdiff_t idx_i0 = off + str[1] * j + str[2] * k;
        ptrdiff_t idx_fy_i0 = off_fy + str_fy[1] * j + str_fy[2] * k;
        for (int i = imin_fy[0]; i < imax_fy[0]; i += VS) {
          ptrdiff_t idx = idx_i0 + i;
          ptrdiff_t idx_fy = idx_fy_i0 + i;

          CCTK_VEC_REAL lambda_m(1.0);
          CCTK_VEC_REAL lambda_p(-1.0);
          auto var_m = CCTK_VEC_REAL::loadu(u_[idx - str[1]]);
          auto var_p = CCTK_VEC_REAL::loadu(u_[idx]);
          auto flux_m = flux_y(idx - str[1]);
          auto flux_p = flux_y(idx);
          auto fy =
              vec_dAy * llf(lambda_m, lambda_p, var_m, var_p, flux_m, flux_p);
          fy.storeu_partial(fy_[idx_fy], i, i, imax_fy[0]);
        }
      }
    }
  };

  calcflux_y(fyrho, rho, [&](ptrdiff_t idx) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    return CCTK_VEC_REAL::loadu(rho[idx]) * CCTK_VEC_REAL::loadu(vely[idx]);
  });
  calcflux_y(fymomx, momx, [&](ptrdiff_t idx) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    return CCTK_VEC_REAL::loadu(momx[idx]) * CCTK_VEC_REAL::loadu(vely[idx]);
  });
  calcflux_y(fymomy, momy, [&](ptrdiff_t idx) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    return CCTK_VEC_REAL::loadu(momy[idx]) * CCTK_VEC_REAL::loadu(vely[idx]) +
           CCTK_VEC_REAL::loadu(press[idx]);
  });
  calcflux_y(fymomz, momz, [&](ptrdiff_t idx) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    return CCTK_VEC_REAL::loadu(momz[idx]) * CCTK_VEC_REAL::loadu(vely[idx]);
  });
  calcflux_y(fyetot, etot, [&](ptrdiff_t idx) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    return (CCTK_VEC_REAL::loadu(etot[idx]) +
            CCTK_VEC_REAL::loadu(press[idx])) *
           CCTK_VEC_REAL::loadu(vely[idx]);
  });

  const auto calcflux_z = [&](auto &fz_, const auto &u_, const auto &flux_z) {
    for (int k = imin_fz[2]; k < imax_fz[2]; ++k) {
      for (int j = imin_fz[1]; j < imax_fz[1]; ++j) {
        ptrdiff_t idx_i0 = off + str[1] * j + str[2] * k;
        ptrdiff_t idx_fz_i0 = off_fz + str_fz[1] * j + str_fz[2] * k;
        for (int i = imin_fz[0]; i < imax_fz[0]; i += VS) {
          ptrdiff_t idx = idx_i0 + i;
          ptrdiff_t idx_fz = idx_fz_i0 + i;

          CCTK_VEC_REAL lambda_m(1.0);
          CCTK_VEC_REAL lambda_p(-1.0);
          auto var_m = CCTK_VEC_REAL::loadu(u_[idx - str[2]]);
          auto var_p = CCTK_VEC_REAL::loadu(u_[idx]);
          auto flux_m = flux_z(idx - str[2]);
          auto flux_p = flux_z(idx);
          auto fz =
              vec_dAz * llf(lambda_m, lambda_p, var_m, var_p, flux_m, flux_p);
          fz.storeu_partial(fz_[idx_fz], i, i, imax_fz[0]);
        }
      }
    }
  };

  calcflux_z(fzrho, rho, [&](ptrdiff_t idx) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    return CCTK_VEC_REAL::loadu(rho[idx]) * CCTK_VEC_REAL::loadu(velz[idx]);
  });
  calcflux_z(fzmomx, momx, [&](ptrdiff_t idx) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    return CCTK_VEC_REAL::loadu(momx[idx]) * CCTK_VEC_REAL::loadu(velz[idx]);
  });
  calcflux_z(fzmomy, momy, [&](ptrdiff_t idx) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    return CCTK_VEC_REAL::loadu(momy[idx]) * CCTK_VEC_REAL::loadu(velz[idx]);
  });
  calcflux_z(fzmomz, momz, [&](ptrdiff_t idx) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    return CCTK_VEC_REAL::loadu(momz[idx]) * CCTK_VEC_REAL::loadu(velz[idx]) +
           CCTK_VEC_REAL::loadu(press[idx]);
  });
  calcflux_z(fzetot, etot, [&](ptrdiff_t idx) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    return (CCTK_VEC_REAL::loadu(etot[idx]) +
            CCTK_VEC_REAL::loadu(press[idx])) *
           CCTK_VEC_REAL::loadu(velz[idx]);
  });
}

} // namespace HydroToyCarpetX
