#include <loop.hxx>

#include <vectors.h>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>

namespace HydroToyCarpetX {
using namespace std;

constexpr int dim = 3;

extern "C" void HydroToyCarpetX_Evolve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyCarpetX_Evolve;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);
  const CCTK_REAL dV1 = 1 / (dx * dy * dz);

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

  const Loop::vect<int, dim> imin = max(tmin, nghostzones);
  const Loop::vect<int, dim> imax = min(tmax, lsh - nghostzones);

  // regular, cell centred variables
  const Loop::vect<int, dim> ash{{cctk_ash[0], cctk_ash[1], cctk_ash[2]}};
  // variables without ghosts
  const Loop::vect<int, dim> ash_ng = ash - 2 * nghostzones;

  // fluxes :face-centred without ghosts
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
  // assert((imax[0] - imin[0]) % VS == 0);

  const auto vec_dV1 = CCTK_VEC_REAL::set1(dV1);

  // Transport
  // dt rho + d_i (rho vel^i) = 0
  // dt mom_j + d_i (mom_j vel^i) = 0
  // dt etot + d_i (etot vel^i) = 0

  const auto calcupdate = [&](auto &u_, const auto &u_p_, const auto &fx_,
                              const auto &fy_, const auto &fz_) {
    for (int k = imin[2]; k < imax[2]; ++k) {
      for (int j = imin[1]; j < imax[1]; ++j) {
        ptrdiff_t idx_i0 = off + str[1] * j + str[2] * k;
        ptrdiff_t idx_fx_i0 = off_fx + str_fx[1] * j + str_fx[2] * k;
        ptrdiff_t idx_fy_i0 = off_fy + str_fy[1] * j + str_fy[2] * k;
        ptrdiff_t idx_fz_i0 = off_fz + str_fz[1] * j + str_fz[2] * k;
        for (int i = imin[0]; i < imax[0]; i += VS) {
          ptrdiff_t idx = idx_i0 + i;
          ptrdiff_t idx_fx = idx_fx_i0 + i;
          ptrdiff_t idx_fy = idx_fy_i0 + i;
          ptrdiff_t idx_fz = idx_fz_i0 + i;
          auto u_p = CCTK_VEC_REAL::loadu(u_p_[idx]);
          auto fx1 = CCTK_VEC_REAL::loadu(fx_[idx_fx + str_fx[0]]);
          auto fx0 = CCTK_VEC_REAL::loadu(fx_[idx_fx]);
          auto fy1 = CCTK_VEC_REAL::loadu(fy_[idx_fy + str_fy[1]]);
          auto fy0 = CCTK_VEC_REAL::loadu(fy_[idx_fy]);
          auto fz1 = CCTK_VEC_REAL::loadu(fz_[idx_fz + str_fz[2]]);
          auto fz0 = CCTK_VEC_REAL::loadu(fz_[idx_fz]);
          auto u = u_p - vec_dV1 * ((fx1 - fx0) + (fy1 - fy0) + (fz1 - fz0));
          u.storeu_partial(u_[idx], i, i, imax[0]);
        }
      }
    }
  };

  calcupdate(rho, rho_p, fxrho, fyrho, fzrho);
  calcupdate(momx, momx_p, fxmomx, fymomx, fzmomx);
  calcupdate(momy, momy_p, fxmomy, fymomy, fzmomy);
  calcupdate(momz, momz_p, fxmomz, fymomz, fzmomz);
  calcupdate(etot, etot_p, fxetot, fyetot, fzetot);
}

} // namespace HydroToyCarpetX
