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

namespace {
template <typename T> T pow2(T x) { return x * x; }
} // namespace

extern "C" void HydroToyCarpetX_Pressure(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyCarpetX_Pressure;
  DECLARE_CCTK_PARAMETERS;

  const Loop::vect<int, dim> zero{{0, 0, 0}};

  const Loop::vect<int, dim> lsh{{cctk_lsh[0], cctk_lsh[1], cctk_lsh[2]}};

  array<CCTK_INT, dim> tmin_, tmax_;
  GetTileExtent(cctkGH, tmin_.data(), tmax_.data());
  const Loop::vect<int, dim> tmin(tmin_);
  const Loop::vect<int, dim> tmax(tmax_);

  const Loop::vect<int, dim> imin = max(tmin, zero);
  const Loop::vect<int, dim> imax = min(tmax, lsh);

  // regular, cell centred variables
  const Loop::vect<int, dim> ash{{cctk_ash[0], cctk_ash[1], cctk_ash[2]}};

  const auto mkstr{[](const Loop::vect<int, dim> &ash) {
    Loop::vect<ptrdiff_t, dim> str;
    ptrdiff_t s = 1;
    for (int d = 0; d < dim; ++d) {
      str[d] = s;
      s *= ash[d];
    }
    return str;
  }};

  const Loop::vect<ptrdiff_t, dim> str = mkstr(ash);

  constexpr ptrdiff_t off = 0;

  typedef vectype<CCTK_REAL> CCTK_VEC_REAL;
  constexpr int VS = vecprops<CCTK_REAL>::size();
  // assert((imax[0] - imin[0]) % VS == 0);

  // Equation of state: p = (gamma - 1) e

  // vel^j = delta^j_i mom_i / rho

  for (int k = imin[2]; k < imax[2]; ++k) {
    for (int j = imin[1]; j < imax[1]; ++j) {
      ptrdiff_t idx_i0 = off + str[1] * j + str[2] * k;
      for (int i = imin[0]; i < imax[0]; i += VS) {
        ptrdiff_t idx = idx_i0 + i;

        auto rho_ = CCTK_VEC_REAL::loadu(rho[idx]);
        auto momx_ = CCTK_VEC_REAL::loadu(momx[idx]);
        auto momy_ = CCTK_VEC_REAL::loadu(momy[idx]);
        auto momz_ = CCTK_VEC_REAL::loadu(momz[idx]);
        auto etot_ = CCTK_VEC_REAL::loadu(etot[idx]);

        auto rho_inv = CCTK_VEC_REAL::set1(1.0) / rho_;

        auto ekin = CCTK_VEC_REAL::set1(0.5) * rho_inv *
                    sqrt(pow2(momx_) + pow2(momy_) + pow2(momz_));
        auto eint_ = etot_ - ekin;

        auto press_ = CCTK_VEC_REAL::set1(gamma - 1) * eint_;

        auto velx_ = rho_inv * momx_;
        auto vely_ = rho_inv * momy_;
        auto velz_ = rho_inv * momz_;

        press_.storeu_partial(press[idx], i, i, imax[0]);
        velx_.storeu_partial(velx[idx], i, i, imax[0]);
        vely_.storeu_partial(vely[idx], i, i, imax[0]);
        velz_.storeu_partial(velz[idx], i, i, imax[0]);
        eint_.storeu_partial(eint[idx], i, i, imax[0]);
      }
    }
  }
}

} // namespace HydroToyCarpetX
