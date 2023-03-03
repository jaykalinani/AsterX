#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <nsimd/nsimd-all.hpp>

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

  typedef CCTK_INT8 CCTK_INTEGER;
  static_assert(sizeof(CCTK_INTEGER) == sizeof(CCTK_REAL), "");
  typedef nsimd::pack<CCTK_REAL> CCTK_VEC_REAL;
  typedef nsimd::pack<CCTK_INTEGER> CCTK_VEC_INTEGER;
  typedef nsimd::packl<CCTK_REAL> CCTK_VEC_BOOLEAN;
  constexpr int VS = sizeof(CCTK_VEC_REAL) / sizeof(CCTK_REAL);

  // Equation of state: p = (gamma - 1) e

  // vel^j = delta^j_i mom_i / rho

  for (int k = imin[2]; k < imax[2]; ++k) {
    for (int j = imin[1]; j < imax[1]; ++j) {
      ptrdiff_t idx_i0 = off + str[1] * j + str[2] * k;
      for (int i = imin[0]; i < imax[0]; i += VS) {
        ptrdiff_t idx = idx_i0 + i;

        auto rho_ = nsimd::loadu<CCTK_VEC_REAL>(&rho[idx]);
        auto momx_ = nsimd::loadu<CCTK_VEC_REAL>(&momx[idx]);
        auto momy_ = nsimd::loadu<CCTK_VEC_REAL>(&momy[idx]);
        auto momz_ = nsimd::loadu<CCTK_VEC_REAL>(&momz[idx]);
        auto etot_ = nsimd::loadu<CCTK_VEC_REAL>(&etot[idx]);

        auto rho_inv = nsimd::set1<CCTK_VEC_REAL>(1.0) / rho_;

        auto ekin = nsimd::set1<CCTK_VEC_REAL>(0.5) * rho_inv *
                    sqrt(pow2(momx_) + pow2(momy_) + pow2(momz_));
        auto eint_ = etot_ - ekin;

        auto press_ = nsimd::set1<CCTK_VEC_REAL>(gamma - 1) * eint_;

        auto velx_ = rho_inv * momx_;
        auto vely_ = rho_inv * momy_;
        auto velz_ = rho_inv * momz_;

        if (i + VS >= imax[0]) {
          auto iota = nsimd::iota<CCTK_VEC_INTEGER>();
          auto imask = iota < nsimd::set1<CCTK_VEC_INTEGER>(
                                  CCTK_INTEGER(imax[0] - (i + VS)));
          auto mask = reinterpretl(CCTK_VEC_BOOLEAN(), imask);
          nsimd::storeu_masked(&press[idx], press_, mask);
          nsimd::storeu_masked(&velx[idx], velx_, mask);
          nsimd::storeu_masked(&vely[idx], vely_, mask);
          nsimd::storeu_masked(&velz[idx], velz_, mask);
          nsimd::storeu_masked(&eint[idx], eint_, mask);
          break;
        }
        nsimd::storeu(&press[idx], press_);
        nsimd::storeu(&velx[idx], velx_);
        nsimd::storeu(&vely[idx], vely_);
        nsimd::storeu(&velz[idx], velz_);
        nsimd::storeu(&eint[idx], eint_);
      }
    }
  }
}

} // namespace HydroToyCarpetX
