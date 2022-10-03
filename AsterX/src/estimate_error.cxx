#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cmath>

namespace AsterX {
using namespace std;
using namespace Loop;

enum class regrid_t {
  first_deriv_all,
  first_deriv_rho,
  second_deriv_rho,
  first_grad_rho
};

/* calculate max of d(var)/(dx) * hx in all dirs */
template <typename T>
CCTK_DEVICE CCTK_HOST const T calc_grad_1st(const GF3D2<const T> &gf,
                                            const PointDesc &p) {
  constexpr auto DI = PointDesc::DI;
  CCTK_REAL err{0}, errp{0}, errm{0};
  for (int d = 0; d < dim; ++d) {
    auto varm = gf(p.I - DI[d]);
    auto var0 = gf(p.I);
    auto varp = gf(p.I + DI[d]);
    errp += (varp - var0) * (varp - var0);
    errm += (var0 - varm) * (var0 - varm);
  }
  err = max({errp, errm});
  return sqrt(err);
}

/* calculate max of d(var)/(dx) * hx in all dirs */
template <typename T>
CCTK_DEVICE CCTK_HOST const T calc_deriv_1st(const GF3D2<const T> &gf,
                                             const PointDesc &p) {
  constexpr auto DI = PointDesc::DI;
  CCTK_REAL err{0};
  for (int d = 0; d < dim; ++d) {
    auto varm = gf(p.I - DI[d]);
    auto var0 = gf(p.I);
    auto varp = gf(p.I + DI[d]);
    err = max({err, fabs(var0 - varm), fabs(varp - var0)});
  }
  return err;
}

/* calculate max of d^2(var)/(dx^2) * hx^2 in all dirs */
template <typename T>
CCTK_DEVICE CCTK_HOST const T calc_deriv_2nd(const GF3D2<const T> &gf,
                                             const PointDesc &p) {
  constexpr auto DI = PointDesc::DI;
  CCTK_REAL err{0};
  for (int d = 0; d < dim; ++d) {
    auto varm = gf(p.I - DI[d]);
    auto var0 = gf(p.I);
    auto varp = gf(p.I + DI[d]);
    err = max({err, fabs(varp + varm - 2.0 * var0)});
  }
  return err;
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void AsterX_EstimateError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_EstimateError;
  DECLARE_CCTK_PARAMETERS;

  regrid_t regrid;
  if (CCTK_EQUALS(regrid_method, "1st deriv"))
    regrid = regrid_t::first_deriv_all;
  else if (CCTK_EQUALS(regrid_method, "1st deriv of rho"))
    regrid = regrid_t::first_deriv_rho;
  else if (CCTK_EQUALS(regrid_method, "2nd deriv of rho"))
    regrid = regrid_t::second_deriv_rho;
  else if (CCTK_EQUALS(regrid_method, "1st grad of rho"))
    regrid = regrid_t::first_grad_rho;
  else
    CCTK_ERROR("Unknown value for parameter \"regrid_method\"");

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        switch (regrid) {

        case regrid_t::first_deriv_all: {
          regrid_error(p.I) =
              max({calc_deriv_1st(rho, p), calc_deriv_1st(velx, p),
                   calc_deriv_1st(vely, p), calc_deriv_1st(velz, p),
                   calc_deriv_1st(eps, p), calc_deriv_1st(press, p)});
          break;
        }

        case regrid_t::first_deriv_rho: {
          regrid_error(p.I) = calc_deriv_1st(rho, p);
          break;
        }

        case regrid_t::second_deriv_rho: {
          regrid_error(p.I) = calc_deriv_2nd(rho, p);
          break;
        }

        case regrid_t::first_grad_rho: {
          regrid_error(p.I) = calc_grad_1st(rho, p);
          break;
        }

        default:
          assert(0);
        }
      });
}

} // namespace AsterX
