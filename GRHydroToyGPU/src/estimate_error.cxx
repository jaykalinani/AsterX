#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace GRHydroToyGPU {
using namespace std;
using namespace Loop;

////////////////////////////////////////////////////////////////////////////////

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T fmax6(T x0, T x1,
                                                                  T x2, T x3,
                                                                  T x4, T x5) {
  T r01 = fmax(x0, x1);
  T r23 = fmax(x2, x3);
  T r014 = fmax(r01, x4);
  T r01234 = fmax(r014, r23);
  T r012345 = fmax(r01234, x5);
  return r012345;
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void GRHydroToyGPU_EstimateError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHydroToyGPU_EstimateError;
  DECLARE_CCTK_PARAMETERS;

  constexpr auto DI = PointDesc::DI;
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const auto calcerr = [&](auto &gf_var) {
          CCTK_REAL err{0};
          for (int d = 0; d < dim; ++d) {
            auto varm = gf_var(p.I - DI[d]);
            auto var0 = gf_var(p.I);
            auto varp = gf_var(p.I + DI[d]);
            // Calculate derivative
            err = fmax(err, fmax(fabs(var0 - varm), fabs(varp - var0)));
          }
          return err;
        };

        regrid_error(p.I) = fmax6(calcerr(rho), calcerr(velx), calcerr(vely),
                                  calcerr(velz), calcerr(eps), calcerr(press));
      });
}

} // namespace GRHydroToyGPU
