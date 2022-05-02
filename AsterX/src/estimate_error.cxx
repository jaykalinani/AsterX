#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cmath>

namespace AsterX {
using namespace std;
using namespace Loop;

////////////////////////////////////////////////////////////////////////////////

extern "C" void AsterX_EstimateError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_EstimateError;
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
            err = max({err, fabs(var0 - varm), fabs(varp - var0)});
          }
          return err;
        };

        regrid_error(p.I) = max({calcerr(rho), calcerr(velx), calcerr(vely),
                                 calcerr(velz), calcerr(eps), calcerr(press)});
      });
}

} // namespace AsterX
