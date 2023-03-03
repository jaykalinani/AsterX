#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cmath>

namespace HydroToyGPU {
using namespace std;
using namespace Loop;

////////////////////////////////////////////////////////////////////////////////

extern "C" void HydroToyGPU_EstimateError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyGPU_EstimateError;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  constexpr array<int, dim> cell_centred = {1, 1, 1};
  const GF3D2layout gf_layout(cctkGH, cell_centred);

  const GF3D2<const CCTK_REAL> gf_rho(gf_layout, rho);
  const GF3D2<const CCTK_REAL> gf_momx(gf_layout, momx);
  const GF3D2<const CCTK_REAL> gf_momy(gf_layout, momy);
  const GF3D2<const CCTK_REAL> gf_momz(gf_layout, momz);
  const GF3D2<const CCTK_REAL> gf_etot(gf_layout, etot);

  const GF3D2<CCTK_REAL> gf_regrid_error(gf_layout, regrid_error);

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

        gf_regrid_error(p.I) =
            max({calcerr(gf_rho), calcerr(gf_momx), calcerr(gf_momy),
                 calcerr(gf_momz), calcerr(gf_etot)});
      });
}

} // namespace HydroToyGPU
