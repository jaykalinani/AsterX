#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
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
  DECLARE_CCTK_ARGUMENTS_GRHydroToyGPU_EstimateError;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  constexpr array<int, dim> cell_centred = {1, 1, 1};
  const GF3D2layout gf_layout(cctkGH, cell_centred);

  const GF3D2<const CCTK_REAL> gf_rho(gf_layout, rho);
  const GF3D2<const CCTK_REAL> gf_velx(gf_layout, velx);
  const GF3D2<const CCTK_REAL> gf_vely(gf_layout, vely);
  const GF3D2<const CCTK_REAL> gf_velz(gf_layout, velz);
  const GF3D2<const CCTK_REAL> gf_eps(gf_layout, eps);
  const GF3D2<const CCTK_REAL> gf_press(gf_layout, press);

  const GF3D2<CCTK_REAL> gf_regrid_error(gf_layout, regrid_error);

  constexpr auto DI = PointDesc::DI;
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE CCTK_HOST(
                            const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
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

        gf_regrid_error(p.I) =
            fmax6(calcerr(gf_rho), calcerr(gf_velx), calcerr(gf_vely),
                  calcerr(gf_velz), calcerr(gf_eps), calcerr(gf_press));
      });
}

} // namespace GRHydroToyGPU
