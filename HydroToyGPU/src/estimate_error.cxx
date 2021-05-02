#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace HydroToyGPU {
using namespace std;
using namespace Loop;

////////////////////////////////////////////////////////////////////////////////

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE T fmax5(T x0, T x1,
                                                                  T x2, T x3,
                                                                  T x4) {
  T r01 = fmax(x0, x1);
  T r23 = fmax(x2, x3);
  T r014 = fmax(r01, x4);
  T r01234 = fmax(r014, r23);
  return r01234;
}

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
      [=](const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE {
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
                fmax5(calcerr(gf_rho), calcerr(gf_momx), calcerr(gf_momy),
                      calcerr(gf_momz), calcerr(gf_etot));
          });
}

} // namespace HydroToyGPU
