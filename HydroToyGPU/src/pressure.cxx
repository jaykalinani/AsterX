#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace HydroToyGPU {
using namespace std;
using namespace Loop;

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T pow2(T x) {
  return x * x;
}

extern "C" void HydroToyGPU_Pressure(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyGPU_Pressure;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  constexpr array<int, dim> cell_centred = {1, 1, 1};
  const GF3D2layout gf_layout(cctkGH, cell_centred);

  const GF3D2<const CCTK_REAL> gf_rho(gf_layout, rho);
  const GF3D2<const CCTK_REAL> gf_momx(gf_layout, momx);
  const GF3D2<const CCTK_REAL> gf_momy(gf_layout, momy);
  const GF3D2<const CCTK_REAL> gf_momz(gf_layout, momz);
  const GF3D2<const CCTK_REAL> gf_etot(gf_layout, etot);

  const GF3D2<CCTK_REAL> gf_press(gf_layout, press);
  const GF3D2<CCTK_REAL> gf_velx(gf_layout, velx);
  const GF3D2<CCTK_REAL> gf_vely(gf_layout, vely);
  const GF3D2<CCTK_REAL> gf_velz(gf_layout, velz);
  const GF3D2<CCTK_REAL> gf_eint(gf_layout, eint);

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE CCTK_HOST(
                            const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL rho_inv = 1.0 / (gf_rho(p.I) + 1.0e-20);

        CCTK_REAL ekin =
            0.5 * rho_inv *
            sqrt(pow2(gf_momx(p.I)) + pow2(gf_momy(p.I)) + pow2(gf_momz(p.I)));
        gf_eint(p.I) = gf_etot(p.I) - ekin;

        gf_press(p.I) = (gamma - 1) * gf_eint(p.I);

        gf_velx(p.I) = rho_inv * gf_momx(p.I);
        gf_vely(p.I) = rho_inv * gf_momy(p.I);
        gf_velz(p.I) = rho_inv * gf_momz(p.I);
      });
}

} // namespace HydroToyGPU
