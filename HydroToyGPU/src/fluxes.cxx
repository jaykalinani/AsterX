#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cmath>

namespace HydroToyGPU {
using namespace std;
using namespace Loop;

// Calculate the fluxes in direction `dir`. This function is more
// complex because it has to handle any direction, but as reward,
// there is only one function, not three.
template <int dir> void CalcFlux(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyGPU_Fluxes;
  DECLARE_CCTK_PARAMETERS;

  static_assert(dir >= 0 && dir < 3, "");

  // const array<CCTK_REAL, dim> dx = {CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1),
  //                                   CCTK_DELTA_SPACE(2)};
  // const CCTK_REAL dV = dx[0] * dx[1] * dx[2]; // cell volume
  // const CCTK_REAL dA = dV / dx[dir];          // face area

  // Cell centred grid functions
  const GridDescBaseDevice grid(cctkGH);
  constexpr array<int, dim> cell_centred = {1, 1, 1};
  const GF3D2layout gf_layout(cctkGH, cell_centred);

  const GF3D2<const CCTK_REAL> gf_rho(gf_layout, rho);
  const GF3D2<const CCTK_REAL> gf_momx(gf_layout, momx);
  const GF3D2<const CCTK_REAL> gf_momy(gf_layout, momy);
  const GF3D2<const CCTK_REAL> gf_momz(gf_layout, momz);
  const GF3D2<const CCTK_REAL> gf_etot(gf_layout, etot);

  // Get the grid function pointer for velocity in direction `dir`
  const array<const CCTK_REAL *, dim> vels = {velx, vely, velz};
  const GF3D2<const CCTK_REAL> gf_vel(gf_layout, vels[dir]);
  const GF3D2<const CCTK_REAL> gf_press(gf_layout, press);

  // Face-centred grid functions (in direction `dir`)
  constexpr array<int, dim> face_centred = {!(dir == 0), !(dir == 1),
                                            !(dir == 2)};
  const GF3D2layout gf_fluxlayout(cctkGH, face_centred);

  // Get the grid function pointers for fluxes in direction `dir`
  const array<CCTK_REAL *, dim> fluxrhos = {fxrho, fyrho, fzrho};
  const array<CCTK_REAL *, dim> fluxmomxs = {fxmomx, fymomx, fzmomx};
  const array<CCTK_REAL *, dim> fluxmomys = {fxmomy, fymomy, fzmomy};
  const array<CCTK_REAL *, dim> fluxmomzs = {fxmomz, fymomz, fzmomz};
  const array<CCTK_REAL *, dim> fluxetots = {fxetot, fyetot, fzetot};
  const GF3D2<CCTK_REAL> gf_fluxrho(gf_fluxlayout, fluxrhos[dir]);
  const GF3D2<CCTK_REAL> gf_fluxmomx(gf_fluxlayout, fluxmomxs[dir]);
  const GF3D2<CCTK_REAL> gf_fluxmomy(gf_fluxlayout, fluxmomys[dir]);
  const GF3D2<CCTK_REAL> gf_fluxmomz(gf_fluxlayout, fluxmomzs[dir]);
  const GF3D2<CCTK_REAL> gf_fluxetot(gf_fluxlayout, fluxetots[dir]);

  // frho^i = rho vel^i
  // fmom^i_j = mom_j vel^i + delta^i_j press
  // fetot^i = (etot + press) vel^i

  const auto calcflux =
      [=](CCTK_REAL var_m, CCTK_REAL var_p, CCTK_REAL flux_m, CCTK_REAL flux_p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE {
            CCTK_REAL lambda_m = 1.0;
            CCTK_REAL lambda_p = -1.0;
            CCTK_REAL llf =
                0.5 * ((flux_m + flux_p) -
                       fmax(fabs(lambda_m), fabs(lambda_p)) * (var_p - var_m));
            // return dA * llf;
            return llf;
          };

  constexpr auto DI = PointDesc::DI;
  grid.loop_int_device<face_centred[0], face_centred[1], face_centred[2]>(
      grid.nghostzones,
      [=](const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE {
            // Neighbouring "plus" and "minus" cell indices
            const auto Im = p.I - DI[dir];
            const auto Ip = p.I;
            gf_fluxrho(p.I) =
                calcflux(gf_rho(Im), gf_rho(Ip), gf_rho(Im) * gf_vel(Im),
                         gf_rho(Ip) * gf_vel(Ip));
            gf_fluxmomx(p.I) =
                calcflux(gf_momx(Im), gf_momx(Ip),
                         gf_momx(Im) * gf_vel(Im) + (dir == 0) * gf_press(Im),
                         gf_momx(Ip) * gf_vel(Ip) + (dir == 0) * gf_press(Ip));
            gf_fluxmomy(p.I) =
                calcflux(gf_momy(Im), gf_momy(Ip),
                         gf_momy(Im) * gf_vel(Im) + (dir == 1) * gf_press(Im),
                         gf_momy(Ip) * gf_vel(Ip) + (dir == 1) * gf_press(Ip));
            gf_fluxmomz(p.I) =
                calcflux(gf_momz(Im), gf_momz(Ip),
                         gf_momz(Im) * gf_vel(Im) + (dir == 2) * gf_press(Im),
                         gf_momz(Ip) * gf_vel(Ip) + (dir == 2) * gf_press(Ip));
            gf_fluxetot(p.I) =
                calcflux(gf_etot(Im), gf_etot(Ip),
                         (gf_etot(Im) + gf_press(Im)) * gf_vel(Im),
                         (gf_etot(Ip) + gf_press(Ip)) * gf_vel(Ip));
          });
}

extern "C" void HydroToyGPU_Fluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyGPU_Fluxes;
  DECLARE_CCTK_PARAMETERS;

  CalcFlux<0>(cctkGH);
  CalcFlux<1>(cctkGH);
  CalcFlux<2>(cctkGH);
}

} // namespace HydroToyGPU
