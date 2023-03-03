#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

// #ifdef AMREX_USE_GPU
// #include <AMReX_GpuDevice.H>
// #endif

#include <array>

namespace GRHydroToyGPU {
using namespace Loop;

extern "C" void GRHydroToyGPU_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_GRHydroToyGPU_RHS;
  DECLARE_CCTK_PARAMETERS;

  const std::array dx{CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1),
                      CCTK_DELTA_SPACE(2)};
  const std::array dx1{1 / dx[0], 1 / dx[1], 1 / dx[2]};
  // const CCTK_REAL dV = dx[0] * dx[1] * dx[2]; // cell volume
  // const CCTK_REAL dV1 = 1 / dV;

  // Transport
  // dt dens + d_i (dens vel^i) = 0
  // dt mom_j + d_i (mom_j vel^i) = 0
  // dt tau + d_i (tau vel^i) = 0

  const auto calcupdate =
      [=] CCTK_DEVICE(CCTK_REAL fx_m, CCTK_REAL fx_p, CCTK_REAL fy_m,
                      CCTK_REAL fy_p, CCTK_REAL fz_m, CCTK_REAL fz_p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
            // return -dV1 * (fx_p - fx_m + fy_p - fy_m + fz_p - fz_m);
            return -dx1[0] * (fx_p - fx_m) - dx1[1] * (fy_p - fy_m) -
                   dx1[2] * (fz_p - fz_m);
          };

  constexpr auto DI = PointDesc::DI;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Neighbouring "plus" and "minus" face indices in
        // the x, y, and z directions
        const auto Imx = p.I;
        const auto Imy = p.I;
        const auto Imz = p.I;
        const auto Ipx = p.I + DI[0];
        const auto Ipy = p.I + DI[1];
        const auto Ipz = p.I + DI[2];

        densrhs(p.I) = calcupdate(fxdens(Imx), fxdens(Ipx), fydens(Imy),
                                  fydens(Ipy), fzdens(Imz), fzdens(Ipz));
      });

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Neighbouring "plus" and "minus" face indices in
        // the x, y, and z directions
        const auto Imx = p.I;
        const auto Imy = p.I;
        const auto Imz = p.I;
        const auto Ipx = p.I + DI[0];
        const auto Ipy = p.I + DI[1];
        const auto Ipz = p.I + DI[2];

        momxrhs(p.I) = calcupdate(fxmomx(Imx), fxmomx(Ipx), fymomx(Imy),
                                  fymomx(Ipy), fzmomx(Imz), fzmomx(Ipz));
      });

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Neighbouring "plus" and "minus" face indices in
        // the x, y, and z directions
        const auto Imx = p.I;
        const auto Imy = p.I;
        const auto Imz = p.I;
        const auto Ipx = p.I + DI[0];
        const auto Ipy = p.I + DI[1];
        const auto Ipz = p.I + DI[2];

        momyrhs(p.I) = calcupdate(fxmomy(Imx), fxmomy(Ipx), fymomy(Imy),
                                  fymomy(Ipy), fzmomy(Imz), fzmomy(Ipz));
      });

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Neighbouring "plus" and "minus" face indices in
        // the x, y, and z directions
        const auto Imx = p.I;
        const auto Imy = p.I;
        const auto Imz = p.I;
        const auto Ipx = p.I + DI[0];
        const auto Ipy = p.I + DI[1];
        const auto Ipz = p.I + DI[2];

        momzrhs(p.I) = calcupdate(fxmomz(Imx), fxmomz(Ipx), fymomz(Imy),
                                  fymomz(Ipy), fzmomz(Imz), fzmomz(Ipz));
      });

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Neighbouring "plus" and "minus" face indices in
        // the x, y, and z directions
        const auto Imx = p.I;
        const auto Imy = p.I;
        const auto Imz = p.I;
        const auto Ipx = p.I + DI[0];
        const auto Ipy = p.I + DI[1];
        const auto Ipz = p.I + DI[2];

        taurhs(p.I) = calcupdate(fxtau(Imx), fxtau(Ipx), fytau(Imy), fytau(Ipy),
                                 fztau(Imz), fztau(Ipz));
      });
}

} // namespace GRHydroToyGPU
