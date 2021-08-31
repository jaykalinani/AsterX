#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

// #ifdef AMREX_USE_GPU
// #include <AMReX_GpuDevice.H>
// #endif

#include <array>

namespace GRHydroToyGPU {
using namespace std;
using namespace Loop;

extern "C" void GRHydroToyGPU_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHydroToyGPU_RHS;
  DECLARE_CCTK_PARAMETERS;

  const array<CCTK_REAL, dim> dx = {CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1),
                                    CCTK_DELTA_SPACE(2)};
  const array<CCTK_REAL, dim> dx1 = {1 / dx[0], 1 / dx[1], 1 / dx[2]};
  // const CCTK_REAL dV = dx[0] * dx[1] * dx[2]; // cell volume
  // const CCTK_REAL dV1 = 1 / dV;

  // Cell centred grid functions
  const GridDescBaseDevice grid(cctkGH);
  constexpr array<int, dim> cell_centred = {1, 1, 1};
  const GF3D2layout gf_layout(cctkGH, cell_centred);

  const GF3D2<CCTK_REAL> gf_densrhs(gf_layout, densrhs);
  const GF3D2<CCTK_REAL> gf_momxrhs(gf_layout, momxrhs);
  const GF3D2<CCTK_REAL> gf_momyrhs(gf_layout, momyrhs);
  const GF3D2<CCTK_REAL> gf_momzrhs(gf_layout, momzrhs);
  const GF3D2<CCTK_REAL> gf_taurhs(gf_layout, taurhs);

  // Fluxes in x direction
  const GF3D2layout gf_fxlayout(cctkGH, {0, 1, 1});
  const GF3D2<const CCTK_REAL> gf_fxdens(gf_fxlayout, fxdens);
  const GF3D2<const CCTK_REAL> gf_fxmomx(gf_fxlayout, fxmomx);
  const GF3D2<const CCTK_REAL> gf_fxmomy(gf_fxlayout, fxmomy);
  const GF3D2<const CCTK_REAL> gf_fxmomz(gf_fxlayout, fxmomz);
  const GF3D2<const CCTK_REAL> gf_fxtau(gf_fxlayout, fxtau);

  // Fluxes in y direction
  const GF3D2layout gf_fylayout(cctkGH, {1, 0, 1});
  const GF3D2<const CCTK_REAL> gf_fydens(gf_fylayout, fydens);
  const GF3D2<const CCTK_REAL> gf_fymomx(gf_fylayout, fymomx);
  const GF3D2<const CCTK_REAL> gf_fymomy(gf_fylayout, fymomy);
  const GF3D2<const CCTK_REAL> gf_fymomz(gf_fylayout, fymomz);
  const GF3D2<const CCTK_REAL> gf_fytau(gf_fylayout, fytau);

  // Fluxes in z direction
  const GF3D2layout gf_fzlayout(cctkGH, {1, 1, 0});
  const GF3D2<const CCTK_REAL> gf_fzdens(gf_fzlayout, fzdens);
  const GF3D2<const CCTK_REAL> gf_fzmomx(gf_fzlayout, fzmomx);
  const GF3D2<const CCTK_REAL> gf_fzmomy(gf_fzlayout, fzmomy);
  const GF3D2<const CCTK_REAL> gf_fzmomz(gf_fzlayout, fzmomz);
  const GF3D2<const CCTK_REAL> gf_fztau(gf_fzlayout, fztau);

  // Transport
  // dt dens + d_i (dens vel^i) = 0
  // dt mom_j + d_i (mom_j vel^i) = 0
  // dt tau + d_i (tau vel^i) = 0

  const auto calcupdate =
      [=] CCTK_DEVICE CCTK_HOST(CCTK_REAL fx_m, CCTK_REAL fx_p, CCTK_REAL fy_m,
                                CCTK_REAL fy_p, CCTK_REAL fz_m, CCTK_REAL fz_p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
            // return -dV1 * (fx_p - fx_m + fy_p - fy_m + fz_p - fz_m);
            return -dx1[0] * (fx_p - fx_m) - dx1[1] * (fy_p - fy_m) -
                   dx1[2] * (fz_p - fz_m);
          };

  constexpr auto DI = PointDesc::DI;

#if 0
  // This kernel fails on GPUs for unknown reasons. Maybe it is too
  // complex? The alternative implementation below works fine.

#ifdef AMREX_USE_GPU
  AMREX_GPU_ERROR_CHECK();
  amrex::Gpu::synchronize();
  AMREX_GPU_ERROR_CHECK();
#endif

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE CCTK_HOST(
                            const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE  {
        // Neighbouring "plus" and "minus" face indices in the x, y, and z
        // directions
        const auto Imx = p.I;
        const auto Imy = p.I;
        const auto Imz = p.I;
        const auto Ipx = p.I + DI[0];
        const auto Ipy = p.I + DI[1];
        const auto Ipz = p.I + DI[2];

        gf_rhsdens(p.I) =
            calcupdate(gf_fxdens(Imx), gf_fxrho(Ipx), gf_fyrho(Imy),
                       gf_fydens(Ipy), gf_fzrho(Imz), gf_fzrho(Ipz));
        gf_rhsmomx(p.I) =
            calcupdate(gf_fxmomx(Imx), gf_fxmomx(Ipx), gf_fymomx(Imy),
                       gf_fymomx(Ipy), gf_fzmomx(Imz), gf_fzmomx(Ipz));
        gf_rhsmomy(p.I) =
            calcupdate(gf_fxmomy(Imx), gf_fxmomy(Ipx), gf_fymomy(Imy),
                       gf_fymomy(Ipy), gf_fzmomy(Imz), gf_fzmomy(Ipz));
        gf_rhsmomz(p.I) =
            calcupdate(gf_fxmomz(Imx), gf_fxmomz(Ipx), gf_fymomz(Imy),
                       gf_fymomz(Ipy), gf_fzmomz(Imz), gf_fzmomz(Ipz));
        gf_rhstau(p.I) =
            calcupdate(gf_fxtau(Imx), gf_fxtau(Ipx), gf_fytau(Imy),
                       gf_fytau(Ipy), gf_fztau(Imz), gf_fztau(Ipz));
      });

#ifdef AMREX_USE_GPU
  AMREX_GPU_ERROR_CHECK();
  amrex::Gpu::synchronize();
  AMREX_GPU_ERROR_CHECK();
#endif

#else

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE CCTK_HOST(
                            const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Neighbouring "plus" and "minus" face indices in
        // the x, y, and z directions
        const auto Imx = p.I;
        const auto Imy = p.I;
        const auto Imz = p.I;
        const auto Ipx = p.I + DI[0];
        const auto Ipy = p.I + DI[1];
        const auto Ipz = p.I + DI[2];

        gf_densrhs(p.I) =
            calcupdate(gf_fxdens(Imx), gf_fxdens(Ipx), gf_fydens(Imy),
                       gf_fydens(Ipy), gf_fzdens(Imz), gf_fzdens(Ipz));
      });

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE CCTK_HOST(
                            const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Neighbouring "plus" and "minus" face indices in
        // the x, y, and z directions
        const auto Imx = p.I;
        const auto Imy = p.I;
        const auto Imz = p.I;
        const auto Ipx = p.I + DI[0];
        const auto Ipy = p.I + DI[1];
        const auto Ipz = p.I + DI[2];

        gf_momxrhs(p.I) =
            calcupdate(gf_fxmomx(Imx), gf_fxmomx(Ipx), gf_fymomx(Imy),
                       gf_fymomx(Ipy), gf_fzmomx(Imz), gf_fzmomx(Ipz));
      });

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE CCTK_HOST(
                            const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Neighbouring "plus" and "minus" face indices in
        // the x, y, and z directions
        const auto Imx = p.I;
        const auto Imy = p.I;
        const auto Imz = p.I;
        const auto Ipx = p.I + DI[0];
        const auto Ipy = p.I + DI[1];
        const auto Ipz = p.I + DI[2];

        gf_momyrhs(p.I) =
            calcupdate(gf_fxmomy(Imx), gf_fxmomy(Ipx), gf_fymomy(Imy),
                       gf_fymomy(Ipy), gf_fzmomy(Imz), gf_fzmomy(Ipz));
      });

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE CCTK_HOST(
                            const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Neighbouring "plus" and "minus" face indices in
        // the x, y, and z directions
        const auto Imx = p.I;
        const auto Imy = p.I;
        const auto Imz = p.I;
        const auto Ipx = p.I + DI[0];
        const auto Ipy = p.I + DI[1];
        const auto Ipz = p.I + DI[2];

        gf_momzrhs(p.I) =
            calcupdate(gf_fxmomz(Imx), gf_fxmomz(Ipx), gf_fymomz(Imy),
                       gf_fymomz(Ipy), gf_fzmomz(Imz), gf_fzmomz(Ipz));
      });

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE CCTK_HOST(
                            const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Neighbouring "plus" and "minus" face indices in
        // the x, y, and z directions
        const auto Imx = p.I;
        const auto Imy = p.I;
        const auto Imz = p.I;
        const auto Ipx = p.I + DI[0];
        const auto Ipy = p.I + DI[1];
        const auto Ipz = p.I + DI[2];

        gf_taurhs(p.I) =
            calcupdate(gf_fxtau(Imx), gf_fxtau(Ipx), gf_fytau(Imy),
                       gf_fytau(Ipy), gf_fztau(Imz), gf_fztau(Ipz));
      });

#endif
}

} // namespace GRHydroToyGPU
