#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "utils.hxx"

// #ifdef AMREX_USE_GPU
// #include <AMReX_GpuDevice.H>
// #endif

#include <array>

namespace AsterX {
using namespace Loop;

extern "C" void AsterX_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_RHS;
  DECLARE_CCTK_PARAMETERS;

  const std::array dx{CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1),
                      CCTK_DELTA_SPACE(2)};
  const std::array idx{1 / dx[0], 1 / dx[1], 1 / dx[2]};

  const auto calcupdate =
      [=] CCTK_DEVICE(CCTK_REAL fx_m, CCTK_REAL fx_p, CCTK_REAL fy_m,
                      CCTK_REAL fy_p, CCTK_REAL fz_m, CCTK_REAL fz_p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return -idx[0] * (fx_p - fx_m) - idx[1] * (fy_p - fy_m) -
                   idx[2] * (fz_p - fz_m);
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

        // TODO: add lapse terms, and check calcupdate code
        densrhs(p.I) =
            densrhs(p.I) + calcupdate(fxdens(Imx), fxdens(Ipx), fydens(Imy),
                                      fydens(Ipy), fzdens(Imz), fzdens(Ipz));
        momxrhs(p.I) =
            momxrhs(p.I) + calcupdate(fxmomx(Imx), fxmomx(Ipx), fymomx(Imy),
                                      fymomx(Ipy), fzmomx(Imz), fzmomx(Ipz));
        momyrhs(p.I) =
            momyrhs(p.I) + calcupdate(fxmomy(Imx), fxmomy(Ipx), fymomy(Imy),
                                      fymomy(Ipy), fzmomy(Imz), fzmomy(Ipz));
        momzrhs(p.I) =
            momzrhs(p.I) + calcupdate(fxmomz(Imx), fxmomz(Ipx), fymomz(Imy),
                                      fymomz(Ipy), fzmomz(Imz), fzmomz(Ipz));
        taurhs(p.I) =
            taurhs(p.I) + calcupdate(fxtau(Imx), fxtau(Ipx), fytau(Imy),
                                     fytau(Ipy), fztau(Imz), fztau(Ipz));
      });

  grid.loop_int_device<1, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const auto I1m = p.I - DI[1];
        const auto I2m = p.I - DI[2];
        const auto Ex = 0.25*(- fyBz(p.I) - fyBz(I2m) + fzBy(p.I) + fzBy(I1m));
        Avec_x_rhs(p.I) = -Ex - calc_fd2_v2e(G, p, 0);
      });

  grid.loop_int_device<0, 1, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const auto I2m = p.I - DI[2];
        const auto I0m = p.I - DI[0];
        const auto Ey = 0.25*(- fzBx(p.I) - fzBx(I0m) + fxBz(p.I) + fxBz(I2m));
        Avec_y_rhs(p.I) = -Ey - calc_fd2_v2e(G, p, 1);
      });

  grid.loop_int_device<0, 0, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const auto I0m = p.I - DI[0];
        const auto I1m = p.I - DI[1];
        const auto Ez = 0.25*(- fxBy(p.I) - fxBy(I1m) + fyBx(p.I) + fyBx(I0m));
        Avec_z_rhs(p.I) = -Ez - calc_fd2_v2e(G, p, 2);
      });

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        /* diFi on vertices */
        const CCTK_REAL dF =
            calc_fd2_c2c(Fx, p, 0) + calc_fd2_c2c(Fy, p, 1) +
            calc_fd2_c2c(Fz, p, 2); // should be v2v but c2c works too
        /* diFbetai */
        CCTK_REAL dFbeta = 0.0;
        if (betax(p.I) < 0) {
          dFbeta += calc_fd2_v2v_oneside(Fbetax, p, 0, -1);
        } else {
          dFbeta += calc_fd2_v2v_oneside(Fbetax, p, 0, 1);
        }
        if (betay(p.I) < 0) {
          dFbeta += calc_fd2_v2v_oneside(Fbetay, p, 1, -1);
        } else {
          dFbeta += calc_fd2_v2v_oneside(Fbetay, p, 1, 1);
        }
        if (betaz(p.I) < 0) {
          dFbeta += calc_fd2_v2v_oneside(Fbetaz, p, 2, -1);
        } else {
          dFbeta += calc_fd2_v2v_oneside(Fbetaz, p, 2, 1);
        }

        const CCTK_REAL lorentz_damp_fac = 0.0;
        /* rhs of Psi */
        Psi_rhs(p.I) = -dF + dFbeta - lorentz_damp_fac * alp(p.I) * Psi(p.I);
      });
}

} // namespace AsterX
