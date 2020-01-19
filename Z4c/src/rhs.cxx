#include "physics.hxx"
#include "tensor.hxx"
#include "z4c_vars.hxx"

#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace Z4c {
using namespace Loop;
using namespace std;

template <typename T>
void apply_upwind(const cGH *restrict const cctkGH,
                  const GF3D<const T, 0, 0, 0> gf_,
                  const GF3D<const T, 0, 0, 0> gf_betaGx_,
                  const GF3D<const T, 0, 0, 0> gf_betaGy_,
                  const GF3D<const T, 0, 0, 0> gf_betaGz_,
                  const GF3D<T, 0, 0, 0> gf_rhs_) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const vec3<CCTK_REAL> dx{
      CCTK_DELTA_SPACE(0),
      CCTK_DELTA_SPACE(1),
      CCTK_DELTA_SPACE(2),
  };

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    const vec3<CCTK_REAL> betaG(gf_betaGx_, gf_betaGy_, gf_betaGz_, p.I);
    const vec3<CCTK_REAL> dgf_upwind(deriv_upwind(gf_, p.I, betaG, dx));
    gf_rhs_(p.I) += sum1([&](int x) { return betaG(x) * dgf_upwind(x); });
  });
}

template <typename T>
void apply_diss(const cGH *restrict const cctkGH,
                const GF3D<const T, 0, 0, 0> gf_,
                const GF3D<T, 0, 0, 0> gf_rhs_) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const vec3<CCTK_REAL> dx{
      CCTK_DELTA_SPACE(0),
      CCTK_DELTA_SPACE(1),
      CCTK_DELTA_SPACE(2),
  };

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    gf_rhs_(p.I) += epsdiss * diss(gf_, p.I, dx);
  });
}

extern "C" void Z4c_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_RHS;
  DECLARE_CCTK_PARAMETERS;

  for (int d = 0; d < 3; ++d)
    if (cctk_nghostzones[d] < 2)
      CCTK_VERROR("Need at least 2 ghost zones");

  const vec3<CCTK_REAL> dx{
      CCTK_DELTA_SPACE(0),
      CCTK_DELTA_SPACE(1),
      CCTK_DELTA_SPACE(2),
  };

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_chi_(cctkGH, chi);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatxx_(cctkGH, gammatxx);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatxy_(cctkGH, gammatxy);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatxz_(cctkGH, gammatxz);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatyy_(cctkGH, gammatyy);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatyz_(cctkGH, gammatyz);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatzz_(cctkGH, gammatzz);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Kh_(cctkGH, Kh);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Atxx_(cctkGH, Atxx);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Atxy_(cctkGH, Atxy);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Atxz_(cctkGH, Atxz);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Atyy_(cctkGH, Atyy);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Atyz_(cctkGH, Atyz);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Atzz_(cctkGH, Atzz);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Gamtx_(cctkGH, Gamtx);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Gamty_(cctkGH, Gamty);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Gamtz_(cctkGH, Gamtz);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Theta_(cctkGH, Theta);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_alphaG_(cctkGH, alphaG);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_betaGx_(cctkGH, betaGx);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_betaGy_(cctkGH, betaGy);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_betaGz_(cctkGH, betaGz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_chi_rhs_(cctkGH, chi_rhs);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatxx_rhs_(cctkGH, gammatxx_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatxy_rhs_(cctkGH, gammatxy_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatxz_rhs_(cctkGH, gammatxz_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatyy_rhs_(cctkGH, gammatyy_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatyz_rhs_(cctkGH, gammatyz_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatzz_rhs_(cctkGH, gammatzz_rhs);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_Kh_rhs_(cctkGH, Kh_rhs);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atxx_rhs_(cctkGH, Atxx_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atxy_rhs_(cctkGH, Atxy_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atxz_rhs_(cctkGH, Atxz_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atyy_rhs_(cctkGH, Atyy_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atyz_rhs_(cctkGH, Atyz_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atzz_rhs_(cctkGH, Atzz_rhs);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamtx_rhs_(cctkGH, Gamtx_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamty_rhs_(cctkGH, Gamty_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamtz_rhs_(cctkGH, Gamtz_rhs);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_Theta_rhs_(cctkGH, Theta_rhs);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_alphaG_rhs_(cctkGH, alphaG_rhs);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_betaGx_rhs_(cctkGH, betaGx_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_betaGy_rhs_(cctkGH, betaGy_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_betaGz_rhs_(cctkGH, betaGz_rhs);

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    // Load and calculate
    const z4c_vars<CCTK_REAL> vars(
        kappa1, kappa2, f_mu_L, f_mu_S, eta, //
        gf_chi_, gf_gammatxx_, gf_gammatxy_, gf_gammatxz_, gf_gammatyy_,
        gf_gammatyz_, gf_gammatzz_, gf_Kh_, gf_Atxx_, gf_Atxy_, gf_Atxz_,
        gf_Atyy_, gf_Atyz_, gf_Atzz_, gf_Gamtx_, gf_Gamty_, gf_Gamtz_,
        gf_Theta_, gf_alphaG_, gf_betaGx_, gf_betaGy_, gf_betaGz_, //
        p.I, dx);

    // Store
    gf_chi_rhs_(p.I) = vars.chi_rhs;
    vars.gammat_rhs.store(gf_gammatxx_rhs_, gf_gammatxy_rhs_, gf_gammatxz_rhs_,
                          gf_gammatyy_rhs_, gf_gammatyz_rhs_, gf_gammatzz_rhs_,
                          p);
    gf_Kh_rhs_(p.I) = vars.Kh_rhs;
    vars.At_rhs.store(gf_Atxx_rhs_, gf_Atxy_rhs_, gf_Atxz_rhs_, gf_Atyy_rhs_,
                      gf_Atyz_rhs_, gf_Atzz_rhs_, p);
    vars.Gamt_rhs.store(gf_Gamtx_rhs_, gf_Gamty_rhs_, gf_Gamtz_rhs_, p);
    gf_Theta_rhs_(p.I) = vars.Theta_rhs;
    gf_alphaG_rhs_(p.I) = vars.alphaG_rhs;
    vars.betaG_rhs.store(gf_betaGx_rhs_, gf_betaGy_rhs_, gf_betaGz_rhs_, p);
  });

  // Upwind terms

  apply_upwind(cctkGH, gf_chi_, gf_betaGx_, gf_betaGy_, gf_betaGz_,
               gf_chi_rhs_);

  apply_upwind(cctkGH, gf_gammatxx_, gf_betaGx_, gf_betaGy_, gf_betaGz_,
               gf_gammatxx_rhs_);
  apply_upwind(cctkGH, gf_gammatxy_, gf_betaGx_, gf_betaGy_, gf_betaGz_,
               gf_gammatxy_rhs_);
  apply_upwind(cctkGH, gf_gammatxz_, gf_betaGx_, gf_betaGy_, gf_betaGz_,
               gf_gammatxz_rhs_);
  apply_upwind(cctkGH, gf_gammatyy_, gf_betaGx_, gf_betaGy_, gf_betaGz_,
               gf_gammatyy_rhs_);
  apply_upwind(cctkGH, gf_gammatyz_, gf_betaGx_, gf_betaGy_, gf_betaGz_,
               gf_gammatyz_rhs_);
  apply_upwind(cctkGH, gf_gammatzz_, gf_betaGx_, gf_betaGy_, gf_betaGz_,
               gf_gammatzz_rhs_);

  apply_upwind(cctkGH, gf_Kh_, gf_betaGx_, gf_betaGy_, gf_betaGz_, gf_Kh_rhs_);

  apply_upwind(cctkGH, gf_Atxx_, gf_betaGx_, gf_betaGy_, gf_betaGz_,
               gf_Atxx_rhs_);
  apply_upwind(cctkGH, gf_Atxy_, gf_betaGx_, gf_betaGy_, gf_betaGz_,
               gf_Atxy_rhs_);
  apply_upwind(cctkGH, gf_Atxz_, gf_betaGx_, gf_betaGy_, gf_betaGz_,
               gf_Atxz_rhs_);
  apply_upwind(cctkGH, gf_Atyy_, gf_betaGx_, gf_betaGy_, gf_betaGz_,
               gf_Atyy_rhs_);
  apply_upwind(cctkGH, gf_Atyz_, gf_betaGx_, gf_betaGy_, gf_betaGz_,
               gf_Atyz_rhs_);
  apply_upwind(cctkGH, gf_Atzz_, gf_betaGx_, gf_betaGy_, gf_betaGz_,
               gf_Atzz_rhs_);

  apply_upwind(cctkGH, gf_Gamtx_, gf_betaGx_, gf_betaGy_, gf_betaGz_,
               gf_Gamtx_rhs_);
  apply_upwind(cctkGH, gf_Gamty_, gf_betaGx_, gf_betaGy_, gf_betaGz_,
               gf_Gamty_rhs_);
  apply_upwind(cctkGH, gf_Gamtz_, gf_betaGx_, gf_betaGy_, gf_betaGz_,
               gf_Gamtz_rhs_);

  apply_upwind(cctkGH, gf_Theta_, gf_betaGx_, gf_betaGy_, gf_betaGz_,
               gf_Theta_rhs_);

  apply_upwind(cctkGH, gf_alphaG_, gf_betaGx_, gf_betaGy_, gf_betaGz_,
               gf_alphaG_rhs_);

  apply_upwind(cctkGH, gf_betaGx_, gf_betaGx_, gf_betaGy_, gf_betaGz_,
               gf_betaGx_rhs_);
  apply_upwind(cctkGH, gf_betaGy_, gf_betaGx_, gf_betaGy_, gf_betaGz_,
               gf_betaGy_rhs_);
  apply_upwind(cctkGH, gf_betaGz_, gf_betaGx_, gf_betaGy_, gf_betaGz_,
               gf_betaGz_rhs_);

  // Dissipation

  apply_diss(cctkGH, gf_chi_, gf_chi_rhs_);

  apply_diss(cctkGH, gf_gammatxx_, gf_gammatxx_rhs_);
  apply_diss(cctkGH, gf_gammatxy_, gf_gammatxy_rhs_);
  apply_diss(cctkGH, gf_gammatxz_, gf_gammatxz_rhs_);
  apply_diss(cctkGH, gf_gammatyy_, gf_gammatyy_rhs_);
  apply_diss(cctkGH, gf_gammatyz_, gf_gammatyz_rhs_);
  apply_diss(cctkGH, gf_gammatzz_, gf_gammatzz_rhs_);

  apply_diss(cctkGH, gf_Kh_, gf_Kh_rhs_);

  apply_diss(cctkGH, gf_Atxx_, gf_Atxx_rhs_);
  apply_diss(cctkGH, gf_Atxy_, gf_Atxy_rhs_);
  apply_diss(cctkGH, gf_Atxz_, gf_Atxz_rhs_);
  apply_diss(cctkGH, gf_Atyy_, gf_Atyy_rhs_);
  apply_diss(cctkGH, gf_Atyz_, gf_Atyz_rhs_);
  apply_diss(cctkGH, gf_Atzz_, gf_Atzz_rhs_);

  apply_diss(cctkGH, gf_Gamtx_, gf_Gamtx_rhs_);
  apply_diss(cctkGH, gf_Gamty_, gf_Gamty_rhs_);
  apply_diss(cctkGH, gf_Gamtz_, gf_Gamtz_rhs_);

  apply_diss(cctkGH, gf_Theta_, gf_Theta_rhs_);

  apply_diss(cctkGH, gf_alphaG_, gf_alphaG_rhs_);

  apply_diss(cctkGH, gf_betaGx_, gf_betaGx_rhs_);
  apply_diss(cctkGH, gf_betaGy_, gf_betaGy_rhs_);
  apply_diss(cctkGH, gf_betaGz_, gf_betaGz_rhs_);
}

} // namespace Z4c
