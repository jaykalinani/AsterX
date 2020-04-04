#include "derivs.hxx"
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

extern "C" void Z4c_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_RHS;
  DECLARE_CCTK_PARAMETERS;

  for (int d = 0; d < 3; ++d)
    if (cctk_nghostzones[d] < deriv_order / 2 + 1)
      CCTK_VERROR("Need at least %d ghost zones", deriv_order / 2 + 1);

  const vec3<CCTK_REAL, UP> dx{
      CCTK_DELTA_SPACE(0),
      CCTK_DELTA_SPACE(1),
      CCTK_DELTA_SPACE(2),
  };

  //

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_chi_(cctkGH, chi);

  const mat3<GF3D<const CCTK_REAL, 0, 0, 0>, DN, DN> gf_gammat_(
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, gammatxx),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, gammatxy),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, gammatxz),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, gammatyy),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, gammatyz),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, gammatzz));

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Kh_(cctkGH, Kh);

  const mat3<GF3D<const CCTK_REAL, 0, 0, 0>, DN, DN> gf_At_(
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, Atxx),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, Atxy),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, Atxz),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, Atyy),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, Atyz),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, Atzz));

  const vec3<GF3D<const CCTK_REAL, 0, 0, 0>, UP> gf_Gamt_(
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, Gamtx),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, Gamty),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, Gamtz));

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Theta_(cctkGH, Theta);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_alphaG_(cctkGH, alphaG);

  const vec3<GF3D<const CCTK_REAL, 0, 0, 0>, UP> gf_betaG_(
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, betaGx),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, betaGy),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, betaGz));

  //

  const vec3<GF3D<CCTK_REAL, 0, 0, 0>, DN> gf_dchi_(cctkGH, allocate());
  const mat3<GF3D<CCTK_REAL, 0, 0, 0>, DN, DN> gf_ddchi_(cctkGH, allocate());
  calc_derivs2(cctkGH, gf_chi_, gf_dchi_, gf_ddchi_);

  const mat3<vec3<GF3D<CCTK_REAL, 0, 0, 0>, DN>, DN, DN> gf_dgammat_(
      cctkGH, allocate());
  const mat3<mat3<GF3D<CCTK_REAL, 0, 0, 0>, DN, DN>, DN, DN> gf_ddgammat_(
      cctkGH, allocate());
  calc_derivs2(cctkGH, gf_gammat_, gf_dgammat_, gf_ddgammat_);

  const vec3<GF3D<CCTK_REAL, 0, 0, 0>, DN> gf_dKh_(cctkGH, allocate());
  calc_derivs(cctkGH, gf_Kh_, gf_dKh_);

  const mat3<vec3<GF3D<CCTK_REAL, 0, 0, 0>, DN>, DN, DN> gf_dAt_(cctkGH,
                                                                 allocate());
  calc_derivs(cctkGH, gf_At_, gf_dAt_);

  const vec3<vec3<GF3D<CCTK_REAL, 0, 0, 0>, DN>, UP> gf_dGamt_(cctkGH,
                                                               allocate());
  calc_derivs(cctkGH, gf_Gamt_, gf_dGamt_);

  const vec3<GF3D<CCTK_REAL, 0, 0, 0>, DN> gf_dTheta_(cctkGH, allocate());
  calc_derivs(cctkGH, gf_Theta_, gf_dTheta_);

  const vec3<GF3D<CCTK_REAL, 0, 0, 0>, DN> gf_dalphaG_(cctkGH, allocate());
  const mat3<GF3D<CCTK_REAL, 0, 0, 0>, DN, DN> gf_ddalphaG_(cctkGH, allocate());
  calc_derivs2(cctkGH, gf_alphaG_, gf_dalphaG_, gf_ddalphaG_);

  const vec3<vec3<GF3D<CCTK_REAL, 0, 0, 0>, DN>, UP> gf_dbetaG_(cctkGH,
                                                                allocate());
  const vec3<mat3<GF3D<CCTK_REAL, 0, 0, 0>, DN, DN>, UP> gf_ddbetaG_(
      cctkGH, allocate());
  calc_derivs2(cctkGH, gf_betaG_, gf_dbetaG_, gf_ddbetaG_);

  //

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

  const GF3D<CCTK_REAL, 0, 0, 0> gf_dtkxx_(cctkGH, dtkxx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_dtkxy_(cctkGH, dtkxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_dtkxz_(cctkGH, dtkxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_dtkyy_(cctkGH, dtkyy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_dtkyz_(cctkGH, dtkyz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_dtkzz_(cctkGH, dtkzz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_dt2alp_(cctkGH, dt2alp);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_dt2betax_(cctkGH, dt2betax);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_dt2betay_(cctkGH, dt2betay);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_dt2betaz_(cctkGH, dt2betaz);

  //

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) Z4C_INLINE {
    // Load and calculate

    const CCTK_REAL rho{0};
    const vec3<CCTK_REAL, DN> Si{0, 0, 0};
    const mat3<CCTK_REAL, DN, DN> Sij{0, 0, 0, 0, 0, 0};

    const z4c_vars<CCTK_REAL> vars(
        kappa1, kappa2, f_mu_L, f_mu_S, eta,                  //
        gf_chi_(p.I), gf_dchi_(p.I), gf_ddchi_(p.I),          //
        gf_gammat_(p.I), gf_dgammat_(p.I), gf_ddgammat_(p.I), //
        gf_Kh_(p.I), gf_dKh_(p.I),                            //
        gf_At_(p.I), gf_dAt_(p.I),                            //
        gf_Gamt_(p.I), gf_dGamt_(p.I),                        //
        gf_Theta_(p.I), gf_dTheta_(p.I),                      //
        gf_alphaG_(p.I), gf_dalphaG_(p.I), gf_ddalphaG_(p.I), //
        gf_betaG_(p.I), gf_dbetaG_(p.I), gf_ddbetaG_(p.I),    //
        rho,                                                  //
        Si,                                                   //
        Sij);

    // Store
    gf_chi_rhs_(p.I) = vars.chi_rhs;
    vars.gammat_rhs.store(gf_gammatxx_rhs_, gf_gammatxy_rhs_, gf_gammatxz_rhs_,
                          gf_gammatyy_rhs_, gf_gammatyz_rhs_, gf_gammatzz_rhs_,
                          p.I);
    gf_Kh_rhs_(p.I) = vars.Kh_rhs;
    vars.At_rhs.store(gf_Atxx_rhs_, gf_Atxy_rhs_, gf_Atxz_rhs_, gf_Atyy_rhs_,
                      gf_Atyz_rhs_, gf_Atzz_rhs_, p.I);
    vars.Gamt_rhs.store(gf_Gamtx_rhs_, gf_Gamty_rhs_, gf_Gamtz_rhs_, p.I);
    gf_Theta_rhs_(p.I) = vars.Theta_rhs;
    gf_alphaG_rhs_(p.I) = vars.alphaG_rhs;
    vars.betaG_rhs.store(gf_betaGx_rhs_, gf_betaGy_rhs_, gf_betaGz_rhs_, p.I);

    vars.K_rhs.store(gf_dtkxx_, gf_dtkxy_, gf_dtkxz_, gf_dtkyy_, gf_dtkyz_,
                     gf_dtkzz_, p.I);
    gf_dt2alp_(p.I) = vars.dtalpha_rhs;
    vars.dtbeta_rhs.store(gf_dt2betax_, gf_dt2betay_, gf_dt2betaz_, p.I);
  });

  // Upwind terms

  apply_upwind(cctkGH, gf_chi_, gf_betaG_, gf_chi_rhs_);

  apply_upwind(cctkGH, gf_gammat_(0, 0), gf_betaG_, gf_gammatxx_rhs_);
  apply_upwind(cctkGH, gf_gammat_(0, 1), gf_betaG_, gf_gammatxy_rhs_);
  apply_upwind(cctkGH, gf_gammat_(0, 2), gf_betaG_, gf_gammatxz_rhs_);
  apply_upwind(cctkGH, gf_gammat_(1, 1), gf_betaG_, gf_gammatyy_rhs_);
  apply_upwind(cctkGH, gf_gammat_(1, 2), gf_betaG_, gf_gammatyz_rhs_);
  apply_upwind(cctkGH, gf_gammat_(2, 2), gf_betaG_, gf_gammatzz_rhs_);

  apply_upwind(cctkGH, gf_Kh_, gf_betaG_, gf_Kh_rhs_);

  apply_upwind(cctkGH, gf_At_(0, 0), gf_betaG_, gf_Atxx_rhs_);
  apply_upwind(cctkGH, gf_At_(0, 1), gf_betaG_, gf_Atxy_rhs_);
  apply_upwind(cctkGH, gf_At_(0, 2), gf_betaG_, gf_Atxz_rhs_);
  apply_upwind(cctkGH, gf_At_(1, 1), gf_betaG_, gf_Atyy_rhs_);
  apply_upwind(cctkGH, gf_At_(1, 2), gf_betaG_, gf_Atyz_rhs_);
  apply_upwind(cctkGH, gf_At_(2, 2), gf_betaG_, gf_Atzz_rhs_);

  apply_upwind(cctkGH, gf_Gamt_(0), gf_betaG_, gf_Gamtx_rhs_);
  apply_upwind(cctkGH, gf_Gamt_(1), gf_betaG_, gf_Gamty_rhs_);
  apply_upwind(cctkGH, gf_Gamt_(2), gf_betaG_, gf_Gamtz_rhs_);

  apply_upwind(cctkGH, gf_Theta_, gf_betaG_, gf_Theta_rhs_);

  apply_upwind(cctkGH, gf_alphaG_, gf_betaG_, gf_alphaG_rhs_);

  apply_upwind(cctkGH, gf_betaG_(0), gf_betaG_, gf_betaGx_rhs_);
  apply_upwind(cctkGH, gf_betaG_(1), gf_betaG_, gf_betaGy_rhs_);
  apply_upwind(cctkGH, gf_betaG_(2), gf_betaG_, gf_betaGz_rhs_);

  // Dissipation

  apply_diss(cctkGH, gf_chi_, gf_chi_rhs_);

  apply_diss(cctkGH, gf_gammat_(0, 0), gf_gammatxx_rhs_);
  apply_diss(cctkGH, gf_gammat_(0, 1), gf_gammatxy_rhs_);
  apply_diss(cctkGH, gf_gammat_(0, 2), gf_gammatxz_rhs_);
  apply_diss(cctkGH, gf_gammat_(1, 1), gf_gammatyy_rhs_);
  apply_diss(cctkGH, gf_gammat_(1, 2), gf_gammatyz_rhs_);
  apply_diss(cctkGH, gf_gammat_(2, 2), gf_gammatzz_rhs_);

  apply_diss(cctkGH, gf_Kh_, gf_Kh_rhs_);

  apply_diss(cctkGH, gf_At_(0, 0), gf_Atxx_rhs_);
  apply_diss(cctkGH, gf_At_(0, 1), gf_Atxy_rhs_);
  apply_diss(cctkGH, gf_At_(0, 2), gf_Atxz_rhs_);
  apply_diss(cctkGH, gf_At_(1, 1), gf_Atyy_rhs_);
  apply_diss(cctkGH, gf_At_(1, 2), gf_Atyz_rhs_);
  apply_diss(cctkGH, gf_At_(2, 2), gf_Atzz_rhs_);

  apply_diss(cctkGH, gf_Gamt_(0), gf_Gamtx_rhs_);
  apply_diss(cctkGH, gf_Gamt_(1), gf_Gamty_rhs_);
  apply_diss(cctkGH, gf_Gamt_(2), gf_Gamtz_rhs_);

  apply_diss(cctkGH, gf_Theta_, gf_Theta_rhs_);

  apply_diss(cctkGH, gf_alphaG_, gf_alphaG_rhs_);

  apply_diss(cctkGH, gf_betaG_(0), gf_betaGx_rhs_);
  apply_diss(cctkGH, gf_betaG_(1), gf_betaGy_rhs_);
  apply_diss(cctkGH, gf_betaG_(2), gf_betaGz_rhs_);
}

} // namespace Z4c
