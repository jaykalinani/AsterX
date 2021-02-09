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

extern "C" void Z4c_ADM2(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_ADM2;
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

  vector<vector<CCTK_REAL> > buffers;

  const vec3<GF3D<CCTK_REAL, 0, 0, 0>, DN> gf_dchi_(cctkGH, buffers);
  const mat3<GF3D<CCTK_REAL, 0, 0, 0>, DN, DN> gf_ddchi_(cctkGH, buffers);
  calc_derivs2(cctkGH, gf_chi_, gf_dchi_, gf_ddchi_);

  const mat3<vec3<GF3D<CCTK_REAL, 0, 0, 0>, DN>, DN, DN> gf_dgammat_(cctkGH,
                                                                     buffers);
  const mat3<mat3<GF3D<CCTK_REAL, 0, 0, 0>, DN, DN>, DN, DN> gf_ddgammat_(
      cctkGH, buffers);
  calc_derivs2(cctkGH, gf_gammat_, gf_dgammat_, gf_ddgammat_);

  const vec3<GF3D<CCTK_REAL, 0, 0, 0>, DN> gf_dKh_(cctkGH, buffers);
  calc_derivs(cctkGH, gf_Kh_, gf_dKh_);

  const mat3<vec3<GF3D<CCTK_REAL, 0, 0, 0>, DN>, DN, DN> gf_dAt_(cctkGH,
                                                                 buffers);
  calc_derivs(cctkGH, gf_At_, gf_dAt_);

  const vec3<vec3<GF3D<CCTK_REAL, 0, 0, 0>, DN>, UP> gf_dGamt_(cctkGH, buffers);
  calc_derivs(cctkGH, gf_Gamt_, gf_dGamt_);

  const vec3<GF3D<CCTK_REAL, 0, 0, 0>, DN> gf_dTheta_(cctkGH, buffers);
  calc_derivs(cctkGH, gf_Theta_, gf_dTheta_);

  const vec3<GF3D<CCTK_REAL, 0, 0, 0>, DN> gf_dalphaG_(cctkGH, buffers);
  const mat3<GF3D<CCTK_REAL, 0, 0, 0>, DN, DN> gf_ddalphaG_(cctkGH, buffers);
  calc_derivs2(cctkGH, gf_alphaG_, gf_dalphaG_, gf_ddalphaG_);

  const vec3<vec3<GF3D<CCTK_REAL, 0, 0, 0>, DN>, UP> gf_dbetaG_(cctkGH,
                                                                buffers);
  const vec3<mat3<GF3D<CCTK_REAL, 0, 0, 0>, DN, DN>, UP> gf_ddbetaG_(cctkGH,
                                                                     buffers);
  calc_derivs2(cctkGH, gf_betaG_, gf_dbetaG_, gf_ddbetaG_);

  //

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
    vars.K_rhs.store(gf_dtkxx_, gf_dtkxy_, gf_dtkxz_, gf_dtkyy_, gf_dtkyz_,
                     gf_dtkzz_, p.I);
    gf_dt2alp_(p.I) = vars.dtalpha_rhs;
    vars.dtbeta_rhs.store(gf_dt2betax_, gf_dt2betay_, gf_dt2betaz_, p.I);
  });
}

} // namespace Z4c
