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

extern "C" void Z4c_Constraints(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_Constraints;
  DECLARE_CCTK_PARAMETERS;

  for (int d = 0; d < 3; ++d)
    if (cctk_nghostzones[d] < deriv_order / 2 + 1)
      CCTK_VERROR("Need at least %d ghost zones", deriv_order / 2 + 1);

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

  //

  const GF3D<CCTK_REAL, 0, 0, 0> gf_ZtCx_(cctkGH, ZtCx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_ZtCy_(cctkGH, ZtCy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_ZtCz_(cctkGH, ZtCz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_HC_(cctkGH, HC);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_MtCx_(cctkGH, MtCx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_MtCy_(cctkGH, MtCy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_MtCz_(cctkGH, MtCz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_allC_(cctkGH, allC);

  //

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    // Load and calculate

    const CCTK_REAL alphaG{1};
    const vec3<CCTK_REAL, DN> dalphaG{0, 0, 0};
    const mat3<CCTK_REAL, DN, DN> ddalphaG{0, 0, 0, 0, 0, 0};

    const vec3<CCTK_REAL, UP> betaG{0, 0, 0};
    const vec3<vec3<CCTK_REAL, DN>, UP> dbetaG{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    const vec3<mat3<CCTK_REAL, DN, DN>, UP> ddbetaG{
        {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};

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
        alphaG, dalphaG, ddalphaG,                            //
        betaG, dbetaG, ddbetaG,                               //
        rho,                                                  //
        Si,                                                   //
        Sij);

    // Store
    vars.ZtC.store(gf_ZtCx_, gf_ZtCy_, gf_ZtCz_, p.I);
    gf_HC_(p.I) = vars.HC;
    vars.MtC.store(gf_MtCx_, gf_MtCy_, gf_MtCz_, p.I);
    gf_allC_(p.I) = vars.allC;
  });
}

} // namespace Z4c
