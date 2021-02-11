#include "derivs.hxx"
#include "physics.hxx"
#include "tensor.hxx"
#include "z4c_vars.hxx"

#include <loop.hxx>
#include <mempool.hxx>

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

  const array<int, dim> indextype = {0, 0, 0};
  const array<int, dim> noghosts = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);
  const GF3D2layout layout0(cctkGH, indextype, noghosts);

  const GF3D2<const CCTK_REAL> gf_chi0_(&layout, chi);

  const mat3<GF3D2<const CCTK_REAL>, DN, DN> gf_gammat0_(
      GF3D2<const CCTK_REAL>(&layout, gammatxx),
      GF3D2<const CCTK_REAL>(&layout, gammatxy),
      GF3D2<const CCTK_REAL>(&layout, gammatxz),
      GF3D2<const CCTK_REAL>(&layout, gammatyy),
      GF3D2<const CCTK_REAL>(&layout, gammatyz),
      GF3D2<const CCTK_REAL>(&layout, gammatzz));

  const GF3D2<const CCTK_REAL> gf_Kh0_(&layout, Kh);

  const mat3<GF3D2<const CCTK_REAL>, DN, DN> gf_At0_(
      GF3D2<const CCTK_REAL>(&layout, Atxx),
      GF3D2<const CCTK_REAL>(&layout, Atxy),
      GF3D2<const CCTK_REAL>(&layout, Atxz),
      GF3D2<const CCTK_REAL>(&layout, Atyy),
      GF3D2<const CCTK_REAL>(&layout, Atyz),
      GF3D2<const CCTK_REAL>(&layout, Atzz));

  const vec3<GF3D2<const CCTK_REAL>, UP> gf_Gamt0_(
      GF3D2<const CCTK_REAL>(&layout, Gamtx),
      GF3D2<const CCTK_REAL>(&layout, Gamty),
      GF3D2<const CCTK_REAL>(&layout, Gamtz));

  const GF3D2<const CCTK_REAL> gf_Theta0_(&layout, Theta);

  //

  static mempool_set_t mempools;
  mempool_t &restrict mempool = mempools.get_mempool();

  const auto make_gf = [&]() { return GF3D2<CCTK_REAL>(&layout0, mempool); };
  const auto make_vec_gf = [&](int) { return make_gf(); };
  const auto make_mat_gf = [&](int, int) { return make_gf(); };
  const auto make_vec_vec_gf = [&](int) {
    return [&](int) { return make_gf(); };
  };
  const auto make_mat_vec_gf = [&](int, int) {
    return [&](int) { return make_gf(); };
  };
  const auto make_mat_mat_gf = [&](int, int) {
    return [&](int, int) { return make_gf(); };
  };

  const GF3D2<CCTK_REAL> gf_chi_(make_gf());
  const vec3<GF3D2<CCTK_REAL>, DN> gf_dchi_(make_vec_gf);
  const mat3<GF3D2<CCTK_REAL>, DN, DN> gf_ddchi_(make_mat_gf);
  calc_derivs2(cctkGH, gf_chi0_, gf_chi_, gf_dchi_, gf_ddchi_);

  const mat3<GF3D2<CCTK_REAL>, DN, DN> gf_gammat_(make_mat_gf);
  const mat3<vec3<GF3D2<CCTK_REAL>, DN>, DN, DN> gf_dgammat_(make_mat_vec_gf);
  const mat3<mat3<GF3D2<CCTK_REAL>, DN, DN>, DN, DN> gf_ddgammat_(
      make_mat_mat_gf);
  calc_derivs2(cctkGH, gf_gammat0_, gf_gammat_, gf_dgammat_, gf_ddgammat_);

  const GF3D2<CCTK_REAL> gf_Kh_(make_gf());
  const vec3<GF3D2<CCTK_REAL>, DN> gf_dKh_(make_vec_gf);
  calc_derivs(cctkGH, gf_Kh0_, gf_Kh_, gf_dKh_);

  const mat3<GF3D2<CCTK_REAL>, DN, DN> gf_At_(make_mat_gf);
  const mat3<vec3<GF3D2<CCTK_REAL>, DN>, DN, DN> gf_dAt_(make_mat_vec_gf);
  calc_derivs(cctkGH, gf_At0_, gf_At_, gf_dAt_);

  const vec3<GF3D2<CCTK_REAL>, UP> gf_Gamt_(make_vec_gf);
  const vec3<vec3<GF3D2<CCTK_REAL>, DN>, UP> gf_dGamt_(make_vec_vec_gf);
  calc_derivs(cctkGH, gf_Gamt0_, gf_Gamt_, gf_dGamt_);

  const GF3D2<CCTK_REAL> gf_Theta_(make_gf());
  const vec3<GF3D2<CCTK_REAL>, DN> gf_dTheta_(make_vec_gf);
  calc_derivs(cctkGH, gf_Theta0_, gf_Theta_, gf_dTheta_);

  //

  const GF3D2<CCTK_REAL> gf_ZtCx_(&layout, ZtCx);
  const GF3D2<CCTK_REAL> gf_ZtCy_(&layout, ZtCy);
  const GF3D2<CCTK_REAL> gf_ZtCz_(&layout, ZtCz);

  const GF3D2<CCTK_REAL> gf_HC_(&layout, HC);

  const GF3D2<CCTK_REAL> gf_MtCx_(&layout, MtCx);
  const GF3D2<CCTK_REAL> gf_MtCy_(&layout, MtCy);
  const GF3D2<CCTK_REAL> gf_MtCz_(&layout, MtCz);

  const GF3D2<CCTK_REAL> gf_allC_(&layout, allC);

  //

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    // Load and calculate

    const CCTK_REAL alphaG{1};
    const vec3<CCTK_REAL, DN> dalphaG{0, 0, 0};
    const mat3<CCTK_REAL, DN, DN> ddalphaG{0, 0, 0};

    const vec3<CCTK_REAL, UP> betaG{0, 0, 0};
    const vec3<vec3<CCTK_REAL, DN>, UP> dbetaG{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    const vec3<mat3<CCTK_REAL, DN, DN>, UP> ddbetaG{
        {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

    const CCTK_REAL rho{0};
    const vec3<CCTK_REAL, DN> Si{0, 0, 0};
    const mat3<CCTK_REAL, DN, DN> Sij{0, 0, 0};

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
