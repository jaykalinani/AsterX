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
  const array<int, dim> nghostzones = {cctk_nghostzones[0], cctk_nghostzones[1],
                                       cctk_nghostzones[2]};
  const array<int, dim> noghosts = {0, 0, 0};

  const GF3D1<const CCTK_REAL> gf_chi_(cctkGH, indextype, nghostzones, chi);

  const mat3<GF3D1<const CCTK_REAL>, DN, DN> gf_gammat_(
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, gammatxx),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, gammatxy),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, gammatxz),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, gammatyy),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, gammatyz),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, gammatzz));

  const GF3D1<const CCTK_REAL> gf_Kh_(cctkGH, indextype, nghostzones, Kh);

  const mat3<GF3D1<const CCTK_REAL>, DN, DN> gf_At_(
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, Atxx),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, Atxy),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, Atxz),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, Atyy),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, Atyz),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, Atzz));

  const vec3<GF3D1<const CCTK_REAL>, UP> gf_Gamt_(
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, Gamtx),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, Gamty),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, Gamtz));

  const GF3D1<const CCTK_REAL> gf_Theta_(cctkGH, indextype, nghostzones, Theta);

  //

  static mempool_set_t mempools;
  mempool_t &restrict mempool = mempools.get_mempool();

  const auto make_gf = [&]() {
    return GF3D1<CCTK_REAL>(cctkGH, indextype, noghosts, mempool);
  };
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

  const vec3<GF3D1<CCTK_REAL>, DN> gf_dchi_(make_vec_gf);
  const mat3<GF3D1<CCTK_REAL>, DN, DN> gf_ddchi_(make_mat_gf);
  calc_derivs2(cctkGH, gf_chi_, gf_dchi_, gf_ddchi_);

  const mat3<vec3<GF3D1<CCTK_REAL>, DN>, DN, DN> gf_dgammat_(make_mat_vec_gf);
  const mat3<mat3<GF3D1<CCTK_REAL>, DN, DN>, DN, DN> gf_ddgammat_(
      make_mat_mat_gf);
  calc_derivs2(cctkGH, gf_gammat_, gf_dgammat_, gf_ddgammat_);

  const vec3<GF3D1<CCTK_REAL>, DN> gf_dKh_(make_vec_gf);
  calc_derivs(cctkGH, gf_Kh_, gf_dKh_);

  const mat3<vec3<GF3D1<CCTK_REAL>, DN>, DN, DN> gf_dAt_(make_mat_vec_gf);
  calc_derivs(cctkGH, gf_At_, gf_dAt_);

  const vec3<vec3<GF3D1<CCTK_REAL>, DN>, UP> gf_dGamt_(make_vec_vec_gf);
  calc_derivs(cctkGH, gf_Gamt_, gf_dGamt_);

  const vec3<GF3D1<CCTK_REAL>, DN> gf_dTheta_(make_vec_gf);
  calc_derivs(cctkGH, gf_Theta_, gf_dTheta_);

  //

  const GF3D1<CCTK_REAL> gf_ZtCx_(cctkGH, indextype, nghostzones, ZtCx);
  const GF3D1<CCTK_REAL> gf_ZtCy_(cctkGH, indextype, nghostzones, ZtCy);
  const GF3D1<CCTK_REAL> gf_ZtCz_(cctkGH, indextype, nghostzones, ZtCz);

  const GF3D1<CCTK_REAL> gf_HC_(cctkGH, indextype, nghostzones, HC);

  const GF3D1<CCTK_REAL> gf_MtCx_(cctkGH, indextype, nghostzones, MtCx);
  const GF3D1<CCTK_REAL> gf_MtCy_(cctkGH, indextype, nghostzones, MtCy);
  const GF3D1<CCTK_REAL> gf_MtCz_(cctkGH, indextype, nghostzones, MtCz);

  const GF3D1<CCTK_REAL> gf_allC_(cctkGH, indextype, nghostzones, allC);

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
