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

  const GF3D1<const CCTK_REAL> gf_alphaG_(cctkGH, indextype, nghostzones,
                                          alphaG);

  const vec3<GF3D1<const CCTK_REAL>, UP> gf_betaG_(
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, betaGx),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, betaGy),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, betaGz));

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
  const auto make_vec_mat_gf = [&](int) {
    return [&](int, int) { return make_gf(); };
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

  const vec3<GF3D1<CCTK_REAL>, DN> gf_dalphaG_(make_vec_gf);
  const mat3<GF3D1<CCTK_REAL>, DN, DN> gf_ddalphaG_(make_mat_gf);
  calc_derivs2(cctkGH, gf_alphaG_, gf_dalphaG_, gf_ddalphaG_);

  const vec3<vec3<GF3D1<CCTK_REAL>, DN>, UP> gf_dbetaG_(make_vec_vec_gf);
  const vec3<mat3<GF3D1<CCTK_REAL>, DN, DN>, UP> gf_ddbetaG_(make_vec_mat_gf);
  calc_derivs2(cctkGH, gf_betaG_, gf_dbetaG_, gf_ddbetaG_);

  //

  const GF3D1<const CCTK_REAL> gf_eTtt_(cctkGH, indextype, nghostzones, eTtt);

  const vec3<GF3D1<const CCTK_REAL>, DN> gf_eTti_(
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, eTtx),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, eTty),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, eTtz));

  const mat3<GF3D1<const CCTK_REAL>, DN, DN> gf_eTij_(
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, eTxx),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, eTxy),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, eTxz),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, eTyy),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, eTyz),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, eTzz));

  //

  const GF3D1<CCTK_REAL> gf_chi_rhs_(cctkGH, indextype, noghosts, chi_rhs);

  const GF3D1<CCTK_REAL> gf_gammatxx_rhs_(cctkGH, indextype, noghosts,
                                          gammatxx_rhs);
  const GF3D1<CCTK_REAL> gf_gammatxy_rhs_(cctkGH, indextype, noghosts,
                                          gammatxy_rhs);
  const GF3D1<CCTK_REAL> gf_gammatxz_rhs_(cctkGH, indextype, noghosts,
                                          gammatxz_rhs);
  const GF3D1<CCTK_REAL> gf_gammatyy_rhs_(cctkGH, indextype, noghosts,
                                          gammatyy_rhs);
  const GF3D1<CCTK_REAL> gf_gammatyz_rhs_(cctkGH, indextype, noghosts,
                                          gammatyz_rhs);
  const GF3D1<CCTK_REAL> gf_gammatzz_rhs_(cctkGH, indextype, noghosts,
                                          gammatzz_rhs);

  const GF3D1<CCTK_REAL> gf_Kh_rhs_(cctkGH, indextype, noghosts, Kh_rhs);

  const GF3D1<CCTK_REAL> gf_Atxx_rhs_(cctkGH, indextype, noghosts, Atxx_rhs);
  const GF3D1<CCTK_REAL> gf_Atxy_rhs_(cctkGH, indextype, noghosts, Atxy_rhs);
  const GF3D1<CCTK_REAL> gf_Atxz_rhs_(cctkGH, indextype, noghosts, Atxz_rhs);
  const GF3D1<CCTK_REAL> gf_Atyy_rhs_(cctkGH, indextype, noghosts, Atyy_rhs);
  const GF3D1<CCTK_REAL> gf_Atyz_rhs_(cctkGH, indextype, noghosts, Atyz_rhs);
  const GF3D1<CCTK_REAL> gf_Atzz_rhs_(cctkGH, indextype, noghosts, Atzz_rhs);

  const GF3D1<CCTK_REAL> gf_Gamtx_rhs_(cctkGH, indextype, noghosts, Gamtx_rhs);
  const GF3D1<CCTK_REAL> gf_Gamty_rhs_(cctkGH, indextype, noghosts, Gamty_rhs);
  const GF3D1<CCTK_REAL> gf_Gamtz_rhs_(cctkGH, indextype, noghosts, Gamtz_rhs);

  const GF3D1<CCTK_REAL> gf_Theta_rhs_(cctkGH, indextype, noghosts, Theta_rhs);

  const GF3D1<CCTK_REAL> gf_alphaG_rhs_(cctkGH, indextype, noghosts,
                                        alphaG_rhs);

  const GF3D1<CCTK_REAL> gf_betaGx_rhs_(cctkGH, indextype, noghosts,
                                        betaGx_rhs);
  const GF3D1<CCTK_REAL> gf_betaGy_rhs_(cctkGH, indextype, noghosts,
                                        betaGy_rhs);
  const GF3D1<CCTK_REAL> gf_betaGz_rhs_(cctkGH, indextype, noghosts,
                                        betaGz_rhs);

  //

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) Z4C_INLINE {
    // Load and calculate

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
        gf_eTtt_(p.I), gf_eTti_(p.I), gf_eTij_(p.I));

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
