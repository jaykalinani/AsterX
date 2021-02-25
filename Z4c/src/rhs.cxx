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
#include <memory>

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
  vect<int, dim> imin, imax;
  GridDescBase(cctkGH).box_int<0, 0, 0>(nghostzones, imin, imax);
  // Suffix 1: with ghost zones, suffix 0: without ghost zones
  const GF3D2layout layout1(cctkGH, indextype);
  const GF3D2layout layout0(imin, imax);

  const GF3D2<const CCTK_REAL> gf_chi1(layout1, chi);

  const mat3<GF3D2<const CCTK_REAL>, DN, DN> gf_gammat1(
      GF3D2<const CCTK_REAL>(layout1, gammatxx),
      GF3D2<const CCTK_REAL>(layout1, gammatxy),
      GF3D2<const CCTK_REAL>(layout1, gammatxz),
      GF3D2<const CCTK_REAL>(layout1, gammatyy),
      GF3D2<const CCTK_REAL>(layout1, gammatyz),
      GF3D2<const CCTK_REAL>(layout1, gammatzz));

  const GF3D2<const CCTK_REAL> gf_Kh1(layout1, Kh);

  const mat3<GF3D2<const CCTK_REAL>, DN, DN> gf_At1(
      GF3D2<const CCTK_REAL>(layout1, Atxx),
      GF3D2<const CCTK_REAL>(layout1, Atxy),
      GF3D2<const CCTK_REAL>(layout1, Atxz),
      GF3D2<const CCTK_REAL>(layout1, Atyy),
      GF3D2<const CCTK_REAL>(layout1, Atyz),
      GF3D2<const CCTK_REAL>(layout1, Atzz));

  const vec3<GF3D2<const CCTK_REAL>, UP> gf_Gamt1(
      GF3D2<const CCTK_REAL>(layout1, Gamtx),
      GF3D2<const CCTK_REAL>(layout1, Gamty),
      GF3D2<const CCTK_REAL>(layout1, Gamtz));

  const GF3D2<const CCTK_REAL> gf_Theta1(layout1, Theta);

  const GF3D2<const CCTK_REAL> gf_alphaG1(layout1, alphaG);

  const vec3<GF3D2<const CCTK_REAL>, UP> gf_betaG1(
      GF3D2<const CCTK_REAL>(layout1, betaGx),
      GF3D2<const CCTK_REAL>(layout1, betaGy),
      GF3D2<const CCTK_REAL>(layout1, betaGz));

  //

  // Ideas:
  //
  // - Outline certain functions, e.g. `det` or `raise_index`. Ensure
  //   they are called with floating-point arguments, not tensor
  //   indices.

#if 0

  static mempool_set_t mempools;
  mempool_t &restrict mempool = mempools.get_mempool();

  const auto make_gf = [&]() { return GF3D2<CCTK_REAL>(layout0, mempool); };
  const auto make_vec_gf = [&](int) { return make_gf(); };
  const auto make_mat_gf = [&](int, int) { return make_gf(); };
  const auto make_vec_vec_gf = [&](int) { return make_vec_gf; };
  const auto make_vec_mat_gf = [&](int) { return make_mat_gf; };
  const auto make_mat_vec_gf = [&](int, int) { return make_vec_gf; };
  const auto make_mat_mat_gf = [&](int, int) { return make_mat_gf; };

  const GF3D2<CCTK_REAL> gf_chi0(make_gf());
  const vec3<GF3D2<CCTK_REAL>, DN> gf_dchi0(make_vec_gf);
  const mat3<GF3D2<CCTK_REAL>, DN, DN> gf_ddchi0(make_mat_gf);
  calc_derivs2(cctkGH, gf_chi1, gf_chi0, gf_dchi0, gf_ddchi0);

  const mat3<GF3D2<CCTK_REAL>, DN, DN> gf_gammat0(make_mat_gf);
  const mat3<vec3<GF3D2<CCTK_REAL>, DN>, DN, DN> gf_dgammat0(make_mat_vec_gf);
  const mat3<mat3<GF3D2<CCTK_REAL>, DN, DN>, DN, DN> gf_ddgammat0(
      make_mat_mat_gf);
  calc_derivs2(cctkGH, gf_gammat1, gf_gammat0, gf_dgammat0, gf_ddgammat0);

  const GF3D2<CCTK_REAL> gf_Kh0(make_gf());
  const vec3<GF3D2<CCTK_REAL>, DN> gf_dKh0(make_vec_gf);
  calc_derivs(cctkGH, gf_Kh1, gf_Kh0, gf_dKh0);

  const mat3<GF3D2<CCTK_REAL>, DN, DN> gf_At0(make_mat_gf);
  const mat3<vec3<GF3D2<CCTK_REAL>, DN>, DN, DN> gf_dAt0(make_mat_vec_gf);
  calc_derivs(cctkGH, gf_At1, gf_At0, gf_dAt0);

  const vec3<GF3D2<CCTK_REAL>, UP> gf_Gamt0(make_vec_gf);
  const vec3<vec3<GF3D2<CCTK_REAL>, DN>, UP> gf_dGamt0(make_vec_vec_gf);
  calc_derivs(cctkGH, gf_Gamt1, gf_Gamt0, gf_dGamt0);

  const GF3D2<CCTK_REAL> gf_Theta0(make_gf());
  const vec3<GF3D2<CCTK_REAL>, DN> gf_dTheta0(make_vec_gf);
  calc_derivs(cctkGH, gf_Theta1, gf_Theta0, gf_dTheta0);

  const GF3D2<CCTK_REAL> gf_alphaG0(make_gf());
  const vec3<GF3D2<CCTK_REAL>, DN> gf_dalphaG0(make_vec_gf);
  const mat3<GF3D2<CCTK_REAL>, DN, DN> gf_ddalphaG0(make_mat_gf);
  calc_derivs2(cctkGH, gf_alphaG1, gf_alphaG0, gf_dalphaG0, gf_ddalphaG0);

  const vec3<GF3D2<CCTK_REAL>, UP> gf_betaG0(make_vec_gf);
  const vec3<vec3<GF3D2<CCTK_REAL>, DN>, UP> gf_dbetaG0(make_vec_vec_gf);
  const vec3<mat3<GF3D2<CCTK_REAL>, DN, DN>, UP> gf_ddbetaG0(make_vec_mat_gf);
  calc_derivs2(cctkGH, gf_betaG1, gf_betaG0, gf_dbetaG0, gf_ddbetaG0);

#endif

#if 1

  static mempool_set_t mempools;
  mempool_t &restrict mempool = mempools.get_mempool();

  const auto make_gf = [&]() { return GF3D5<CCTK_REAL>(layout0, mempool); };
  const auto make_vec_gf = [&](int) { return make_gf(); };
  const auto make_mat_gf = [&](int, int) { return make_gf(); };
  const auto make_vec_vec_gf = [&](int) { return make_vec_gf; };
  const auto make_vec_mat_gf = [&](int) { return make_mat_gf; };
  const auto make_mat_vec_gf = [&](int, int) { return make_vec_gf; };
  const auto make_mat_mat_gf = [&](int, int) { return make_mat_gf; };

  const GF3D5<CCTK_REAL> gf_chi0(make_gf());
  const vec3<GF3D5<CCTK_REAL>, DN> gf_dchi0(make_vec_gf);
  const mat3<GF3D5<CCTK_REAL>, DN, DN> gf_ddchi0(make_mat_gf);
  calc_derivs2(cctkGH, gf_chi1, gf_chi0, gf_dchi0, gf_ddchi0, layout0);

  const mat3<GF3D5<CCTK_REAL>, DN, DN> gf_gammat0(make_mat_gf);
  const mat3<vec3<GF3D5<CCTK_REAL>, DN>, DN, DN> gf_dgammat0(make_mat_vec_gf);
  const mat3<mat3<GF3D5<CCTK_REAL>, DN, DN>, DN, DN> gf_ddgammat0(
      make_mat_mat_gf);
  calc_derivs2(cctkGH, gf_gammat1, gf_gammat0, gf_dgammat0, gf_ddgammat0,
               layout0);

  const GF3D5<CCTK_REAL> gf_Kh0(make_gf());
  const vec3<GF3D5<CCTK_REAL>, DN> gf_dKh0(make_vec_gf);
  calc_derivs(cctkGH, gf_Kh1, gf_Kh0, gf_dKh0, layout0);

  const mat3<GF3D5<CCTK_REAL>, DN, DN> gf_At0(make_mat_gf);
  const mat3<vec3<GF3D5<CCTK_REAL>, DN>, DN, DN> gf_dAt0(make_mat_vec_gf);
  calc_derivs(cctkGH, gf_At1, gf_At0, gf_dAt0, layout0);

  const vec3<GF3D5<CCTK_REAL>, UP> gf_Gamt0(make_vec_gf);
  const vec3<vec3<GF3D5<CCTK_REAL>, DN>, UP> gf_dGamt0(make_vec_vec_gf);
  calc_derivs(cctkGH, gf_Gamt1, gf_Gamt0, gf_dGamt0, layout0);

  const GF3D5<CCTK_REAL> gf_Theta0(make_gf());
  const vec3<GF3D5<CCTK_REAL>, DN> gf_dTheta0(make_vec_gf);
  calc_derivs(cctkGH, gf_Theta1, gf_Theta0, gf_dTheta0, layout0);

  const GF3D5<CCTK_REAL> gf_alphaG0(make_gf());
  const vec3<GF3D5<CCTK_REAL>, DN> gf_dalphaG0(make_vec_gf);
  const mat3<GF3D5<CCTK_REAL>, DN, DN> gf_ddalphaG0(make_mat_gf);
  calc_derivs2(cctkGH, gf_alphaG1, gf_alphaG0, gf_dalphaG0, gf_ddalphaG0,
               layout0);

  const vec3<GF3D5<CCTK_REAL>, UP> gf_betaG0(make_vec_gf);
  const vec3<vec3<GF3D5<CCTK_REAL>, DN>, UP> gf_dbetaG0(make_vec_vec_gf);
  const vec3<mat3<GF3D5<CCTK_REAL>, DN, DN>, UP> gf_ddbetaG0(make_vec_mat_gf);
  calc_derivs2(cctkGH, gf_betaG1, gf_betaG0, gf_dbetaG0, gf_ddbetaG0, layout0);

#endif

  //

  const GF3D2<const CCTK_REAL> gf_eTtt1(layout1, eTtt);

  const vec3<GF3D2<const CCTK_REAL>, DN> gf_eTti1(
      GF3D2<const CCTK_REAL>(layout1, eTtx),
      GF3D2<const CCTK_REAL>(layout1, eTty),
      GF3D2<const CCTK_REAL>(layout1, eTtz));

  const mat3<GF3D2<const CCTK_REAL>, DN, DN> gf_eTij1(
      GF3D2<const CCTK_REAL>(layout1, eTxx),
      GF3D2<const CCTK_REAL>(layout1, eTxy),
      GF3D2<const CCTK_REAL>(layout1, eTxz),
      GF3D2<const CCTK_REAL>(layout1, eTyy),
      GF3D2<const CCTK_REAL>(layout1, eTyz),
      GF3D2<const CCTK_REAL>(layout1, eTzz));

  //

  const GF3D2<CCTK_REAL> gf_chi_rhs1(layout1, chi_rhs);

  const GF3D2<CCTK_REAL> gf_gammatxx_rhs1(layout1, gammatxx_rhs);
  const GF3D2<CCTK_REAL> gf_gammatxy_rhs1(layout1, gammatxy_rhs);
  const GF3D2<CCTK_REAL> gf_gammatxz_rhs1(layout1, gammatxz_rhs);
  const GF3D2<CCTK_REAL> gf_gammatyy_rhs1(layout1, gammatyy_rhs);
  const GF3D2<CCTK_REAL> gf_gammatyz_rhs1(layout1, gammatyz_rhs);
  const GF3D2<CCTK_REAL> gf_gammatzz_rhs1(layout1, gammatzz_rhs);

  const GF3D2<CCTK_REAL> gf_Kh_rhs1(layout1, Kh_rhs);

  const GF3D2<CCTK_REAL> gf_Atxx_rhs1(layout1, Atxx_rhs);
  const GF3D2<CCTK_REAL> gf_Atxy_rhs1(layout1, Atxy_rhs);
  const GF3D2<CCTK_REAL> gf_Atxz_rhs1(layout1, Atxz_rhs);
  const GF3D2<CCTK_REAL> gf_Atyy_rhs1(layout1, Atyy_rhs);
  const GF3D2<CCTK_REAL> gf_Atyz_rhs1(layout1, Atyz_rhs);
  const GF3D2<CCTK_REAL> gf_Atzz_rhs1(layout1, Atzz_rhs);

  const GF3D2<CCTK_REAL> gf_Gamtx_rhs1(layout1, Gamtx_rhs);
  const GF3D2<CCTK_REAL> gf_Gamty_rhs1(layout1, Gamty_rhs);
  const GF3D2<CCTK_REAL> gf_Gamtz_rhs1(layout1, Gamtz_rhs);

  const GF3D2<CCTK_REAL> gf_Theta_rhs1(layout1, Theta_rhs);

  const GF3D2<CCTK_REAL> gf_alphaG_rhs1(layout1, alphaG_rhs);

  const GF3D2<CCTK_REAL> gf_betaGx_rhs1(layout1, betaGx_rhs);
  const GF3D2<CCTK_REAL> gf_betaGy_rhs1(layout1, betaGy_rhs);
  const GF3D2<CCTK_REAL> gf_betaGz_rhs1(layout1, betaGz_rhs);

  //

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) Z4C_INLINE {
  // Load and calculate
#if 0
    const z4c_vars<CCTK_REAL> vars(kappa1, kappa2, f_mu_L, f_mu_S, eta, //
                                   gf_chi0(p.I), gf_dchi0(p.I),
                                   gf_ddchi0(p.I), //
                                   gf_gammat0(p.I), gf_dgammat0(p.I),
                                   gf_ddgammat0(p.I),               //
                                   gf_Kh0(p.I), gf_dKh0(p.I),       //
                                   gf_At0(p.I), gf_dAt0(p.I),       //
                                   gf_Gamt0(p.I), gf_dGamt0(p.I),   //
                                   gf_Theta0(p.I), gf_dTheta0(p.I), //
                                   gf_alphaG0(p.I), gf_dalphaG0(p.I),
                                   gf_ddalphaG0(p.I), //
                                   gf_betaG0(p.I), gf_dbetaG0(p.I),
                                   gf_ddbetaG0(p.I), //
                                   gf_eTtt1(p.I), gf_eTti1(p.I), gf_eTij1(p.I));
#endif
#if 1
    const z4c_vars<CCTK_REAL> vars(
        kappa1, kappa2, f_mu_L, f_mu_S, eta, //
        gf_chi0(layout0, p.I), gf_dchi0(layout0, p.I),
        gf_ddchi0(layout0, p.I), //
        gf_gammat0(layout0, p.I), gf_dgammat0(layout0, p.I),
        gf_ddgammat0(layout0, p.I),                        //
        gf_Kh0(layout0, p.I), gf_dKh0(layout0, p.I),       //
        gf_At0(layout0, p.I), gf_dAt0(layout0, p.I),       //
        gf_Gamt0(layout0, p.I), gf_dGamt0(layout0, p.I),   //
        gf_Theta0(layout0, p.I), gf_dTheta0(layout0, p.I), //
        gf_alphaG0(layout0, p.I), gf_dalphaG0(layout0, p.I),
        gf_ddalphaG0(layout0, p.I), //
        gf_betaG0(layout0, p.I), gf_dbetaG0(layout0, p.I),
        gf_ddbetaG0(layout0, p.I), //
        gf_eTtt1(p.I), gf_eTti1(p.I), gf_eTij1(p.I));
#endif

    // Store
    gf_chi_rhs1(p.I) = vars.chi_rhs;
    vars.gammat_rhs.store(gf_gammatxx_rhs1, gf_gammatxy_rhs1, gf_gammatxz_rhs1,
                          gf_gammatyy_rhs1, gf_gammatyz_rhs1, gf_gammatzz_rhs1,
                          p.I);
    gf_Kh_rhs1(p.I) = vars.Kh_rhs;
    vars.At_rhs.store(gf_Atxx_rhs1, gf_Atxy_rhs1, gf_Atxz_rhs1, gf_Atyy_rhs1,
                      gf_Atyz_rhs1, gf_Atzz_rhs1, p.I);
    vars.Gamt_rhs.store(gf_Gamtx_rhs1, gf_Gamty_rhs1, gf_Gamtz_rhs1, p.I);
    gf_Theta_rhs1(p.I) = vars.Theta_rhs;
    gf_alphaG_rhs1(p.I) = vars.alphaG_rhs;
    vars.betaG_rhs.store(gf_betaGx_rhs1, gf_betaGy_rhs1, gf_betaGz_rhs1, p.I);
  });

  // Upwind and dissipation terms

  // TODO: Consider fusing the loops to reduce memory bandwidth

  apply_upwind_diss(cctkGH, gf_chi1, gf_betaG1, gf_chi_rhs1);

  apply_upwind_diss(cctkGH, gf_gammat1(0, 0), gf_betaG1, gf_gammatxx_rhs1);
  apply_upwind_diss(cctkGH, gf_gammat1(0, 1), gf_betaG1, gf_gammatxy_rhs1);
  apply_upwind_diss(cctkGH, gf_gammat1(0, 2), gf_betaG1, gf_gammatxz_rhs1);
  apply_upwind_diss(cctkGH, gf_gammat1(1, 1), gf_betaG1, gf_gammatyy_rhs1);
  apply_upwind_diss(cctkGH, gf_gammat1(1, 2), gf_betaG1, gf_gammatyz_rhs1);
  apply_upwind_diss(cctkGH, gf_gammat1(2, 2), gf_betaG1, gf_gammatzz_rhs1);

  apply_upwind_diss(cctkGH, gf_Kh1, gf_betaG1, gf_Kh_rhs1);

  apply_upwind_diss(cctkGH, gf_At1(0, 0), gf_betaG1, gf_Atxx_rhs1);
  apply_upwind_diss(cctkGH, gf_At1(0, 1), gf_betaG1, gf_Atxy_rhs1);
  apply_upwind_diss(cctkGH, gf_At1(0, 2), gf_betaG1, gf_Atxz_rhs1);
  apply_upwind_diss(cctkGH, gf_At1(1, 1), gf_betaG1, gf_Atyy_rhs1);
  apply_upwind_diss(cctkGH, gf_At1(1, 2), gf_betaG1, gf_Atyz_rhs1);
  apply_upwind_diss(cctkGH, gf_At1(2, 2), gf_betaG1, gf_Atzz_rhs1);

  apply_upwind_diss(cctkGH, gf_Gamt1(0), gf_betaG1, gf_Gamtx_rhs1);
  apply_upwind_diss(cctkGH, gf_Gamt1(1), gf_betaG1, gf_Gamty_rhs1);
  apply_upwind_diss(cctkGH, gf_Gamt1(2), gf_betaG1, gf_Gamtz_rhs1);

  apply_upwind_diss(cctkGH, gf_Theta1, gf_betaG1, gf_Theta_rhs1);

  apply_upwind_diss(cctkGH, gf_alphaG1, gf_betaG1, gf_alphaG_rhs1);

  apply_upwind_diss(cctkGH, gf_betaG1(0), gf_betaG1, gf_betaGx_rhs1);
  apply_upwind_diss(cctkGH, gf_betaG1(1), gf_betaG1, gf_betaGy_rhs1);
  apply_upwind_diss(cctkGH, gf_betaG1(2), gf_betaG1, gf_betaGz_rhs1);
}

} // namespace Z4c
