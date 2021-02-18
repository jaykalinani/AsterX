#include "derivs.hxx"
#include "physics.hxx"
#include "tensor.hxx"
#include "weyl_vars.hxx"

#include <loop.hxx>
#include <mempool.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace Weyl {
using namespace Loop;
using namespace std;

extern "C" void Weyl_Weyl(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Weyl_Weyl;
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
  const array<int, dim> noghosts = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);
  const GF3D2layout layout0(cctkGH, indextype, noghosts);

  const mat3<GF3D2<const CCTK_REAL>, DN, DN> gf_gamma0_(
      GF3D2<const CCTK_REAL>(&layout, gxx),
      GF3D2<const CCTK_REAL>(&layout, gxy),
      GF3D2<const CCTK_REAL>(&layout, gxz),
      GF3D2<const CCTK_REAL>(&layout, gyy),
      GF3D2<const CCTK_REAL>(&layout, gyz),
      GF3D2<const CCTK_REAL>(&layout, gzz));

  const GF3D2<const CCTK_REAL> gf_alpha0_(&layout, alp);

  const vec3<GF3D2<const CCTK_REAL>, UP> gf_beta0_(
      GF3D2<const CCTK_REAL>(&layout, betax),
      GF3D2<const CCTK_REAL>(&layout, betay),
      GF3D2<const CCTK_REAL>(&layout, betaz));

  const mat3<GF3D2<const CCTK_REAL>, DN, DN> gf_k0_(
      GF3D2<const CCTK_REAL>(&layout, kxx),
      GF3D2<const CCTK_REAL>(&layout, kxy),
      GF3D2<const CCTK_REAL>(&layout, kxz),
      GF3D2<const CCTK_REAL>(&layout, kyy),
      GF3D2<const CCTK_REAL>(&layout, kyz),
      GF3D2<const CCTK_REAL>(&layout, kzz));

  const GF3D2<const CCTK_REAL> gf_dtalpha0_(&layout, dtalp);

  const vec3<GF3D2<const CCTK_REAL>, UP> gf_dtbeta0_(
      GF3D2<const CCTK_REAL>(&layout, dtbetax),
      GF3D2<const CCTK_REAL>(&layout, dtbetay),
      GF3D2<const CCTK_REAL>(&layout, dtbetaz));

  const mat3<GF3D2<const CCTK_REAL>, DN, DN> gf_dtk_(
      GF3D2<const CCTK_REAL>(&layout, dtkxx),
      GF3D2<const CCTK_REAL>(&layout, dtkxy),
      GF3D2<const CCTK_REAL>(&layout, dtkxz),
      GF3D2<const CCTK_REAL>(&layout, dtkyy),
      GF3D2<const CCTK_REAL>(&layout, dtkyz),
      GF3D2<const CCTK_REAL>(&layout, dtkzz));

  const GF3D2<const CCTK_REAL> gf_dt2alpha_(&layout, dt2alp);

  const vec3<GF3D2<const CCTK_REAL>, UP> gf_dt2beta_(
      GF3D2<const CCTK_REAL>(&layout, dt2betax),
      GF3D2<const CCTK_REAL>(&layout, dt2betay),
      GF3D2<const CCTK_REAL>(&layout, dt2betaz));

  //

  static mempool_set_t mempools;
  mempool_t &restrict mempool = mempools.get_mempool();

  const auto make_gf = [&]() { return GF3D2<CCTK_REAL>(&layout0, mempool); };
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

  const mat3<GF3D2<CCTK_REAL>, DN, DN> gf_gamma_(make_mat_gf);
  const mat3<vec3<GF3D2<CCTK_REAL>, DN>, DN, DN> gf_dgamma_(make_mat_vec_gf);
  const mat3<mat3<GF3D2<CCTK_REAL>, DN, DN>, DN, DN> gf_ddgamma_(
      make_mat_mat_gf);
  calc_derivs2(cctkGH, gf_gamma0_, gf_gamma_, gf_dgamma_, gf_ddgamma_);

  const GF3D2<CCTK_REAL> gf_alpha_(make_gf());
  const vec3<GF3D2<CCTK_REAL>, DN> gf_dalpha_(make_vec_gf);
  const mat3<GF3D2<CCTK_REAL>, DN, DN> gf_ddalpha_(make_mat_gf);
  calc_derivs2(cctkGH, gf_alpha0_, gf_alpha_, gf_dalpha_, gf_ddalpha_);

  const vec3<GF3D2<CCTK_REAL>, UP> gf_beta_(make_vec_gf);
  const vec3<vec3<GF3D2<CCTK_REAL>, DN>, UP> gf_dbeta_(make_vec_vec_gf);
  const vec3<mat3<GF3D2<CCTK_REAL>, DN, DN>, UP> gf_ddbeta_(make_vec_mat_gf);
  calc_derivs2(cctkGH, gf_beta0_, gf_beta_, gf_dbeta_, gf_ddbeta_);

  const mat3<GF3D2<CCTK_REAL>, DN, DN> gf_k_(make_mat_gf);
  const mat3<vec3<GF3D2<CCTK_REAL>, DN>, DN, DN> gf_dk_(make_mat_vec_gf);
  calc_derivs(cctkGH, gf_k0_, gf_k_, gf_dk_);

  const GF3D2<CCTK_REAL> gf_dtalpha_(make_gf());
  const vec3<GF3D2<CCTK_REAL>, DN> gf_ddtalpha_(make_vec_gf);
  calc_derivs(cctkGH, gf_dtalpha0_, gf_dtalpha_, gf_ddtalpha_);

  const vec3<GF3D2<CCTK_REAL>, UP> gf_dtbeta_(make_vec_gf);
  const vec3<vec3<GF3D2<CCTK_REAL>, DN>, UP> gf_ddtbeta_(make_vec_vec_gf);
  calc_derivs(cctkGH, gf_dtbeta0_, gf_dtbeta_, gf_ddtbeta_);

  //

  const GF3D2<CCTK_REAL> gf_g4tt_(&layout, g4tt);
  const GF3D2<CCTK_REAL> gf_g4tx_(&layout, g4tx);
  const GF3D2<CCTK_REAL> gf_g4ty_(&layout, g4ty);
  const GF3D2<CCTK_REAL> gf_g4tz_(&layout, g4tz);
  const GF3D2<CCTK_REAL> gf_g4xx_(&layout, g4xx);
  const GF3D2<CCTK_REAL> gf_g4xy_(&layout, g4xy);
  const GF3D2<CCTK_REAL> gf_g4xz_(&layout, g4xz);
  const GF3D2<CCTK_REAL> gf_g4yy_(&layout, g4yy);
  const GF3D2<CCTK_REAL> gf_g4yz_(&layout, g4yz);
  const GF3D2<CCTK_REAL> gf_g4zz_(&layout, g4zz);

  // const GF3D2<CCTK_REAL> gf_Gamma4ttt_(&layout, Gamma4ttt);
  // const GF3D2<CCTK_REAL> gf_Gamma4ttx_(&layout, Gamma4ttx);
  // const GF3D2<CCTK_REAL> gf_Gamma4tty_(&layout, Gamma4tty);
  // const GF3D2<CCTK_REAL> gf_Gamma4ttz_(&layout, Gamma4ttz);
  // const GF3D2<CCTK_REAL> gf_Gamma4txx_(&layout, Gamma4txx);
  // const GF3D2<CCTK_REAL> gf_Gamma4txy_(&layout, Gamma4txy);
  // const GF3D2<CCTK_REAL> gf_Gamma4txz_(&layout, Gamma4txz);
  // const GF3D2<CCTK_REAL> gf_Gamma4tyy_(&layout, Gamma4tyy);
  // const GF3D2<CCTK_REAL> gf_Gamma4tyz_(&layout, Gamma4tyz);
  // const GF3D2<CCTK_REAL> gf_Gamma4tzz_(&layout, Gamma4tzz);

  // const GF3D2<CCTK_REAL> gf_Gamma4xtt_(&layout, Gamma4xtt);
  // const GF3D2<CCTK_REAL> gf_Gamma4xtx_(&layout, Gamma4xtx);
  // const GF3D2<CCTK_REAL> gf_Gamma4xty_(&layout, Gamma4xty);
  // const GF3D2<CCTK_REAL> gf_Gamma4xtz_(&layout, Gamma4xtz);
  // const GF3D2<CCTK_REAL> gf_Gamma4xxx_(&layout, Gamma4xxx);
  // const GF3D2<CCTK_REAL> gf_Gamma4xxy_(&layout, Gamma4xxy);
  // const GF3D2<CCTK_REAL> gf_Gamma4xxz_(&layout, Gamma4xxz);
  // const GF3D2<CCTK_REAL> gf_Gamma4xyy_(&layout, Gamma4xyy);
  // const GF3D2<CCTK_REAL> gf_Gamma4xyz_(&layout, Gamma4xyz);
  // const GF3D2<CCTK_REAL> gf_Gamma4xzz_(&layout, Gamma4xzz);

  // const GF3D2<CCTK_REAL> gf_Gamma4ytt_(&layout, Gamma4ytt);
  // const GF3D2<CCTK_REAL> gf_Gamma4ytx_(&layout, Gamma4ytx);
  // const GF3D2<CCTK_REAL> gf_Gamma4yty_(&layout, Gamma4yty);
  // const GF3D2<CCTK_REAL> gf_Gamma4ytz_(&layout, Gamma4ytz);
  // const GF3D2<CCTK_REAL> gf_Gamma4yxx_(&layout, Gamma4yxx);
  // const GF3D2<CCTK_REAL> gf_Gamma4yxy_(&layout, Gamma4yxy);
  // const GF3D2<CCTK_REAL> gf_Gamma4yxz_(&layout, Gamma4yxz);
  // const GF3D2<CCTK_REAL> gf_Gamma4yyy_(&layout, Gamma4yyy);
  // const GF3D2<CCTK_REAL> gf_Gamma4yyz_(&layout, Gamma4yyz);
  // const GF3D2<CCTK_REAL> gf_Gamma4yzz_(&layout, Gamma4yzz);

  // const GF3D2<CCTK_REAL> gf_Gamma4ztt_(&layout, Gamma4ztt);
  // const GF3D2<CCTK_REAL> gf_Gamma4ztx_(&layout, Gamma4ztx);
  // const GF3D2<CCTK_REAL> gf_Gamma4zty_(&layout, Gamma4zty);
  // const GF3D2<CCTK_REAL> gf_Gamma4ztz_(&layout, Gamma4ztz);
  // const GF3D2<CCTK_REAL> gf_Gamma4zxx_(&layout, Gamma4zxx);
  // const GF3D2<CCTK_REAL> gf_Gamma4zxy_(&layout, Gamma4zxy);
  // const GF3D2<CCTK_REAL> gf_Gamma4zxz_(&layout, Gamma4zxz);
  // const GF3D2<CCTK_REAL> gf_Gamma4zyy_(&layout, Gamma4zyy);
  // const GF3D2<CCTK_REAL> gf_Gamma4zyz_(&layout, Gamma4zyz);
  // const GF3D2<CCTK_REAL> gf_Gamma4zzz_(&layout, Gamma4zzz);

  // const GF3D2<CCTK_REAL> gf_rm4txtx_(&layout, rm4txtx);
  // const GF3D2<CCTK_REAL> gf_rm4txty_(&layout, rm4txty);
  // const GF3D2<CCTK_REAL> gf_rm4txtz_(&layout, rm4txtz);
  // const GF3D2<CCTK_REAL> gf_rm4txxy_(&layout, rm4txxy);
  // const GF3D2<CCTK_REAL> gf_rm4txxz_(&layout, rm4txxz);
  // const GF3D2<CCTK_REAL> gf_rm4txyz_(&layout, rm4txyz);

  // const GF3D2<CCTK_REAL> gf_rm4tyty_(&layout, rm4tyty);
  // const GF3D2<CCTK_REAL> gf_rm4tytz_(&layout, rm4tytz);
  // const GF3D2<CCTK_REAL> gf_rm4tyxy_(&layout, rm4tyxy);
  // const GF3D2<CCTK_REAL> gf_rm4tyxz_(&layout, rm4tyxz);
  // const GF3D2<CCTK_REAL> gf_rm4tyyz_(&layout, rm4tyyz);

  // const GF3D2<CCTK_REAL> gf_rm4tztz_(&layout, rm4tztz);
  // const GF3D2<CCTK_REAL> gf_rm4tzxy_(&layout, rm4tzxy);
  // const GF3D2<CCTK_REAL> gf_rm4tzxz_(&layout, rm4tzxz);
  // const GF3D2<CCTK_REAL> gf_rm4tzyz_(&layout, rm4tzyz);

  // const GF3D2<CCTK_REAL> gf_rm4xyxy_(&layout, rm4xyxy);
  // const GF3D2<CCTK_REAL> gf_rm4xyxz_(&layout, rm4xyxz);
  // const GF3D2<CCTK_REAL> gf_rm4xyyz_(&layout, rm4xyyz);

  // const GF3D2<CCTK_REAL> gf_rm4xzxz_(&layout, rm4xzxz);
  // const GF3D2<CCTK_REAL> gf_rm4xzyz_(&layout, rm4xzyz);

  // const GF3D2<CCTK_REAL> gf_rm4yzyz_(&layout, rm4yzyz);

  // const GF3D2<CCTK_REAL> gf_r4tt_(&layout, r4tt);
  // const GF3D2<CCTK_REAL> gf_r4tx_(&layout, r4tx);
  // const GF3D2<CCTK_REAL> gf_r4ty_(&layout, r4ty);
  // const GF3D2<CCTK_REAL> gf_r4tz_(&layout, r4tz);
  // const GF3D2<CCTK_REAL> gf_r4xx_(&layout, r4xx);
  // const GF3D2<CCTK_REAL> gf_r4xy_(&layout, r4xy);
  // const GF3D2<CCTK_REAL> gf_r4xz_(&layout, r4xz);
  // const GF3D2<CCTK_REAL> gf_r4yy_(&layout, r4yy);
  // const GF3D2<CCTK_REAL> gf_r4yz_(&layout, r4yz);
  // const GF3D2<CCTK_REAL> gf_r4zz_(&layout, r4zz);

  // const GF3D2<CCTK_REAL> gf_rsc4_(&layout, rsc4);

  // const GF3D2<CCTK_REAL> gf_c4txtx_(&layout, c4txtx);
  // const GF3D2<CCTK_REAL> gf_c4txty_(&layout, c4txty);
  // const GF3D2<CCTK_REAL> gf_c4txtz_(&layout, c4txtz);
  // const GF3D2<CCTK_REAL> gf_c4txxy_(&layout, c4txxy);
  // const GF3D2<CCTK_REAL> gf_c4txxz_(&layout, c4txxz);
  // const GF3D2<CCTK_REAL> gf_c4txyz_(&layout, c4txyz);

  // const GF3D2<CCTK_REAL> gf_c4tyty_(&layout, c4tyty);
  // const GF3D2<CCTK_REAL> gf_c4tytz_(&layout, c4tytz);
  // const GF3D2<CCTK_REAL> gf_c4tyxy_(&layout, c4tyxy);
  // const GF3D2<CCTK_REAL> gf_c4tyxz_(&layout, c4tyxz);
  // const GF3D2<CCTK_REAL> gf_c4tyyz_(&layout, c4tyyz);

  // const GF3D2<CCTK_REAL> gf_c4tztz_(&layout, c4tztz);
  // const GF3D2<CCTK_REAL> gf_c4tzxy_(&layout, c4tzxy);
  // const GF3D2<CCTK_REAL> gf_c4tzxz_(&layout, c4tzxz);
  // const GF3D2<CCTK_REAL> gf_c4tzyz_(&layout, c4tzyz);

  // const GF3D2<CCTK_REAL> gf_c4xyxy_(&layout, c4xyxy);
  // const GF3D2<CCTK_REAL> gf_c4xyxz_(&layout, c4xyxz);
  // const GF3D2<CCTK_REAL> gf_c4xyyz_(&layout, c4xyyz);

  // const GF3D2<CCTK_REAL> gf_c4xzxz_(&layout, c4xzxz);
  // const GF3D2<CCTK_REAL> gf_c4xzyz_(&layout, c4xzyz);

  // const GF3D2<CCTK_REAL> gf_c4yzyz_(&layout, c4yzyz);

  // //

  // const GF3D2<CCTK_REAL> gf_lt_(&layout, lt);
  // const GF3D2<CCTK_REAL> gf_lx_(&layout, lx);
  // const GF3D2<CCTK_REAL> gf_ly_(&layout, ly);
  // const GF3D2<CCTK_REAL> gf_lz_(&layout, lz);

  // const GF3D2<CCTK_REAL> gf_nt_(&layout, nt);
  // const GF3D2<CCTK_REAL> gf_nx_(&layout, nx);
  // const GF3D2<CCTK_REAL> gf_ny_(&layout, ny);
  // const GF3D2<CCTK_REAL> gf_nz_(&layout, nz);

  // const GF3D2<CCTK_REAL> gf_mret_(&layout, mret);
  // const GF3D2<CCTK_REAL> gf_mrex_(&layout, mrex);
  // const GF3D2<CCTK_REAL> gf_mrey_(&layout, mrey);
  // const GF3D2<CCTK_REAL> gf_mrez_(&layout, mrez);

  // const GF3D2<CCTK_REAL> gf_mimt_(&layout, mimt);
  // const GF3D2<CCTK_REAL> gf_mimx_(&layout, mimx);
  // const GF3D2<CCTK_REAL> gf_mimy_(&layout, mimy);
  // const GF3D2<CCTK_REAL> gf_mimz_(&layout, mimz);

  // //

  // const GF3D2<CCTK_REAL> gf_Lambda_(&layout, Lambda);
  // const GF3D2<CCTK_REAL> gf_Phi00_(&layout, Phi00);
  // const GF3D2<CCTK_REAL> gf_Phi11_(&layout, Phi11);
  // const GF3D2<CCTK_REAL> gf_Phi22_(&layout, Phi22);
  // const GF3D2<CCTK_REAL> gf_Phi10re_(&layout, Phi10re);
  // const GF3D2<CCTK_REAL> gf_Phi10im_(&layout, Phi10im);
  // const GF3D2<CCTK_REAL> gf_Phi20re_(&layout, Phi20re);
  // const GF3D2<CCTK_REAL> gf_Phi20im_(&layout, Phi20im);
  // const GF3D2<CCTK_REAL> gf_Phi21re_(&layout, Phi21re);
  // const GF3D2<CCTK_REAL> gf_Phi21im_(&layout, Phi21im);

  // //

  const GF3D2<CCTK_REAL> gf_Psi0re_(&layout, Psi0re);
  const GF3D2<CCTK_REAL> gf_Psi0im_(&layout, Psi0im);
  const GF3D2<CCTK_REAL> gf_Psi1re_(&layout, Psi1re);
  const GF3D2<CCTK_REAL> gf_Psi1im_(&layout, Psi1im);
  const GF3D2<CCTK_REAL> gf_Psi2re_(&layout, Psi2re);
  const GF3D2<CCTK_REAL> gf_Psi2im_(&layout, Psi2im);
  const GF3D2<CCTK_REAL> gf_Psi3re_(&layout, Psi3re);
  const GF3D2<CCTK_REAL> gf_Psi3im_(&layout, Psi3im);
  const GF3D2<CCTK_REAL> gf_Psi4re_(&layout, Psi4re);
  const GF3D2<CCTK_REAL> gf_Psi4im_(&layout, Psi4im);

  // //

  // const GF3D2<CCTK_REAL> gf_npkappare_(&layout, npkappare);
  // const GF3D2<CCTK_REAL> gf_npkappaim_(&layout, npkappaim);
  // const GF3D2<CCTK_REAL> gf_npsigmare_(&layout, npsigmare);
  // const GF3D2<CCTK_REAL> gf_npsigmaim_(&layout, npsigmaim);
  // const GF3D2<CCTK_REAL> gf_nprhore_(&layout, nprhore);
  // const GF3D2<CCTK_REAL> gf_nprhoim_(&layout, nprhoim);
  // const GF3D2<CCTK_REAL> gf_nptaure_(&layout, nptaure);
  // const GF3D2<CCTK_REAL> gf_nptauim_(&layout, nptauim);
  // const GF3D2<CCTK_REAL> gf_npepsilonre_(&layout, npepsilonre);
  // const GF3D2<CCTK_REAL> gf_npepsilonim_(&layout, npepsilonim);
  // const GF3D2<CCTK_REAL> gf_npbetare_(&layout, npbetare);
  // const GF3D2<CCTK_REAL> gf_npbetaim_(&layout, npbetaim);
  // const GF3D2<CCTK_REAL> gf_npalphare_(&layout, npalphare);
  // const GF3D2<CCTK_REAL> gf_npalphaim_(&layout, npalphaim);
  // const GF3D2<CCTK_REAL> gf_npgammare_(&layout, npgammare);
  // const GF3D2<CCTK_REAL> gf_npgammaim_(&layout, npgammaim);
  // const GF3D2<CCTK_REAL> gf_nppire_(&layout, nppire);
  // const GF3D2<CCTK_REAL> gf_nppiim_(&layout, nppiim);
  // const GF3D2<CCTK_REAL> gf_npmure_(&layout, npmure);
  // const GF3D2<CCTK_REAL> gf_npmuim_(&layout, npmuim);
  // const GF3D2<CCTK_REAL> gf_nplambdare_(&layout, nplambdare);
  // const GF3D2<CCTK_REAL> gf_nplambdaim_(&layout, nplambdaim);
  // const GF3D2<CCTK_REAL> gf_npnure_(&layout, npnure);
  // const GF3D2<CCTK_REAL> gf_npnuim_(&layout, npnuim);

  //

  loop_int<0, 0, 0>(
      cctkGH, [&](const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Load and calculate

        const vec3<CCTK_REAL, UP> coord3{p.x, p.y, p.z};

        const weyl_vars<CCTK_REAL> vars(
            cctk_time, coord3,                                 //
            gf_gamma_(p.I), gf_alpha_(p.I), gf_beta_(p.I),     //
            gf_k_(p.I), gf_dtalpha_(p.I), gf_dtbeta_(p.I),     //
            gf_dgamma_(p.I), gf_dalpha_(p.I), gf_dbeta_(p.I),  //
            gf_dtk_(p.I), gf_dt2alpha_(p.I), gf_dt2beta_(p.I), //
            gf_dk_(p.I), gf_ddtalpha_(p.I), gf_ddtbeta_(p.I),  //
            gf_ddgamma_(p.I), gf_ddalpha_(p.I), gf_ddbeta_(p.I));

        // Store
        vars.g.store(gf_g4tt_, gf_g4tx_, gf_g4ty_, gf_g4tz_, gf_g4xx_, gf_g4xy_,
                     gf_g4xz_, gf_g4yy_, gf_g4yz_, gf_g4zz_, p.I);

        // gf_Gamma4ttt_(p.I) = vars.Gamma(0)(0, 0);
        // gf_Gamma4ttx_(p.I) = vars.Gamma(0)(0, 1);
        // gf_Gamma4tty_(p.I) = vars.Gamma(0)(0, 2);
        // gf_Gamma4ttz_(p.I) = vars.Gamma(0)(0, 3);
        // gf_Gamma4txx_(p.I) = vars.Gamma(0)(1, 1);
        // gf_Gamma4txy_(p.I) = vars.Gamma(0)(1, 2);
        // gf_Gamma4txz_(p.I) = vars.Gamma(0)(1, 3);
        // gf_Gamma4tyy_(p.I) = vars.Gamma(0)(2, 2);
        // gf_Gamma4tyz_(p.I) = vars.Gamma(0)(2, 3);
        // gf_Gamma4tzz_(p.I) = vars.Gamma(0)(3, 3);

        // gf_Gamma4xtt_(p.I) = vars.Gamma(1)(0, 0);
        // gf_Gamma4xtx_(p.I) = vars.Gamma(1)(0, 1);
        // gf_Gamma4xty_(p.I) = vars.Gamma(1)(0, 2);
        // gf_Gamma4xtz_(p.I) = vars.Gamma(1)(0, 3);
        // gf_Gamma4xxx_(p.I) = vars.Gamma(1)(1, 1);
        // gf_Gamma4xxy_(p.I) = vars.Gamma(1)(1, 2);
        // gf_Gamma4xxz_(p.I) = vars.Gamma(1)(1, 3);
        // gf_Gamma4xyy_(p.I) = vars.Gamma(1)(2, 2);
        // gf_Gamma4xyz_(p.I) = vars.Gamma(1)(2, 3);
        // gf_Gamma4xzz_(p.I) = vars.Gamma(1)(3, 3);

        // gf_Gamma4ytt_(p.I) = vars.Gamma(2)(0, 0);
        // gf_Gamma4ytx_(p.I) = vars.Gamma(2)(0, 1);
        // gf_Gamma4yty_(p.I) = vars.Gamma(2)(0, 2);
        // gf_Gamma4ytz_(p.I) = vars.Gamma(2)(0, 3);
        // gf_Gamma4yxx_(p.I) = vars.Gamma(2)(1, 1);
        // gf_Gamma4yxy_(p.I) = vars.Gamma(2)(1, 2);
        // gf_Gamma4yxz_(p.I) = vars.Gamma(2)(1, 3);
        // gf_Gamma4yyy_(p.I) = vars.Gamma(2)(2, 2);
        // gf_Gamma4yyz_(p.I) = vars.Gamma(2)(2, 3);
        // gf_Gamma4yzz_(p.I) = vars.Gamma(2)(3, 3);

        // gf_Gamma4ztt_(p.I) = vars.Gamma(3)(0, 0);
        // gf_Gamma4ztx_(p.I) = vars.Gamma(3)(0, 1);
        // gf_Gamma4zty_(p.I) = vars.Gamma(3)(0, 2);
        // gf_Gamma4ztz_(p.I) = vars.Gamma(3)(0, 3);
        // gf_Gamma4zxx_(p.I) = vars.Gamma(3)(1, 1);
        // gf_Gamma4zxy_(p.I) = vars.Gamma(3)(1, 2);
        // gf_Gamma4zxz_(p.I) = vars.Gamma(3)(1, 3);
        // gf_Gamma4zyy_(p.I) = vars.Gamma(3)(2, 2);
        // gf_Gamma4zyz_(p.I) = vars.Gamma(3)(2, 3);
        // gf_Gamma4zzz_(p.I) = vars.Gamma(3)(3, 3);

        // gf_rm4txtx_(p.I) = vars.Rm(0, 1)(0, 1);
        // gf_rm4txty_(p.I) = vars.Rm(0, 1)(0, 2);
        // gf_rm4txtz_(p.I) = vars.Rm(0, 1)(0, 3);
        // gf_rm4txxy_(p.I) = vars.Rm(0, 1)(1, 2);
        // gf_rm4txxz_(p.I) = vars.Rm(0, 1)(1, 3);
        // gf_rm4txyz_(p.I) = vars.Rm(0, 1)(2, 3);

        // gf_rm4tyty_(p.I) = vars.Rm(0, 2)(0, 2);
        // gf_rm4tytz_(p.I) = vars.Rm(0, 2)(0, 3);
        // gf_rm4tyxy_(p.I) = vars.Rm(0, 2)(1, 2);
        // gf_rm4tyxz_(p.I) = vars.Rm(0, 2)(1, 3);
        // gf_rm4tyyz_(p.I) = vars.Rm(0, 2)(2, 3);

        // gf_rm4tztz_(p.I) = vars.Rm(0, 3)(0, 3);
        // gf_rm4tzxy_(p.I) = vars.Rm(0, 3)(1, 2);
        // gf_rm4tzxz_(p.I) = vars.Rm(0, 3)(1, 3);
        // gf_rm4tzyz_(p.I) = vars.Rm(0, 3)(2, 3);

        // gf_rm4xyxy_(p.I) = vars.Rm(1, 2)(1, 2);
        // gf_rm4xyxz_(p.I) = vars.Rm(1, 2)(1, 3);
        // gf_rm4xyyz_(p.I) = vars.Rm(1, 2)(2, 3);

        // gf_rm4xzxz_(p.I) = vars.Rm(1, 3)(1, 3);
        // gf_rm4xzyz_(p.I) = vars.Rm(1, 3)(2, 3);

        // gf_rm4yzyz_(p.I) = vars.Rm(2, 3)(2, 3);

        // vars.R.store(gf_r4tt_, gf_r4tx_, gf_r4ty_, gf_r4tz_, gf_r4xx_,
        // gf_r4xy_,
        //              gf_r4xz_, gf_r4yy_, gf_r4yz_, gf_r4zz_, p.I);

        // gf_rsc4_(p.I) = vars.Rsc;

        // gf_c4txtx_(p.I) = vars.C(0, 1)(0, 1);
        // gf_c4txty_(p.I) = vars.C(0, 1)(0, 2);
        // gf_c4txtz_(p.I) = vars.C(0, 1)(0, 3);
        // gf_c4txxy_(p.I) = vars.C(0, 1)(1, 2);
        // gf_c4txxz_(p.I) = vars.C(0, 1)(1, 3);
        // gf_c4txyz_(p.I) = vars.C(0, 1)(2, 3);

        // gf_c4tyty_(p.I) = vars.C(0, 2)(0, 2);
        // gf_c4tytz_(p.I) = vars.C(0, 2)(0, 3);
        // gf_c4tyxy_(p.I) = vars.C(0, 2)(1, 2);
        // gf_c4tyxz_(p.I) = vars.C(0, 2)(1, 3);
        // gf_c4tyyz_(p.I) = vars.C(0, 2)(2, 3);

        // gf_c4tztz_(p.I) = vars.C(0, 3)(0, 3);
        // gf_c4tzxy_(p.I) = vars.C(0, 3)(1, 2);
        // gf_c4tzxz_(p.I) = vars.C(0, 3)(1, 3);
        // gf_c4tzyz_(p.I) = vars.C(0, 3)(2, 3);

        // gf_c4xyxy_(p.I) = vars.C(1, 2)(1, 2);
        // gf_c4xyxz_(p.I) = vars.C(1, 2)(1, 3);
        // gf_c4xyyz_(p.I) = vars.C(1, 2)(2, 3);

        // gf_c4xzxz_(p.I) = vars.C(1, 3)(1, 3);
        // gf_c4xzyz_(p.I) = vars.C(1, 3)(2, 3);

        // gf_c4yzyz_(p.I) = vars.C(2, 3)(2, 3);

        // vars.l.store(gf_lt_, gf_lx_, gf_ly_, gf_lz_, p.I);
        // vars.n.store(gf_nt_, gf_nx_, gf_ny_, gf_nz_, p.I);
        // gf_mret_(p.I) = real(vars.m(0));
        // gf_mrex_(p.I) = real(vars.m(1));
        // gf_mrey_(p.I) = real(vars.m(2));
        // gf_mrez_(p.I) = real(vars.m(3));
        // gf_mimt_(p.I) = imag(vars.m(0));
        // gf_mimx_(p.I) = imag(vars.m(1));
        // gf_mimy_(p.I) = imag(vars.m(2));
        // gf_mimz_(p.I) = imag(vars.m(3));

        // gf_Lambda_(p.I) = vars.Lambda;
        // gf_Phi00_(p.I) = vars.Phi00;
        // gf_Phi11_(p.I) = vars.Phi11;
        // gf_Phi22_(p.I) = vars.Phi22;
        // gf_Phi10re_(p.I) = real(vars.Phi10);
        // gf_Phi10im_(p.I) = imag(vars.Phi10);
        // gf_Phi20re_(p.I) = real(vars.Phi20);
        // gf_Phi20im_(p.I) = imag(vars.Phi20);
        // gf_Phi21re_(p.I) = real(vars.Phi21);
        // gf_Phi21im_(p.I) = imag(vars.Phi21);

        gf_Psi0re_(p.I) = real(vars.Psi0);
        gf_Psi0im_(p.I) = imag(vars.Psi0);
        gf_Psi1re_(p.I) = real(vars.Psi1);
        gf_Psi1im_(p.I) = imag(vars.Psi1);
        gf_Psi2re_(p.I) = real(vars.Psi2);
        gf_Psi2im_(p.I) = imag(vars.Psi2);
        gf_Psi3re_(p.I) = real(vars.Psi3);
        gf_Psi3im_(p.I) = imag(vars.Psi3);
        gf_Psi4re_(p.I) = real(vars.Psi4);
        gf_Psi4im_(p.I) = imag(vars.Psi4);

        // gf_npkappare_(p.I) = real(vars.npkappa);
        // gf_npkappaim_(p.I) = imag(vars.npkappa);
        // gf_npsigmare_(p.I) = real(vars.npsigma);
        // gf_npsigmaim_(p.I) = imag(vars.npsigma);
        // gf_nprhore_(p.I) = real(vars.nprho);
        // gf_nprhoim_(p.I) = imag(vars.nprho);
        // gf_nptaure_(p.I) = real(vars.nptau);
        // gf_nptauim_(p.I) = imag(vars.nptau);
        // gf_npepsilonre_(p.I) = real(vars.npepsilon);
        // gf_npepsilonim_(p.I) = imag(vars.npepsilon);
        // gf_npbetare_(p.I) = real(vars.npbeta);
        // gf_npbetaim_(p.I) = imag(vars.npbeta);
        // gf_npalphare_(p.I) = real(vars.npalpha);
        // gf_npalphaim_(p.I) = imag(vars.npalpha);
        // gf_npgammare_(p.I) = real(vars.npgamma);
        // gf_npgammaim_(p.I) = imag(vars.npgamma);
        // gf_nppire_(p.I) = real(vars.nppi);
        // gf_nppiim_(p.I) = imag(vars.nppi);
        // gf_npmure_(p.I) = real(vars.npmu);
        // gf_npmuim_(p.I) = imag(vars.npmu);
        // gf_nplambdare_(p.I) = real(vars.nplambda);
        // gf_nplambdaim_(p.I) = imag(vars.nplambda);
        // gf_npnure_(p.I) = real(vars.npnu);
        // gf_npnuim_(p.I) = imag(vars.npnu);
      });
}

} // namespace Weyl
