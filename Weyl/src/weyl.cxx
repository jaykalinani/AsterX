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
  const array<int, dim> nghostzones = {cctk_nghostzones[0], cctk_nghostzones[1],
                                       cctk_nghostzones[2]};
  const array<int, dim> noghosts = {0, 0, 0};

  const mat3<GF3D1<const CCTK_REAL>, DN, DN> gf_gamma_(
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, gxx),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, gxy),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, gxz),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, gyy),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, gyz),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, gzz));

  const GF3D1<const CCTK_REAL> gf_alpha_(cctkGH, indextype, nghostzones, alp);

  const vec3<GF3D1<const CCTK_REAL>, UP> gf_beta_(
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, betax),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, betay),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, betaz));

  const mat3<GF3D1<const CCTK_REAL>, DN, DN> gf_k_(
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, kxx),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, kxy),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, kxz),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, kyy),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, kyz),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, kzz));

  const GF3D1<const CCTK_REAL> gf_dtalpha_(cctkGH, indextype, nghostzones,
                                           dtalp);

  const vec3<GF3D1<const CCTK_REAL>, UP> gf_dtbeta_(
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, dtbetax),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, dtbetay),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, dtbetaz));

  const mat3<GF3D1<const CCTK_REAL>, DN, DN> gf_dtk_(
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, dtkxx),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, dtkxy),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, dtkxz),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, dtkyy),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, dtkyz),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, dtkzz));

  const GF3D1<const CCTK_REAL> gf_dt2alpha_(cctkGH, indextype, nghostzones,
                                            dt2alp);

  const vec3<GF3D1<const CCTK_REAL>, UP> gf_dt2beta_(
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, dt2betax),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, dt2betay),
      GF3D1<const CCTK_REAL>(cctkGH, indextype, nghostzones, dt2betaz));

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

  const mat3<vec3<GF3D1<CCTK_REAL>, DN>, DN, DN> gf_dgamma_(make_mat_vec_gf);
  const mat3<mat3<GF3D1<CCTK_REAL>, DN, DN>, DN, DN> gf_ddgamma_(
      make_mat_mat_gf);
  calc_derivs2(cctkGH, gf_gamma_, gf_dgamma_, gf_ddgamma_);

  const vec3<GF3D1<CCTK_REAL>, DN> gf_dalpha_(make_vec_gf);
  const mat3<GF3D1<CCTK_REAL>, DN, DN> gf_ddalpha_(make_mat_gf);
  calc_derivs2(cctkGH, gf_alpha_, gf_dalpha_, gf_ddalpha_);

  const vec3<vec3<GF3D1<CCTK_REAL>, DN>, UP> gf_dbeta_(make_vec_vec_gf);
  const vec3<mat3<GF3D1<CCTK_REAL>, DN, DN>, UP> gf_ddbeta_(make_vec_mat_gf);
  calc_derivs2(cctkGH, gf_beta_, gf_dbeta_, gf_ddbeta_);

  const mat3<vec3<GF3D1<CCTK_REAL>, DN>, DN, DN> gf_dk_(make_mat_vec_gf);
  calc_derivs(cctkGH, gf_k_, gf_dk_);

  const vec3<GF3D1<CCTK_REAL>, DN> gf_ddtalpha_(make_vec_gf);
  calc_derivs(cctkGH, gf_dtalpha_, gf_ddtalpha_);

  const vec3<vec3<GF3D1<CCTK_REAL>, DN>, UP> gf_ddtbeta_(make_vec_vec_gf);
  calc_derivs(cctkGH, gf_dtbeta_, gf_ddtbeta_);

  //

  const GF3D1<CCTK_REAL> gf_g4tt_(cctkGH, indextype, nghostzones, g4tt);
  const GF3D1<CCTK_REAL> gf_g4tx_(cctkGH, indextype, nghostzones, g4tx);
  const GF3D1<CCTK_REAL> gf_g4ty_(cctkGH, indextype, nghostzones, g4ty);
  const GF3D1<CCTK_REAL> gf_g4tz_(cctkGH, indextype, nghostzones, g4tz);
  const GF3D1<CCTK_REAL> gf_g4xx_(cctkGH, indextype, nghostzones, g4xx);
  const GF3D1<CCTK_REAL> gf_g4xy_(cctkGH, indextype, nghostzones, g4xy);
  const GF3D1<CCTK_REAL> gf_g4xz_(cctkGH, indextype, nghostzones, g4xz);
  const GF3D1<CCTK_REAL> gf_g4yy_(cctkGH, indextype, nghostzones, g4yy);
  const GF3D1<CCTK_REAL> gf_g4yz_(cctkGH, indextype, nghostzones, g4yz);
  const GF3D1<CCTK_REAL> gf_g4zz_(cctkGH, indextype, nghostzones, g4zz);

  const GF3D1<CCTK_REAL> gf_Gamma4ttt_(cctkGH, indextype, nghostzones,
                                       Gamma4ttt);
  const GF3D1<CCTK_REAL> gf_Gamma4ttx_(cctkGH, indextype, nghostzones,
                                       Gamma4ttx);
  const GF3D1<CCTK_REAL> gf_Gamma4tty_(cctkGH, indextype, nghostzones,
                                       Gamma4tty);
  const GF3D1<CCTK_REAL> gf_Gamma4ttz_(cctkGH, indextype, nghostzones,
                                       Gamma4ttz);
  const GF3D1<CCTK_REAL> gf_Gamma4txx_(cctkGH, indextype, nghostzones,
                                       Gamma4txx);
  const GF3D1<CCTK_REAL> gf_Gamma4txy_(cctkGH, indextype, nghostzones,
                                       Gamma4txy);
  const GF3D1<CCTK_REAL> gf_Gamma4txz_(cctkGH, indextype, nghostzones,
                                       Gamma4txz);
  const GF3D1<CCTK_REAL> gf_Gamma4tyy_(cctkGH, indextype, nghostzones,
                                       Gamma4tyy);
  const GF3D1<CCTK_REAL> gf_Gamma4tyz_(cctkGH, indextype, nghostzones,
                                       Gamma4tyz);
  const GF3D1<CCTK_REAL> gf_Gamma4tzz_(cctkGH, indextype, nghostzones,
                                       Gamma4tzz);

  const GF3D1<CCTK_REAL> gf_Gamma4xtt_(cctkGH, indextype, nghostzones,
                                       Gamma4xtt);
  const GF3D1<CCTK_REAL> gf_Gamma4xtx_(cctkGH, indextype, nghostzones,
                                       Gamma4xtx);
  const GF3D1<CCTK_REAL> gf_Gamma4xty_(cctkGH, indextype, nghostzones,
                                       Gamma4xty);
  const GF3D1<CCTK_REAL> gf_Gamma4xtz_(cctkGH, indextype, nghostzones,
                                       Gamma4xtz);
  const GF3D1<CCTK_REAL> gf_Gamma4xxx_(cctkGH, indextype, nghostzones,
                                       Gamma4xxx);
  const GF3D1<CCTK_REAL> gf_Gamma4xxy_(cctkGH, indextype, nghostzones,
                                       Gamma4xxy);
  const GF3D1<CCTK_REAL> gf_Gamma4xxz_(cctkGH, indextype, nghostzones,
                                       Gamma4xxz);
  const GF3D1<CCTK_REAL> gf_Gamma4xyy_(cctkGH, indextype, nghostzones,
                                       Gamma4xyy);
  const GF3D1<CCTK_REAL> gf_Gamma4xyz_(cctkGH, indextype, nghostzones,
                                       Gamma4xyz);
  const GF3D1<CCTK_REAL> gf_Gamma4xzz_(cctkGH, indextype, nghostzones,
                                       Gamma4xzz);

  const GF3D1<CCTK_REAL> gf_Gamma4ytt_(cctkGH, indextype, nghostzones,
                                       Gamma4ytt);
  const GF3D1<CCTK_REAL> gf_Gamma4ytx_(cctkGH, indextype, nghostzones,
                                       Gamma4ytx);
  const GF3D1<CCTK_REAL> gf_Gamma4yty_(cctkGH, indextype, nghostzones,
                                       Gamma4yty);
  const GF3D1<CCTK_REAL> gf_Gamma4ytz_(cctkGH, indextype, nghostzones,
                                       Gamma4ytz);
  const GF3D1<CCTK_REAL> gf_Gamma4yxx_(cctkGH, indextype, nghostzones,
                                       Gamma4yxx);
  const GF3D1<CCTK_REAL> gf_Gamma4yxy_(cctkGH, indextype, nghostzones,
                                       Gamma4yxy);
  const GF3D1<CCTK_REAL> gf_Gamma4yxz_(cctkGH, indextype, nghostzones,
                                       Gamma4yxz);
  const GF3D1<CCTK_REAL> gf_Gamma4yyy_(cctkGH, indextype, nghostzones,
                                       Gamma4yyy);
  const GF3D1<CCTK_REAL> gf_Gamma4yyz_(cctkGH, indextype, nghostzones,
                                       Gamma4yyz);
  const GF3D1<CCTK_REAL> gf_Gamma4yzz_(cctkGH, indextype, nghostzones,
                                       Gamma4yzz);

  const GF3D1<CCTK_REAL> gf_Gamma4ztt_(cctkGH, indextype, nghostzones,
                                       Gamma4ztt);
  const GF3D1<CCTK_REAL> gf_Gamma4ztx_(cctkGH, indextype, nghostzones,
                                       Gamma4ztx);
  const GF3D1<CCTK_REAL> gf_Gamma4zty_(cctkGH, indextype, nghostzones,
                                       Gamma4zty);
  const GF3D1<CCTK_REAL> gf_Gamma4ztz_(cctkGH, indextype, nghostzones,
                                       Gamma4ztz);
  const GF3D1<CCTK_REAL> gf_Gamma4zxx_(cctkGH, indextype, nghostzones,
                                       Gamma4zxx);
  const GF3D1<CCTK_REAL> gf_Gamma4zxy_(cctkGH, indextype, nghostzones,
                                       Gamma4zxy);
  const GF3D1<CCTK_REAL> gf_Gamma4zxz_(cctkGH, indextype, nghostzones,
                                       Gamma4zxz);
  const GF3D1<CCTK_REAL> gf_Gamma4zyy_(cctkGH, indextype, nghostzones,
                                       Gamma4zyy);
  const GF3D1<CCTK_REAL> gf_Gamma4zyz_(cctkGH, indextype, nghostzones,
                                       Gamma4zyz);
  const GF3D1<CCTK_REAL> gf_Gamma4zzz_(cctkGH, indextype, nghostzones,
                                       Gamma4zzz);

  const GF3D1<CCTK_REAL> gf_rm4txtx_(cctkGH, indextype, nghostzones, rm4txtx);
  const GF3D1<CCTK_REAL> gf_rm4txty_(cctkGH, indextype, nghostzones, rm4txty);
  const GF3D1<CCTK_REAL> gf_rm4txtz_(cctkGH, indextype, nghostzones, rm4txtz);
  const GF3D1<CCTK_REAL> gf_rm4txxy_(cctkGH, indextype, nghostzones, rm4txxy);
  const GF3D1<CCTK_REAL> gf_rm4txxz_(cctkGH, indextype, nghostzones, rm4txxz);
  const GF3D1<CCTK_REAL> gf_rm4txyz_(cctkGH, indextype, nghostzones, rm4txyz);

  const GF3D1<CCTK_REAL> gf_rm4tyty_(cctkGH, indextype, nghostzones, rm4tyty);
  const GF3D1<CCTK_REAL> gf_rm4tytz_(cctkGH, indextype, nghostzones, rm4tytz);
  const GF3D1<CCTK_REAL> gf_rm4tyxy_(cctkGH, indextype, nghostzones, rm4tyxy);
  const GF3D1<CCTK_REAL> gf_rm4tyxz_(cctkGH, indextype, nghostzones, rm4tyxz);
  const GF3D1<CCTK_REAL> gf_rm4tyyz_(cctkGH, indextype, nghostzones, rm4tyyz);

  const GF3D1<CCTK_REAL> gf_rm4tztz_(cctkGH, indextype, nghostzones, rm4tztz);
  const GF3D1<CCTK_REAL> gf_rm4tzxy_(cctkGH, indextype, nghostzones, rm4tzxy);
  const GF3D1<CCTK_REAL> gf_rm4tzxz_(cctkGH, indextype, nghostzones, rm4tzxz);
  const GF3D1<CCTK_REAL> gf_rm4tzyz_(cctkGH, indextype, nghostzones, rm4tzyz);

  const GF3D1<CCTK_REAL> gf_rm4xyxy_(cctkGH, indextype, nghostzones, rm4xyxy);
  const GF3D1<CCTK_REAL> gf_rm4xyxz_(cctkGH, indextype, nghostzones, rm4xyxz);
  const GF3D1<CCTK_REAL> gf_rm4xyyz_(cctkGH, indextype, nghostzones, rm4xyyz);

  const GF3D1<CCTK_REAL> gf_rm4xzxz_(cctkGH, indextype, nghostzones, rm4xzxz);
  const GF3D1<CCTK_REAL> gf_rm4xzyz_(cctkGH, indextype, nghostzones, rm4xzyz);

  const GF3D1<CCTK_REAL> gf_rm4yzyz_(cctkGH, indextype, nghostzones, rm4yzyz);

  const GF3D1<CCTK_REAL> gf_r4tt_(cctkGH, indextype, nghostzones, r4tt);
  const GF3D1<CCTK_REAL> gf_r4tx_(cctkGH, indextype, nghostzones, r4tx);
  const GF3D1<CCTK_REAL> gf_r4ty_(cctkGH, indextype, nghostzones, r4ty);
  const GF3D1<CCTK_REAL> gf_r4tz_(cctkGH, indextype, nghostzones, r4tz);
  const GF3D1<CCTK_REAL> gf_r4xx_(cctkGH, indextype, nghostzones, r4xx);
  const GF3D1<CCTK_REAL> gf_r4xy_(cctkGH, indextype, nghostzones, r4xy);
  const GF3D1<CCTK_REAL> gf_r4xz_(cctkGH, indextype, nghostzones, r4xz);
  const GF3D1<CCTK_REAL> gf_r4yy_(cctkGH, indextype, nghostzones, r4yy);
  const GF3D1<CCTK_REAL> gf_r4yz_(cctkGH, indextype, nghostzones, r4yz);
  const GF3D1<CCTK_REAL> gf_r4zz_(cctkGH, indextype, nghostzones, r4zz);

  const GF3D1<CCTK_REAL> gf_rsc4_(cctkGH, indextype, nghostzones, rsc4);

  const GF3D1<CCTK_REAL> gf_c4txtx_(cctkGH, indextype, nghostzones, c4txtx);
  const GF3D1<CCTK_REAL> gf_c4txty_(cctkGH, indextype, nghostzones, c4txty);
  const GF3D1<CCTK_REAL> gf_c4txtz_(cctkGH, indextype, nghostzones, c4txtz);
  const GF3D1<CCTK_REAL> gf_c4txxy_(cctkGH, indextype, nghostzones, c4txxy);
  const GF3D1<CCTK_REAL> gf_c4txxz_(cctkGH, indextype, nghostzones, c4txxz);
  const GF3D1<CCTK_REAL> gf_c4txyz_(cctkGH, indextype, nghostzones, c4txyz);

  const GF3D1<CCTK_REAL> gf_c4tyty_(cctkGH, indextype, nghostzones, c4tyty);
  const GF3D1<CCTK_REAL> gf_c4tytz_(cctkGH, indextype, nghostzones, c4tytz);
  const GF3D1<CCTK_REAL> gf_c4tyxy_(cctkGH, indextype, nghostzones, c4tyxy);
  const GF3D1<CCTK_REAL> gf_c4tyxz_(cctkGH, indextype, nghostzones, c4tyxz);
  const GF3D1<CCTK_REAL> gf_c4tyyz_(cctkGH, indextype, nghostzones, c4tyyz);

  const GF3D1<CCTK_REAL> gf_c4tztz_(cctkGH, indextype, nghostzones, c4tztz);
  const GF3D1<CCTK_REAL> gf_c4tzxy_(cctkGH, indextype, nghostzones, c4tzxy);
  const GF3D1<CCTK_REAL> gf_c4tzxz_(cctkGH, indextype, nghostzones, c4tzxz);
  const GF3D1<CCTK_REAL> gf_c4tzyz_(cctkGH, indextype, nghostzones, c4tzyz);

  const GF3D1<CCTK_REAL> gf_c4xyxy_(cctkGH, indextype, nghostzones, c4xyxy);
  const GF3D1<CCTK_REAL> gf_c4xyxz_(cctkGH, indextype, nghostzones, c4xyxz);
  const GF3D1<CCTK_REAL> gf_c4xyyz_(cctkGH, indextype, nghostzones, c4xyyz);

  const GF3D1<CCTK_REAL> gf_c4xzxz_(cctkGH, indextype, nghostzones, c4xzxz);
  const GF3D1<CCTK_REAL> gf_c4xzyz_(cctkGH, indextype, nghostzones, c4xzyz);

  const GF3D1<CCTK_REAL> gf_c4yzyz_(cctkGH, indextype, nghostzones, c4yzyz);

  //

  const GF3D1<CCTK_REAL> gf_lt_(cctkGH, indextype, nghostzones, lt);
  const GF3D1<CCTK_REAL> gf_lx_(cctkGH, indextype, nghostzones, lx);
  const GF3D1<CCTK_REAL> gf_ly_(cctkGH, indextype, nghostzones, ly);
  const GF3D1<CCTK_REAL> gf_lz_(cctkGH, indextype, nghostzones, lz);

  const GF3D1<CCTK_REAL> gf_nt_(cctkGH, indextype, nghostzones, nt);
  const GF3D1<CCTK_REAL> gf_nx_(cctkGH, indextype, nghostzones, nx);
  const GF3D1<CCTK_REAL> gf_ny_(cctkGH, indextype, nghostzones, ny);
  const GF3D1<CCTK_REAL> gf_nz_(cctkGH, indextype, nghostzones, nz);

  const GF3D1<CCTK_REAL> gf_mret_(cctkGH, indextype, nghostzones, mret);
  const GF3D1<CCTK_REAL> gf_mrex_(cctkGH, indextype, nghostzones, mrex);
  const GF3D1<CCTK_REAL> gf_mrey_(cctkGH, indextype, nghostzones, mrey);
  const GF3D1<CCTK_REAL> gf_mrez_(cctkGH, indextype, nghostzones, mrez);

  const GF3D1<CCTK_REAL> gf_mimt_(cctkGH, indextype, nghostzones, mimt);
  const GF3D1<CCTK_REAL> gf_mimx_(cctkGH, indextype, nghostzones, mimx);
  const GF3D1<CCTK_REAL> gf_mimy_(cctkGH, indextype, nghostzones, mimy);
  const GF3D1<CCTK_REAL> gf_mimz_(cctkGH, indextype, nghostzones, mimz);

  //

  const GF3D1<CCTK_REAL> gf_Lambda_(cctkGH, indextype, nghostzones, Lambda);
  const GF3D1<CCTK_REAL> gf_Phi00_(cctkGH, indextype, nghostzones, Phi00);
  const GF3D1<CCTK_REAL> gf_Phi11_(cctkGH, indextype, nghostzones, Phi11);
  const GF3D1<CCTK_REAL> gf_Phi22_(cctkGH, indextype, nghostzones, Phi22);
  const GF3D1<CCTK_REAL> gf_Phi10re_(cctkGH, indextype, nghostzones, Phi10re);
  const GF3D1<CCTK_REAL> gf_Phi10im_(cctkGH, indextype, nghostzones, Phi10im);
  const GF3D1<CCTK_REAL> gf_Phi20re_(cctkGH, indextype, nghostzones, Phi20re);
  const GF3D1<CCTK_REAL> gf_Phi20im_(cctkGH, indextype, nghostzones, Phi20im);
  const GF3D1<CCTK_REAL> gf_Phi21re_(cctkGH, indextype, nghostzones, Phi21re);
  const GF3D1<CCTK_REAL> gf_Phi21im_(cctkGH, indextype, nghostzones, Phi21im);

  //

  const GF3D1<CCTK_REAL> gf_Psi0re_(cctkGH, indextype, nghostzones, Psi0re);
  const GF3D1<CCTK_REAL> gf_Psi0im_(cctkGH, indextype, nghostzones, Psi0im);
  const GF3D1<CCTK_REAL> gf_Psi1re_(cctkGH, indextype, nghostzones, Psi1re);
  const GF3D1<CCTK_REAL> gf_Psi1im_(cctkGH, indextype, nghostzones, Psi1im);
  const GF3D1<CCTK_REAL> gf_Psi2re_(cctkGH, indextype, nghostzones, Psi2re);
  const GF3D1<CCTK_REAL> gf_Psi2im_(cctkGH, indextype, nghostzones, Psi2im);
  const GF3D1<CCTK_REAL> gf_Psi3re_(cctkGH, indextype, nghostzones, Psi3re);
  const GF3D1<CCTK_REAL> gf_Psi3im_(cctkGH, indextype, nghostzones, Psi3im);
  const GF3D1<CCTK_REAL> gf_Psi4re_(cctkGH, indextype, nghostzones, Psi4re);
  const GF3D1<CCTK_REAL> gf_Psi4im_(cctkGH, indextype, nghostzones, Psi4im);

  //

  const GF3D1<CCTK_REAL> gf_npkappare_(cctkGH, indextype, nghostzones,
                                       npkappare);
  const GF3D1<CCTK_REAL> gf_npkappaim_(cctkGH, indextype, nghostzones,
                                       npkappaim);
  const GF3D1<CCTK_REAL> gf_npsigmare_(cctkGH, indextype, nghostzones,
                                       npsigmare);
  const GF3D1<CCTK_REAL> gf_npsigmaim_(cctkGH, indextype, nghostzones,
                                       npsigmaim);
  const GF3D1<CCTK_REAL> gf_nprhore_(cctkGH, indextype, nghostzones, nprhore);
  const GF3D1<CCTK_REAL> gf_nprhoim_(cctkGH, indextype, nghostzones, nprhoim);
  const GF3D1<CCTK_REAL> gf_nptaure_(cctkGH, indextype, nghostzones, nptaure);
  const GF3D1<CCTK_REAL> gf_nptauim_(cctkGH, indextype, nghostzones, nptauim);
  const GF3D1<CCTK_REAL> gf_npepsilonre_(cctkGH, indextype, nghostzones,
                                         npepsilonre);
  const GF3D1<CCTK_REAL> gf_npepsilonim_(cctkGH, indextype, nghostzones,
                                         npepsilonim);
  const GF3D1<CCTK_REAL> gf_npbetare_(cctkGH, indextype, nghostzones, npbetare);
  const GF3D1<CCTK_REAL> gf_npbetaim_(cctkGH, indextype, nghostzones, npbetaim);
  const GF3D1<CCTK_REAL> gf_npalphare_(cctkGH, indextype, nghostzones,
                                       npalphare);
  const GF3D1<CCTK_REAL> gf_npalphaim_(cctkGH, indextype, nghostzones,
                                       npalphaim);
  const GF3D1<CCTK_REAL> gf_npgammare_(cctkGH, indextype, nghostzones,
                                       npgammare);
  const GF3D1<CCTK_REAL> gf_npgammaim_(cctkGH, indextype, nghostzones,
                                       npgammaim);
  const GF3D1<CCTK_REAL> gf_nppire_(cctkGH, indextype, nghostzones, nppire);
  const GF3D1<CCTK_REAL> gf_nppiim_(cctkGH, indextype, nghostzones, nppiim);
  const GF3D1<CCTK_REAL> gf_npmure_(cctkGH, indextype, nghostzones, npmure);
  const GF3D1<CCTK_REAL> gf_npmuim_(cctkGH, indextype, nghostzones, npmuim);
  const GF3D1<CCTK_REAL> gf_nplambdare_(cctkGH, indextype, nghostzones,
                                        nplambdare);
  const GF3D1<CCTK_REAL> gf_nplambdaim_(cctkGH, indextype, nghostzones,
                                        nplambdaim);
  const GF3D1<CCTK_REAL> gf_npnure_(cctkGH, indextype, nghostzones, npnure);
  const GF3D1<CCTK_REAL> gf_npnuim_(cctkGH, indextype, nghostzones, npnuim);

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

        gf_Gamma4ttt_(p.I) = vars.Gamma(0)(0, 0);
        gf_Gamma4ttx_(p.I) = vars.Gamma(0)(0, 1);
        gf_Gamma4tty_(p.I) = vars.Gamma(0)(0, 2);
        gf_Gamma4ttz_(p.I) = vars.Gamma(0)(0, 3);
        gf_Gamma4txx_(p.I) = vars.Gamma(0)(1, 1);
        gf_Gamma4txy_(p.I) = vars.Gamma(0)(1, 2);
        gf_Gamma4txz_(p.I) = vars.Gamma(0)(1, 3);
        gf_Gamma4tyy_(p.I) = vars.Gamma(0)(2, 2);
        gf_Gamma4tyz_(p.I) = vars.Gamma(0)(2, 3);
        gf_Gamma4tzz_(p.I) = vars.Gamma(0)(3, 3);

        gf_Gamma4xtt_(p.I) = vars.Gamma(1)(0, 0);
        gf_Gamma4xtx_(p.I) = vars.Gamma(1)(0, 1);
        gf_Gamma4xty_(p.I) = vars.Gamma(1)(0, 2);
        gf_Gamma4xtz_(p.I) = vars.Gamma(1)(0, 3);
        gf_Gamma4xxx_(p.I) = vars.Gamma(1)(1, 1);
        gf_Gamma4xxy_(p.I) = vars.Gamma(1)(1, 2);
        gf_Gamma4xxz_(p.I) = vars.Gamma(1)(1, 3);
        gf_Gamma4xyy_(p.I) = vars.Gamma(1)(2, 2);
        gf_Gamma4xyz_(p.I) = vars.Gamma(1)(2, 3);
        gf_Gamma4xzz_(p.I) = vars.Gamma(1)(3, 3);

        gf_Gamma4ytt_(p.I) = vars.Gamma(2)(0, 0);
        gf_Gamma4ytx_(p.I) = vars.Gamma(2)(0, 1);
        gf_Gamma4yty_(p.I) = vars.Gamma(2)(0, 2);
        gf_Gamma4ytz_(p.I) = vars.Gamma(2)(0, 3);
        gf_Gamma4yxx_(p.I) = vars.Gamma(2)(1, 1);
        gf_Gamma4yxy_(p.I) = vars.Gamma(2)(1, 2);
        gf_Gamma4yxz_(p.I) = vars.Gamma(2)(1, 3);
        gf_Gamma4yyy_(p.I) = vars.Gamma(2)(2, 2);
        gf_Gamma4yyz_(p.I) = vars.Gamma(2)(2, 3);
        gf_Gamma4yzz_(p.I) = vars.Gamma(2)(3, 3);

        gf_Gamma4ztt_(p.I) = vars.Gamma(3)(0, 0);
        gf_Gamma4ztx_(p.I) = vars.Gamma(3)(0, 1);
        gf_Gamma4zty_(p.I) = vars.Gamma(3)(0, 2);
        gf_Gamma4ztz_(p.I) = vars.Gamma(3)(0, 3);
        gf_Gamma4zxx_(p.I) = vars.Gamma(3)(1, 1);
        gf_Gamma4zxy_(p.I) = vars.Gamma(3)(1, 2);
        gf_Gamma4zxz_(p.I) = vars.Gamma(3)(1, 3);
        gf_Gamma4zyy_(p.I) = vars.Gamma(3)(2, 2);
        gf_Gamma4zyz_(p.I) = vars.Gamma(3)(2, 3);
        gf_Gamma4zzz_(p.I) = vars.Gamma(3)(3, 3);

        gf_rm4txtx_(p.I) = vars.Rm(0, 1)(0, 1);
        gf_rm4txty_(p.I) = vars.Rm(0, 1)(0, 2);
        gf_rm4txtz_(p.I) = vars.Rm(0, 1)(0, 3);
        gf_rm4txxy_(p.I) = vars.Rm(0, 1)(1, 2);
        gf_rm4txxz_(p.I) = vars.Rm(0, 1)(1, 3);
        gf_rm4txyz_(p.I) = vars.Rm(0, 1)(2, 3);

        gf_rm4tyty_(p.I) = vars.Rm(0, 2)(0, 2);
        gf_rm4tytz_(p.I) = vars.Rm(0, 2)(0, 3);
        gf_rm4tyxy_(p.I) = vars.Rm(0, 2)(1, 2);
        gf_rm4tyxz_(p.I) = vars.Rm(0, 2)(1, 3);
        gf_rm4tyyz_(p.I) = vars.Rm(0, 2)(2, 3);

        gf_rm4tztz_(p.I) = vars.Rm(0, 3)(0, 3);
        gf_rm4tzxy_(p.I) = vars.Rm(0, 3)(1, 2);
        gf_rm4tzxz_(p.I) = vars.Rm(0, 3)(1, 3);
        gf_rm4tzyz_(p.I) = vars.Rm(0, 3)(2, 3);

        gf_rm4xyxy_(p.I) = vars.Rm(1, 2)(1, 2);
        gf_rm4xyxz_(p.I) = vars.Rm(1, 2)(1, 3);
        gf_rm4xyyz_(p.I) = vars.Rm(1, 2)(2, 3);

        gf_rm4xzxz_(p.I) = vars.Rm(1, 3)(1, 3);
        gf_rm4xzyz_(p.I) = vars.Rm(1, 3)(2, 3);

        gf_rm4yzyz_(p.I) = vars.Rm(2, 3)(2, 3);

        vars.R.store(gf_r4tt_, gf_r4tx_, gf_r4ty_, gf_r4tz_, gf_r4xx_, gf_r4xy_,
                     gf_r4xz_, gf_r4yy_, gf_r4yz_, gf_r4zz_, p.I);

        gf_rsc4_(p.I) = vars.Rsc;

        gf_c4txtx_(p.I) = vars.C(0, 1)(0, 1);
        gf_c4txty_(p.I) = vars.C(0, 1)(0, 2);
        gf_c4txtz_(p.I) = vars.C(0, 1)(0, 3);
        gf_c4txxy_(p.I) = vars.C(0, 1)(1, 2);
        gf_c4txxz_(p.I) = vars.C(0, 1)(1, 3);
        gf_c4txyz_(p.I) = vars.C(0, 1)(2, 3);

        gf_c4tyty_(p.I) = vars.C(0, 2)(0, 2);
        gf_c4tytz_(p.I) = vars.C(0, 2)(0, 3);
        gf_c4tyxy_(p.I) = vars.C(0, 2)(1, 2);
        gf_c4tyxz_(p.I) = vars.C(0, 2)(1, 3);
        gf_c4tyyz_(p.I) = vars.C(0, 2)(2, 3);

        gf_c4tztz_(p.I) = vars.C(0, 3)(0, 3);
        gf_c4tzxy_(p.I) = vars.C(0, 3)(1, 2);
        gf_c4tzxz_(p.I) = vars.C(0, 3)(1, 3);
        gf_c4tzyz_(p.I) = vars.C(0, 3)(2, 3);

        gf_c4xyxy_(p.I) = vars.C(1, 2)(1, 2);
        gf_c4xyxz_(p.I) = vars.C(1, 2)(1, 3);
        gf_c4xyyz_(p.I) = vars.C(1, 2)(2, 3);

        gf_c4xzxz_(p.I) = vars.C(1, 3)(1, 3);
        gf_c4xzyz_(p.I) = vars.C(1, 3)(2, 3);

        gf_c4yzyz_(p.I) = vars.C(2, 3)(2, 3);

        vars.l.store(gf_lt_, gf_lx_, gf_ly_, gf_lz_, p.I);
        vars.n.store(gf_nt_, gf_nx_, gf_ny_, gf_nz_, p.I);
        gf_mret_(p.I) = real(vars.m(0));
        gf_mrex_(p.I) = real(vars.m(1));
        gf_mrey_(p.I) = real(vars.m(2));
        gf_mrez_(p.I) = real(vars.m(3));
        gf_mimt_(p.I) = imag(vars.m(0));
        gf_mimx_(p.I) = imag(vars.m(1));
        gf_mimy_(p.I) = imag(vars.m(2));
        gf_mimz_(p.I) = imag(vars.m(3));

        gf_Lambda_(p.I) = vars.Lambda;
        gf_Phi00_(p.I) = vars.Phi00;
        gf_Phi11_(p.I) = vars.Phi11;
        gf_Phi22_(p.I) = vars.Phi22;
        gf_Phi10re_(p.I) = real(vars.Phi10);
        gf_Phi10im_(p.I) = imag(vars.Phi10);
        gf_Phi20re_(p.I) = real(vars.Phi20);
        gf_Phi20im_(p.I) = imag(vars.Phi20);
        gf_Phi21re_(p.I) = real(vars.Phi21);
        gf_Phi21im_(p.I) = imag(vars.Phi21);

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

        gf_npkappare_(p.I) = real(vars.npkappa);
        gf_npkappaim_(p.I) = imag(vars.npkappa);
        gf_npsigmare_(p.I) = real(vars.npsigma);
        gf_npsigmaim_(p.I) = imag(vars.npsigma);
        gf_nprhore_(p.I) = real(vars.nprho);
        gf_nprhoim_(p.I) = imag(vars.nprho);
        gf_nptaure_(p.I) = real(vars.nptau);
        gf_nptauim_(p.I) = imag(vars.nptau);
        gf_npepsilonre_(p.I) = real(vars.npepsilon);
        gf_npepsilonim_(p.I) = imag(vars.npepsilon);
        gf_npbetare_(p.I) = real(vars.npbeta);
        gf_npbetaim_(p.I) = imag(vars.npbeta);
        gf_npalphare_(p.I) = real(vars.npalpha);
        gf_npalphaim_(p.I) = imag(vars.npalpha);
        gf_npgammare_(p.I) = real(vars.npgamma);
        gf_npgammaim_(p.I) = imag(vars.npgamma);
        gf_nppire_(p.I) = real(vars.nppi);
        gf_nppiim_(p.I) = imag(vars.nppi);
        gf_npmure_(p.I) = real(vars.npmu);
        gf_npmuim_(p.I) = imag(vars.npmu);
        gf_nplambdare_(p.I) = real(vars.nplambda);
        gf_nplambdaim_(p.I) = imag(vars.nplambda);
        gf_npnure_(p.I) = real(vars.npnu);
        gf_npnuim_(p.I) = imag(vars.npnu);
      });
}

} // namespace Weyl
