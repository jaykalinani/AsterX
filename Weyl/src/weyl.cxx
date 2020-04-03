#include "derivs.hxx"
#include "physics.hxx"
#include "tensor.hxx"
#include "weyl_vars.hxx"

#include <loop.hxx>

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

  const mat3<GF3D<const CCTK_REAL, 0, 0, 0>, DN, DN> gf_gamma_(
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, gxx),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, gxy),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, gxz),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, gyy),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, gyz),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, gzz));

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_alpha_(cctkGH, alp);

  const vec3<GF3D<const CCTK_REAL, 0, 0, 0>, UP> gf_beta_(
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, betax),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, betay),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, betaz));

  const mat3<GF3D<const CCTK_REAL, 0, 0, 0>, DN, DN> gf_k_(
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, kxx),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, kxy),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, kxz),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, kyy),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, kyz),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, kzz));

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_dtalpha_(cctkGH, dtalp);

  const vec3<GF3D<const CCTK_REAL, 0, 0, 0>, UP> gf_dtbeta_(
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, dtbetax),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, dtbetay),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, dtbetaz));

  const mat3<GF3D<const CCTK_REAL, 0, 0, 0>, DN, DN> gf_dtk_(
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, dtkxx),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, dtkxy),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, dtkxz),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, dtkyy),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, dtkyz),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, dtkzz));

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_dt2alpha_(cctkGH, dt2alp);

  const vec3<GF3D<const CCTK_REAL, 0, 0, 0>, UP> gf_dt2beta_(
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, dt2betax),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, dt2betay),
      GF3D<const CCTK_REAL, 0, 0, 0>(cctkGH, dt2betaz));

  //

  const mat3<vec3<GF3D<CCTK_REAL, 0, 0, 0>, DN>, DN, DN> gf_dgamma_(cctkGH,
                                                                    allocate());
  const mat3<mat3<GF3D<CCTK_REAL, 0, 0, 0>, DN, DN>, DN, DN> gf_ddgamma_(
      cctkGH, allocate());
  calc_derivs2(cctkGH, gf_gamma_, gf_dgamma_, gf_ddgamma_);

  const vec3<GF3D<CCTK_REAL, 0, 0, 0>, DN> gf_dalpha_(cctkGH, allocate());
  const mat3<GF3D<CCTK_REAL, 0, 0, 0>, DN, DN> gf_ddalpha_(cctkGH, allocate());
  calc_derivs2(cctkGH, gf_alpha_, gf_dalpha_, gf_ddalpha_);

  const vec3<vec3<GF3D<CCTK_REAL, 0, 0, 0>, DN>, UP> gf_dbeta_(cctkGH,
                                                               allocate());
  const vec3<mat3<GF3D<CCTK_REAL, 0, 0, 0>, DN, DN>, UP> gf_ddbeta_(cctkGH,
                                                                    allocate());
  calc_derivs2(cctkGH, gf_beta_, gf_dbeta_, gf_ddbeta_);

  const mat3<vec3<GF3D<CCTK_REAL, 0, 0, 0>, DN>, DN, DN> gf_dk_(cctkGH,
                                                                allocate());
  calc_derivs(cctkGH, gf_k_, gf_dk_);

  const vec3<GF3D<CCTK_REAL, 0, 0, 0>, DN> gf_ddtalpha_(cctkGH, allocate());
  calc_derivs(cctkGH, gf_dtalpha_, gf_ddtalpha_);

  const vec3<vec3<GF3D<CCTK_REAL, 0, 0, 0>, DN>, UP> gf_ddtbeta_(cctkGH,
                                                                 allocate());
  calc_derivs(cctkGH, gf_dtbeta_, gf_ddtbeta_);

  //

  const GF3D<CCTK_REAL, 0, 0, 0> gf_g4tt_(cctkGH, g4tt);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_g4tx_(cctkGH, g4tx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_g4ty_(cctkGH, g4ty);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_g4tz_(cctkGH, g4tz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_g4xx_(cctkGH, g4xx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_g4xy_(cctkGH, g4xy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_g4xz_(cctkGH, g4xz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_g4yy_(cctkGH, g4yy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_g4yz_(cctkGH, g4yz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_g4zz_(cctkGH, g4zz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4ttt_(cctkGH, Gamma4ttt);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4ttx_(cctkGH, Gamma4ttx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4tty_(cctkGH, Gamma4tty);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4ttz_(cctkGH, Gamma4ttz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4txx_(cctkGH, Gamma4txx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4txy_(cctkGH, Gamma4txy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4txz_(cctkGH, Gamma4txz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4tyy_(cctkGH, Gamma4tyy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4tyz_(cctkGH, Gamma4tyz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4tzz_(cctkGH, Gamma4tzz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4xtt_(cctkGH, Gamma4xtt);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4xtx_(cctkGH, Gamma4xtx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4xty_(cctkGH, Gamma4xty);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4xtz_(cctkGH, Gamma4xtz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4xxx_(cctkGH, Gamma4xxx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4xxy_(cctkGH, Gamma4xxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4xxz_(cctkGH, Gamma4xxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4xyy_(cctkGH, Gamma4xyy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4xyz_(cctkGH, Gamma4xyz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4xzz_(cctkGH, Gamma4xzz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4ytt_(cctkGH, Gamma4ytt);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4ytx_(cctkGH, Gamma4ytx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4yty_(cctkGH, Gamma4yty);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4ytz_(cctkGH, Gamma4ytz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4yxx_(cctkGH, Gamma4yxx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4yxy_(cctkGH, Gamma4yxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4yxz_(cctkGH, Gamma4yxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4yyy_(cctkGH, Gamma4yyy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4yyz_(cctkGH, Gamma4yyz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4yzz_(cctkGH, Gamma4yzz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4ztt_(cctkGH, Gamma4ztt);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4ztx_(cctkGH, Gamma4ztx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4zty_(cctkGH, Gamma4zty);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4ztz_(cctkGH, Gamma4ztz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4zxx_(cctkGH, Gamma4zxx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4zxy_(cctkGH, Gamma4zxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4zxz_(cctkGH, Gamma4zxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4zyy_(cctkGH, Gamma4zyy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4zyz_(cctkGH, Gamma4zyz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamma4zzz_(cctkGH, Gamma4zzz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_rm4txtx_(cctkGH, rm4txtx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_rm4txty_(cctkGH, rm4txty);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_rm4txtz_(cctkGH, rm4txtz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_rm4txxy_(cctkGH, rm4txxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_rm4txxz_(cctkGH, rm4txxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_rm4txyz_(cctkGH, rm4txyz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_rm4tyty_(cctkGH, rm4tyty);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_rm4tytz_(cctkGH, rm4tytz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_rm4tyxy_(cctkGH, rm4tyxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_rm4tyxz_(cctkGH, rm4tyxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_rm4tyyz_(cctkGH, rm4tyyz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_rm4tztz_(cctkGH, rm4tztz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_rm4tzxy_(cctkGH, rm4tzxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_rm4tzxz_(cctkGH, rm4tzxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_rm4tzyz_(cctkGH, rm4tzyz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_rm4xyxy_(cctkGH, rm4xyxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_rm4xyxz_(cctkGH, rm4xyxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_rm4xyyz_(cctkGH, rm4xyyz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_rm4xzxz_(cctkGH, rm4xzxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_rm4xzyz_(cctkGH, rm4xzyz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_rm4yzyz_(cctkGH, rm4yzyz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_r4tt_(cctkGH, r4tt);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_r4tx_(cctkGH, r4tx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_r4ty_(cctkGH, r4ty);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_r4tz_(cctkGH, r4tz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_r4xx_(cctkGH, r4xx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_r4xy_(cctkGH, r4xy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_r4xz_(cctkGH, r4xz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_r4yy_(cctkGH, r4yy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_r4yz_(cctkGH, r4yz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_r4zz_(cctkGH, r4zz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_rsc4_(cctkGH, rsc4);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_c4txtx_(cctkGH, c4txtx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_c4txty_(cctkGH, c4txty);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_c4txtz_(cctkGH, c4txtz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_c4txxy_(cctkGH, c4txxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_c4txxz_(cctkGH, c4txxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_c4txyz_(cctkGH, c4txyz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_c4tyty_(cctkGH, c4tyty);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_c4tytz_(cctkGH, c4tytz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_c4tyxy_(cctkGH, c4tyxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_c4tyxz_(cctkGH, c4tyxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_c4tyyz_(cctkGH, c4tyyz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_c4tztz_(cctkGH, c4tztz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_c4tzxy_(cctkGH, c4tzxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_c4tzxz_(cctkGH, c4tzxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_c4tzyz_(cctkGH, c4tzyz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_c4xyxy_(cctkGH, c4xyxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_c4xyxz_(cctkGH, c4xyxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_c4xyyz_(cctkGH, c4xyyz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_c4xzxz_(cctkGH, c4xzxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_c4xzyz_(cctkGH, c4xzyz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_c4yzyz_(cctkGH, c4yzyz);

  //

  loop_int<0, 0, 0>(
      cctkGH, [&](const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Load and calculate

        const weyl_vars<CCTK_REAL> vars(
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
      });
}

} // namespace Weyl
