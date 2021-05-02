#include "tensor.hxx"
#include "z4c_vars.hxx"

#include <loop_device.hxx>
#include <simd.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace Z4c {
using namespace Arith;
using namespace Loop;
using namespace std;

extern "C" void Z4c_ADM(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_ADM;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout1(cctkGH, indextype);

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

  const mat3<GF3D2<CCTK_REAL>, DN, DN> gf_g1(
      GF3D2<CCTK_REAL>(layout1, gxx), GF3D2<CCTK_REAL>(layout1, gxy),
      GF3D2<CCTK_REAL>(layout1, gxz), GF3D2<CCTK_REAL>(layout1, gyy),
      GF3D2<CCTK_REAL>(layout1, gyz), GF3D2<CCTK_REAL>(layout1, gzz));

  const mat3<GF3D2<CCTK_REAL>, DN, DN> gf_K1(
      GF3D2<CCTK_REAL>(layout1, kxx), GF3D2<CCTK_REAL>(layout1, kxy),
      GF3D2<CCTK_REAL>(layout1, kxz), GF3D2<CCTK_REAL>(layout1, kyy),
      GF3D2<CCTK_REAL>(layout1, kyz), GF3D2<CCTK_REAL>(layout1, kzz));

  const GF3D2<CCTK_REAL> gf_alp1(layout1, alp);

  const GF3D2<CCTK_REAL> gf_dtalp1(layout1, dtalp);

  const vec3<GF3D2<CCTK_REAL>, UP> gf_beta1(GF3D2<CCTK_REAL>(layout1, betax),
                                            GF3D2<CCTK_REAL>(layout1, betay),
                                            GF3D2<CCTK_REAL>(layout1, betaz));

  const vec3<GF3D2<CCTK_REAL>, UP> gf_dtbeta1(
      GF3D2<CCTK_REAL>(layout1, dtbetax), GF3D2<CCTK_REAL>(layout1, dtbetay),
      GF3D2<CCTK_REAL>(layout1, dtbetaz));

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  const Loop::GridDescBaseDevice grid(cctkGH);
  grid.loop_all_device<0, 0, 0, vsize>(
      grid.nghostzones, [=](const PointDesc &p) Z4C_INLINE Z4C_GPU {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);

        // Load and calculate
        const z4c_vars_noderivs<vreal> vars(
            kappa1, kappa2, f_mu_L, f_mu_S,
            eta, //
            gf_chi1(mask, p.I, 1), gf_gammat1(mask, p.I, 1), gf_Kh1(mask, p.I),
            gf_At1(mask, p.I), gf_Gamt1(mask, p.I), gf_Theta1(mask, p.I),
            gf_alphaG1(mask, p.I, 1), gf_betaG1(mask, p.I), //
            gf_eTtt1(mask, p.I), gf_eTti1(mask, p.I), gf_eTij1(mask, p.I));

        // Store
        gf_g1.store(mask, p.I, vars.g);
        gf_K1.store(mask, p.I, vars.K);
        gf_alp1.store(mask, p.I, vars.alp);
        gf_dtalp1.store(mask, p.I, vars.dtalp);
        gf_beta1.store(mask, p.I, vars.beta);
        gf_dtbeta1.store(mask, p.I, vars.dtbeta);
      });
}

} // namespace Z4c
