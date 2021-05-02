#include "tensor.hxx"

#include <loop_device.hxx>
#include <simd.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

#include <cmath>

namespace Z4c {
using namespace Arith;
using namespace Loop;
using namespace std;

extern "C" void Z4c_Initial1(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_Initial1;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout1(cctkGH, indextype);

  const mat3<GF3D2<const CCTK_REAL>, DN, DN> gf_g1(
      GF3D2<const CCTK_REAL>(layout1, gxx),
      GF3D2<const CCTK_REAL>(layout1, gxy),
      GF3D2<const CCTK_REAL>(layout1, gxz),
      GF3D2<const CCTK_REAL>(layout1, gyy),
      GF3D2<const CCTK_REAL>(layout1, gyz),
      GF3D2<const CCTK_REAL>(layout1, gzz));

  const mat3<GF3D2<const CCTK_REAL>, DN, DN> gf_K1(
      GF3D2<const CCTK_REAL>(layout1, kxx),
      GF3D2<const CCTK_REAL>(layout1, kxy),
      GF3D2<const CCTK_REAL>(layout1, kxz),
      GF3D2<const CCTK_REAL>(layout1, kyy),
      GF3D2<const CCTK_REAL>(layout1, kyz),
      GF3D2<const CCTK_REAL>(layout1, kzz));

  const GF3D2<const CCTK_REAL> gf_alp1(layout1, alp);

  const vec3<GF3D2<const CCTK_REAL>, UP> gf_beta1(
      GF3D2<const CCTK_REAL>(layout1, betax),
      GF3D2<const CCTK_REAL>(layout1, betay),
      GF3D2<const CCTK_REAL>(layout1, betaz));

  const GF3D2<CCTK_REAL> gf_chi1(layout1, chi);

  const mat3<GF3D2<CCTK_REAL>, DN, DN> gf_gammat1(
      GF3D2<CCTK_REAL>(layout1, gammatxx), GF3D2<CCTK_REAL>(layout1, gammatxy),
      GF3D2<CCTK_REAL>(layout1, gammatxz), GF3D2<CCTK_REAL>(layout1, gammatyy),
      GF3D2<CCTK_REAL>(layout1, gammatyz), GF3D2<CCTK_REAL>(layout1, gammatzz));

  const GF3D2<CCTK_REAL> gf_Kh1(layout1, Kh);

  const mat3<GF3D2<CCTK_REAL>, DN, DN> gf_At1(
      GF3D2<CCTK_REAL>(layout1, Atxx), GF3D2<CCTK_REAL>(layout1, Atxy),
      GF3D2<CCTK_REAL>(layout1, Atxz), GF3D2<CCTK_REAL>(layout1, Atyy),
      GF3D2<CCTK_REAL>(layout1, Atyz), GF3D2<CCTK_REAL>(layout1, Atzz));

  const GF3D2<CCTK_REAL> gf_Theta1(layout1, Theta);

  const GF3D2<CCTK_REAL> gf_alphaG1(layout1, alphaG);

  const vec3<GF3D2<CCTK_REAL>, UP> gf_betaG1(GF3D2<CCTK_REAL>(layout1, betaGx),
                                             GF3D2<CCTK_REAL>(layout1, betaGy),
                                             GF3D2<CCTK_REAL>(layout1, betaGz));

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  const Loop::GridDescBaseDevice grid(cctkGH);
  grid.loop_all_device<0, 0, 0, vsize>(
      grid.nghostzones, [=](const PointDesc &p) Z4C_INLINE Z4C_GPU {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);

        // Load
        const mat3<vreal, DN, DN> g = gf_g1(mask, p.I, 1);
        const mat3<vreal, DN, DN> K = gf_K1(mask, p.I);
        const vreal alp = gf_alp1(mask, p.I, 1);
        const vec3<vreal, UP> beta = gf_beta1(mask, p.I);

        // Calculate Z4c variables (all except Gammat)
        const vreal detg = g.det();
        const mat3<vreal, UP, UP> gu = g.inv(detg);

        const vreal chi = 1 / cbrt(detg);

        const mat3<vreal, DN, DN> gammat(
            [&](int a, int b) Z4C_INLINE { return chi * g(a, b); });

        const vreal trK = sum2sym(
            [&](int x, int y) Z4C_INLINE { return gu(x, y) * K(x, y); });

        const vreal Theta = 0;

        const vreal Kh = trK - 2 * Theta;

        const mat3<vreal, DN, DN> At([&](int a, int b) Z4C_INLINE {
          return chi * (K(a, b) - trK / 3 * g(a, b));
        });

        const vreal alphaG = alp;

        const vec3<vreal, UP> betaG([&](int a) Z4C_INLINE { return beta(a); });

        // Store
        gf_chi1.store(mask, p.I, chi);
        gf_gammat1.store(mask, p.I, gammat);
        gf_Kh1.store(mask, p.I, Kh);
        gf_At1.store(mask, p.I, At);
        gf_Theta1.store(mask, p.I, Theta);
        gf_alphaG1.store(mask, p.I, alphaG);
        gf_betaG1.store(mask, p.I, betaG);
      });
}

} // namespace Z4c
