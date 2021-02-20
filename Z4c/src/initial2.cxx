#include "derivs.hxx"
#include "field.hxx"
#include "physics.hxx"
#include "tensor.hxx"

#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

#include <cmath>

namespace Z4c {
using namespace Loop;
using namespace std;

extern "C" void Z4c_Initial2(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_Initial2;

  const vec3<CCTK_REAL, UP> dx{
      CCTK_DELTA_SPACE(0),
      CCTK_DELTA_SPACE(1),
      CCTK_DELTA_SPACE(2),
  };

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout1(cctkGH, indextype);

  const GF3D2<const CCTK_REAL> gf_gammatxx1(layout1, gammatxx);
  const GF3D2<const CCTK_REAL> gf_gammatxy1(layout1, gammatxy);
  const GF3D2<const CCTK_REAL> gf_gammatxz1(layout1, gammatxz);
  const GF3D2<const CCTK_REAL> gf_gammatyy1(layout1, gammatyy);
  const GF3D2<const CCTK_REAL> gf_gammatyz1(layout1, gammatyz);
  const GF3D2<const CCTK_REAL> gf_gammatzz1(layout1, gammatzz);

  const GF3D2<CCTK_REAL> gf_Gamtx1(layout1, Gamtx);
  const GF3D2<CCTK_REAL> gf_Gamty1(layout1, Gamty);
  const GF3D2<CCTK_REAL> gf_Gamtz1(layout1, Gamtz);

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    // Load
    const mat3<CCTK_REAL, DN, DN> gammat(gf_gammatxx1, gf_gammatxy1,
                                         gf_gammatxz1, gf_gammatyy1,
                                         gf_gammatyz1, gf_gammatzz1, p.I);

    // Calculate Z4c variables (only Gamt)
    const mat3<CCTK_REAL, UP, UP> gammatu = gammat.inv(1);

    const mat3<vec3<CCTK_REAL, DN>, DN, DN> dgammat{
        deriv(gf_gammatxx1, p.I, dx), deriv(gf_gammatxy1, p.I, dx),
        deriv(gf_gammatxz1, p.I, dx), deriv(gf_gammatyy1, p.I, dx),
        deriv(gf_gammatyz1, p.I, dx), deriv(gf_gammatzz1, p.I, dx),
    };

    const vec3<mat3<CCTK_REAL, DN, DN>, DN> Gammatl = calc_gammal(dgammat);
    const vec3<mat3<CCTK_REAL, DN, DN>, UP> Gammat =
        calc_gamma(gammatu, Gammatl);
    const vec3<CCTK_REAL, UP> Gamt([&](int a) {
      return sum2(
          [&](int x, int y) { return gammatu(x, y) * Gammat(a)(x, y); });
    });

    // Store
    Gamt.store(gf_Gamtx1, gf_Gamty1, gf_Gamtz1, p.I);
  });
}

} // namespace Z4c
