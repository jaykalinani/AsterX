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
  const array<int, dim> nghostzones = {cctk_nghostzones[0], cctk_nghostzones[1],
                                       cctk_nghostzones[2]};

  const GF3D1<const CCTK_REAL> gf_gammatxx_(cctkGH, indextype, nghostzones,
                                            gammatxx);
  const GF3D1<const CCTK_REAL> gf_gammatxy_(cctkGH, indextype, nghostzones,
                                            gammatxy);
  const GF3D1<const CCTK_REAL> gf_gammatxz_(cctkGH, indextype, nghostzones,
                                            gammatxz);
  const GF3D1<const CCTK_REAL> gf_gammatyy_(cctkGH, indextype, nghostzones,
                                            gammatyy);
  const GF3D1<const CCTK_REAL> gf_gammatyz_(cctkGH, indextype, nghostzones,
                                            gammatyz);
  const GF3D1<const CCTK_REAL> gf_gammatzz_(cctkGH, indextype, nghostzones,
                                            gammatzz);

  const GF3D1<CCTK_REAL> gf_Gamtx_(cctkGH, indextype, nghostzones, Gamtx);
  const GF3D1<CCTK_REAL> gf_Gamty_(cctkGH, indextype, nghostzones, Gamty);
  const GF3D1<CCTK_REAL> gf_Gamtz_(cctkGH, indextype, nghostzones, Gamtz);

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    // Load
    const mat3<CCTK_REAL, DN, DN> gammat(gf_gammatxx_, gf_gammatxy_,
                                         gf_gammatxz_, gf_gammatyy_,
                                         gf_gammatyz_, gf_gammatzz_, p.I);

    // Calculate Z4c variables (only Gamt)
    const mat3<CCTK_REAL, UP, UP> gammatu = gammat.inv(1);

    const mat3<vec3<CCTK_REAL, DN>, DN, DN> dgammat{
        deriv(gf_gammatxx_, p.I, dx), deriv(gf_gammatxy_, p.I, dx),
        deriv(gf_gammatxz_, p.I, dx), deriv(gf_gammatyy_, p.I, dx),
        deriv(gf_gammatyz_, p.I, dx), deriv(gf_gammatzz_, p.I, dx),
    };

    const vec3<mat3<CCTK_REAL, DN, DN>, DN> Gammatl = calc_gammal(dgammat);
    const vec3<mat3<CCTK_REAL, DN, DN>, UP> Gammat =
        calc_gamma(gammatu, Gammatl);
    const vec3<CCTK_REAL, UP> Gamt([&](int a) {
      return sum2(
          [&](int x, int y) { return gammatu(x, y) * Gammat(a)(x, y); });
    });

    // Store
    Gamt.store(gf_Gamtx_, gf_Gamty_, gf_Gamtz_, p.I);
  });
}

} // namespace Z4c
