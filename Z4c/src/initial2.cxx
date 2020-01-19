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

  const vec3<CCTK_REAL> dx{
      CCTK_DELTA_SPACE(0),
      CCTK_DELTA_SPACE(1),
      CCTK_DELTA_SPACE(2),
  };

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatxx_(cctkGH, gammatxx);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatxy_(cctkGH, gammatxy);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatxz_(cctkGH, gammatxz);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatyy_(cctkGH, gammatyy);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatyz_(cctkGH, gammatyz);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatzz_(cctkGH, gammatzz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamtx_(cctkGH, Gamtx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamty_(cctkGH, Gamty);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamtz_(cctkGH, Gamtz);

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    // Load
    const mat3<CCTK_REAL> gammat(gf_gammatxx_, gf_gammatxy_, gf_gammatxz_,
                                 gf_gammatyy_, gf_gammatyz_, gf_gammatzz_, p.I);

    // Calculate Z4c variables (only Gamt)
    const mat3<CCTK_REAL> gammatu = gammat.inv(1);

    const mat3<vec3<CCTK_REAL> > dgammat{
        deriv(gf_gammatxx_, p.I, dx), deriv(gf_gammatxy_, p.I, dx),
        deriv(gf_gammatxz_, p.I, dx), deriv(gf_gammatyy_, p.I, dx),
        deriv(gf_gammatyz_, p.I, dx), deriv(gf_gammatzz_, p.I, dx),
    };

    const vec3<mat3<CCTK_REAL> > Gammatl = calc_gammal(dgammat);
    const vec3<mat3<CCTK_REAL> > Gammat = calc_gamma(gammatu, Gammatl);
    const vec3<CCTK_REAL> Gamt([&](int a) {
      return sum2(
          [&](int x, int y) { return gammatu(x, y) * Gammat(a)(x, y); });
    });

    // Store
    Gamt.store(gf_Gamtx_, gf_Gamty_, gf_Gamtz_, p);
  });
}

} // namespace Z4c
