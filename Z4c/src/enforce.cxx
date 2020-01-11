#include "physics.hxx"
#include "tensor.hxx"

#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

#include <cmath>

namespace Z4c {
using namespace Loop;
using namespace std;

extern "C" void Z4c_Enforce(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_Enforce;

  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatxx_(cctkGH, gammatxx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatxy_(cctkGH, gammatxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatxz_(cctkGH, gammatxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatyy_(cctkGH, gammatyy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatyz_(cctkGH, gammatyz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatzz_(cctkGH, gammatzz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atxx_(cctkGH, Atxx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atxy_(cctkGH, Atxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atxz_(cctkGH, Atxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atyy_(cctkGH, Atyy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atyz_(cctkGH, Atyz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atzz_(cctkGH, Atzz);

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    // Load
    mat3<CCTK_REAL> gammat(gf_gammatxx_, gf_gammatxy_, gf_gammatxz_,
                           gf_gammatyy_, gf_gammatyz_, gf_gammatzz_, p);
    mat3<CCTK_REAL> At(gf_Atxx_, gf_Atxy_, gf_Atxz_, gf_Atyy_, gf_Atyz_,
                       gf_Atzz_, p);

    // Enforce algebraic constraints
    // See arXiv:1212.2901 [gr-qc].

    const CCTK_REAL det_gammat = gammat.det();
    const CCTK_REAL psi1 = 1 / cbrt(det_gammat);
    for (int a = 0; a < 3; ++a)
      for (int b = a; b < 3; ++b)
        gammat(a, b) *= psi1;

    const mat3<CCTK_REAL> gammatu = gammat.inv(1);

    CCTK_REAL traceAt = 0;
    for (int x = 0; x < 3; ++x)
      for (int y = 0; y < 3; ++y)
        traceAt += gammatu(x, y) * At(x, y);

    for (int a = 0; a < 3; ++a)
      for (int b = a; b < 3; ++b)
        At(a, b) -= 1 / 3 * traceAt * gammat(a, b);

    // Store
    gammat.store(gf_gammatxx_, gf_gammatxy_, gf_gammatxz_, gf_gammatyy_,
                 gf_gammatyz_, gf_gammatzz_, p);
    At.store(gf_Atxx_, gf_Atxy_, gf_Atxz_, gf_Atyy_, gf_Atyz_, gf_Atzz_, p);
  });
}

} // namespace Z4c
