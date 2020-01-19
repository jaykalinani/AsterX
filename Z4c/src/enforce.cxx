#include "physics.hxx"
#include "tensor.hxx"

#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

#include <cmath>
#include <sstream>

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
                           gf_gammatyy_, gf_gammatyz_, gf_gammatzz_, p.I);
    mat3<CCTK_REAL> At(gf_Atxx_, gf_Atxy_, gf_Atxz_, gf_Atyy_, gf_Atyz_,
                       gf_Atzz_, p.I);

    // Enforce algebraic constraints
    // See arXiv:1212.2901 [gr-qc].

    const CCTK_REAL detgammat = gammat.det();
    const CCTK_REAL chi_new = 1 / cbrt(detgammat);
    for (int a = 0; a < 3; ++a)
      for (int b = a; b < 3; ++b)
        gammat(a, b) *= chi_new;
#ifdef CCTK_DEBUG
    const CCTK_REAL detgammat_new = gammat.det();
    const CCTK_REAL gammat_norm = gammat.maxabs();
    const CCTK_REAL gammat_scale = gammat_norm;
    if (!(fabs(detgammat_new - 1) <= 1.0e-12 * gammat_scale)) {
      ostringstream buf;
      buf << "det gammat is not one: gammat=" << gammat
          << " det(gammat)=" << detgammat_new;
      CCTK_VERROR("%s", buf.str().c_str());
    }
    assert(fabs(detgammat_new - 1) <= 1.0e-12 * gammat_scale);
#endif

    const mat3<CCTK_REAL> gammatu = gammat.inv(1);

    const CCTK_REAL traceAt =
        sum2([&](int x, int y) { return gammatu(x, y) * At(x, y); });
    for (int a = 0; a < 3; ++a)
      for (int b = a; b < 3; ++b)
        At(a, b) -= traceAt / 3 * gammat(a, b);
#ifdef CCTK_DEBUG
    const CCTK_REAL traceAt_new =
        sum2([&](int x, int y) { return gammatu(x, y) * At(x, y); });
    const CCTK_REAL gammatu_norm = gammatu.maxabs();
    const CCTK_REAL At_norm = At.maxabs();
    const CCTK_REAL At_scale = fmax(fmax(gammat_norm, gammatu_norm), At_norm);
    if (!(fabs(traceAt_new) <= 1.0e-12 * At_scale)) {
      ostringstream buf;
      buf << "tr At: At=" << At << " tr(At)=" << traceAt_new;
      CCTK_VERROR("%s", buf.str().c_str());
    }
    assert(fabs(traceAt_new) <= 1.0e-12 * At_scale);
#endif

    // Store
    gammat.store(gf_gammatxx_, gf_gammatxy_, gf_gammatxz_, gf_gammatyy_,
                 gf_gammatyz_, gf_gammatzz_, p);
    At.store(gf_Atxx_, gf_Atxy_, gf_Atxz_, gf_Atyy_, gf_Atyz_, gf_Atzz_, p);
  });
}

} // namespace Z4c
