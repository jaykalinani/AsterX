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

  loop_all<0, 0, 0>(
      cctkGH, [&](const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Load
        const mat3<CCTK_REAL, DN, DN> gammat_old(
            gf_gammatxx_, gf_gammatxy_, gf_gammatxz_, gf_gammatyy_,
            gf_gammatyz_, gf_gammatzz_, p.I);
        const mat3<CCTK_REAL, DN, DN> At_old(gf_Atxx_, gf_Atxy_, gf_Atxz_,
                                             gf_Atyy_, gf_Atyz_, gf_Atzz_, p.I);

        // Enforce algebraic constraints
        // See arXiv:1212.2901 [gr-qc].

        const CCTK_REAL detgammat_old = gammat_old.det();
        const CCTK_REAL chi_old = 1 / cbrt(detgammat_old);
        const mat3<CCTK_REAL, DN, DN> gammat(
            [&](int a, int b) { return chi_old * gammat_old(a, b); });
#ifdef CCTK_DEBUG
        const CCTK_REAL detgammat = gammat.det();
        const CCTK_REAL gammat_norm = gammat.maxabs();
        const CCTK_REAL gammat_scale = gammat_norm;
        if (!(fabs(detgammat - 1) <= 1.0e-12 * gammat_scale)) {
          ostringstream buf;
          buf << "det gammat is not one: gammat=" << gammat
              << " det(gammat)=" << detgammat;
          CCTK_VERROR("%s", buf.str().c_str());
        }
        assert(fabs(detgammat - 1) <= 1.0e-12 * gammat_scale);
#endif

        const mat3<CCTK_REAL, UP, UP> gammatu = gammat.inv(1);

        const CCTK_REAL traceAt_old =
            sum2([&](int x, int y) { return gammatu(x, y) * At_old(x, y); });
        const mat3<CCTK_REAL, DN, DN> At([&](int a, int b) {
          return At_old(a, b) - traceAt_old / 3 * gammat(a, b);
        });
#ifdef CCTK_DEBUG
        const CCTK_REAL traceAt =
            sum2([&](int x, int y) { return gammatu(x, y) * At(x, y); });
        const CCTK_REAL gammatu_norm = gammatu.maxabs();
        const CCTK_REAL At_norm = At.maxabs();
        const CCTK_REAL At_scale =
            fmax(fmax(gammat_norm, gammatu_norm), At_norm);
        if (!(fabs(traceAt) <= 1.0e-12 * At_scale)) {
          ostringstream buf;
          buf << "tr At: At=" << At << " tr(At)=" << traceAt;
          CCTK_VERROR("%s", buf.str().c_str());
        }
        assert(fabs(traceAt) <= 1.0e-12 * At_scale);
#endif

        // Store
        gammat.store(gf_gammatxx_, gf_gammatxy_, gf_gammatxz_, gf_gammatyy_,
                     gf_gammatyz_, gf_gammatzz_, p.I);
        At.store(gf_Atxx_, gf_Atxy_, gf_Atxz_, gf_Atyy_, gf_Atyz_, gf_Atzz_,
                 p.I);
      });
}

} // namespace Z4c
