#include "physics.hxx"
#include "tensor.hxx"

#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cmath>
#include <sstream>

namespace Z4c {
using namespace Loop;
using namespace std;

extern "C" void Z4c_Enforce(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_Enforce;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<CCTK_REAL> gf_chi_(&layout, chi);

  const GF3D2<CCTK_REAL> gf_gammatxx_(&layout, gammatxx);
  const GF3D2<CCTK_REAL> gf_gammatxy_(&layout, gammatxy);
  const GF3D2<CCTK_REAL> gf_gammatxz_(&layout, gammatxz);
  const GF3D2<CCTK_REAL> gf_gammatyy_(&layout, gammatyy);
  const GF3D2<CCTK_REAL> gf_gammatyz_(&layout, gammatyz);
  const GF3D2<CCTK_REAL> gf_gammatzz_(&layout, gammatzz);

  const GF3D2<CCTK_REAL> gf_Atxx_(&layout, Atxx);
  const GF3D2<CCTK_REAL> gf_Atxy_(&layout, Atxy);
  const GF3D2<CCTK_REAL> gf_Atxz_(&layout, Atxz);
  const GF3D2<CCTK_REAL> gf_Atyy_(&layout, Atyy);
  const GF3D2<CCTK_REAL> gf_Atyz_(&layout, Atyz);
  const GF3D2<CCTK_REAL> gf_Atzz_(&layout, Atzz);

  const GF3D2<CCTK_REAL> gf_alphaG_(&layout, alphaG);

  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) Z4C_INLINE {
    // Load
    const CCTK_REAL chi_old = gf_chi_(p.I);
    const CCTK_REAL alphaG_old = gf_alphaG_(p.I);

    const mat3<CCTK_REAL, DN, DN> gammat_old(gf_gammatxx_, gf_gammatxy_,
                                             gf_gammatxz_, gf_gammatyy_,
                                             gf_gammatyz_, gf_gammatzz_, p.I);
    const mat3<CCTK_REAL, DN, DN> At_old(gf_Atxx_, gf_Atxy_, gf_Atxz_, gf_Atyy_,
                                         gf_Atyz_, gf_Atzz_, p.I);

    // Enforce floors

    const CCTK_REAL chi = fmax(chi_floor, chi_old);
    const CCTK_REAL alphaG = fmax(alphaG_floor, alphaG_old);

    // Enforce algebraic constraints
    // See arXiv:1212.2901 [gr-qc].

    const CCTK_REAL detgammat_old = gammat_old.det();
    const CCTK_REAL chi1_old = 1 / cbrt(detgammat_old);
    const mat3<CCTK_REAL, DN, DN> gammat(
        [&](int a, int b) { return chi1_old * gammat_old(a, b); });
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
    const CCTK_REAL At_scale = fmax(fmax(gammat_norm, gammatu_norm), At_norm);
    if (!(fabs(traceAt) <= 1.0e-12 * At_scale)) {
      ostringstream buf;
      buf << "tr At: At=" << At << " tr(At)=" << traceAt;
      CCTK_VERROR("%s", buf.str().c_str());
    }
    assert(fabs(traceAt) <= 1.0e-12 * At_scale);
#endif

    // Store
    gf_chi_(p.I) = chi;
    gammat.store(gf_gammatxx_, gf_gammatxy_, gf_gammatxz_, gf_gammatyy_,
                 gf_gammatyz_, gf_gammatzz_, p.I);
    At.store(gf_Atxx_, gf_Atxy_, gf_Atxz_, gf_Atyy_, gf_Atyz_, gf_Atzz_, p.I);
    gf_alphaG_(p.I) = alphaG;
  });
}

} // namespace Z4c
