#include "physics.hxx"
#include "tensor.hxx"

#include <loop_device.hxx>
#include <simd.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cmath>
#include <sstream>

namespace Z4c {
using namespace Arith;
using namespace Loop;
using namespace std;

extern "C" void Z4c_Enforce(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_Enforce;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout1(cctkGH, indextype);

  const GF3D2<CCTK_REAL> gf_chi1(layout1, chi);

  const mat3<GF3D2<CCTK_REAL>, DN, DN> gf_gammat1(
      GF3D2<CCTK_REAL>(layout1, gammatxx), GF3D2<CCTK_REAL>(layout1, gammatxy),
      GF3D2<CCTK_REAL>(layout1, gammatxz), GF3D2<CCTK_REAL>(layout1, gammatyy),
      GF3D2<CCTK_REAL>(layout1, gammatyz), GF3D2<CCTK_REAL>(layout1, gammatzz));

  const mat3<GF3D2<CCTK_REAL>, DN, DN> gf_At1(
      GF3D2<CCTK_REAL>(layout1, Atxx), GF3D2<CCTK_REAL>(layout1, Atxy),
      GF3D2<CCTK_REAL>(layout1, Atxz), GF3D2<CCTK_REAL>(layout1, Atyy),
      GF3D2<CCTK_REAL>(layout1, Atyz), GF3D2<CCTK_REAL>(layout1, Atzz));

  const GF3D2<CCTK_REAL> gf_alphaG1(layout1, alphaG);

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  const Loop::GridDescBaseDevice grid(cctkGH);
  grid.loop_all_device<0, 0, 0, vsize>(
      grid.nghostzones, [=] Z4C_INLINE Z4C_GPU(const PointDesc &p) {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);

        // Load
        const vreal chi_old = gf_chi1(mask, p.I, 1);
        const vreal alphaG_old = gf_alphaG1(mask, p.I, 1);

        const mat3<vreal, DN, DN> gammat_old = gf_gammat1(mask, p.I, 1);
        const mat3<vreal, DN, DN> At_old = gf_At1(mask, p.I);

        // Enforce floors

        const vreal chi = fmax(vreal(chi_floor), chi_old);
        const vreal alphaG = fmax(vreal(alphaG_floor), alphaG_old);

        // Enforce algebraic constraints
        // See arXiv:1212.2901 [gr-qc].

        const vreal detgammat_old = gammat_old.det();
        const vreal chi1_old = 1 / cbrt(detgammat_old);
        const mat3<vreal, DN, DN> gammat([&] Z4C_INLINE(int a, int b) {
          return chi1_old * gammat_old(a, b);
        });
#ifdef CCTK_DEBUG
        const vreal detgammat = gammat.det();
        const vreal gammat_norm = gammat.maxabs();
        const vreal gammat_scale = gammat_norm;
        if (!(all(fabs(detgammat - 1) <= 1.0e-12 * gammat_scale))) {
          ostringstream buf;
          buf << "det gammat is not one: gammat=" << gammat
              << " det(gammat)=" << detgammat;
          CCTK_VERROR("%s", buf.str().c_str());
        }
        assert(all(fabs(detgammat - 1) <= 1.0e-12 * gammat_scale));
#endif

        const mat3<vreal, UP, UP> gammatu = gammat.inv(1);

        const vreal traceAt_old = sum2sym([&] Z4C_INLINE(int x, int y) {
          return gammatu(x, y) * At_old(x, y);
        });
        const mat3<vreal, DN, DN> At([&] Z4C_INLINE(int a, int b) {
          return At_old(a, b) - traceAt_old / 3 * gammat(a, b);
        });
#ifdef CCTK_DEBUG
        const vreal traceAt = sum2sym(
            [&] Z4C_INLINE(int x, int y) { return gammatu(x, y) * At(x, y); });
        const vreal gammatu_norm = gammatu.maxabs();
        const vreal At_norm = At.maxabs();
        const vreal At_scale = fmax(fmax(gammat_norm, gammatu_norm), At_norm);
        if (!(all(fabs(traceAt) <= 1.0e-12 * At_scale))) {
          ostringstream buf;
          buf << "tr At: At=" << At << " tr(At)=" << traceAt;
          CCTK_VERROR("%s", buf.str().c_str());
        }
        assert(all(fabs(traceAt) <= 1.0e-12 * At_scale));
#endif

        // Store
        gf_chi1.store(mask, p.I, chi);
        gf_gammat1.store(mask, p.I, gammat);
        gf_At1.store(mask, p.I, At);
        gf_alphaG1.store(mask, p.I, alphaG);
      });
}

} // namespace Z4c
