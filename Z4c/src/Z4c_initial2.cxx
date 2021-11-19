#include "derivs.hxx"
#include "physics.hxx"

#include <loop_device.hxx>
#include <simd.hxx>

#include <fixmath.hxx> // include this before <cctk.h>
#include <cctk.h>
#include <cctk_Arguments_Checked.h>

namespace Z4c {
using namespace Arith;
using namespace Loop;
using namespace std;

extern "C" void Z4c_Initial2(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_Initial2;

  const vec<CCTK_REAL, 3, UP> dx{
      CCTK_DELTA_SPACE(0),
      CCTK_DELTA_SPACE(1),
      CCTK_DELTA_SPACE(2),
  };

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout1(cctkGH, indextype);

  const smat<GF3D2<const CCTK_REAL>, 3, DN, DN> gf_gammat1{
      GF3D2<const CCTK_REAL>(layout1, gammatxx),
      GF3D2<const CCTK_REAL>(layout1, gammatxy),
      GF3D2<const CCTK_REAL>(layout1, gammatxz),
      GF3D2<const CCTK_REAL>(layout1, gammatyy),
      GF3D2<const CCTK_REAL>(layout1, gammatyz),
      GF3D2<const CCTK_REAL>(layout1, gammatzz)};

  const vec<GF3D2<CCTK_REAL>, 3, UP> gf_Gamt1{GF3D2<CCTK_REAL>(layout1, Gamtx),
                                              GF3D2<CCTK_REAL>(layout1, Gamty),
                                              GF3D2<CCTK_REAL>(layout1, Gamtz)};

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  const Loop::GridDescBaseDevice grid(cctkGH);
  grid.loop_int_device<0, 0, 0, vsize>(
      grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
        const GF3D2index index1(layout1, p.I);

        // Load
        const smat<vreal, 3, DN, DN> gammat =
            gf_gammat1(mask, index1, one<smat<int, 3, DN, DN> >()());

        // Calculate Z4c variables (only Gamt)
        const smat<vreal, 3, UP, UP> gammatu = calc_inv(gammat, vreal(1));

        const smat<vec<vreal, 3, DN>, 3, DN, DN> dgammat([&](int a, int b) {
          return deriv(mask, gf_gammat1(a, b), p.I, dx);
        });

        const vec<smat<vreal, 3, DN, DN>, 3, DN> Gammatl = calc_gammal(dgammat);
        const vec<smat<vreal, 3, DN, DN>, 3, UP> Gammat =
            calc_gamma(gammatu, Gammatl);
        const vec<vreal, 3, UP> Gamt([&](int a) ARITH_INLINE {
          return sum_symm<3>([&](int x, int y) ARITH_INLINE {
            return gammatu(x, y) * Gammat(a)(x, y);
          });
        });

        // Store
        gf_Gamt1.store(mask, index1, Gamt);
      });
}

} // namespace Z4c
