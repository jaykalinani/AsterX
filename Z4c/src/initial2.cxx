#include "derivs.hxx"
#include "field.hxx"
#include "physics.hxx"
#include "tensor.hxx"

#include <loop_device.hxx>
#include <simd.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

namespace Z4c {
using namespace Arith;
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

  const mat3<GF3D2<const CCTK_REAL>, DN, DN> gf_gammat1(
      GF3D2<const CCTK_REAL>(layout1, gammatxx),
      GF3D2<const CCTK_REAL>(layout1, gammatxy),
      GF3D2<const CCTK_REAL>(layout1, gammatxz),
      GF3D2<const CCTK_REAL>(layout1, gammatyy),
      GF3D2<const CCTK_REAL>(layout1, gammatyz),
      GF3D2<const CCTK_REAL>(layout1, gammatzz));

  const vec3<GF3D2<CCTK_REAL>, UP> gf_Gamt1(GF3D2<CCTK_REAL>(layout1, Gamtx),
                                            GF3D2<CCTK_REAL>(layout1, Gamty),
                                            GF3D2<CCTK_REAL>(layout1, Gamtz));

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  const Loop::GridDescBaseDevice grid(cctkGH);
  grid.loop_int_device<0, 0, 0, vsize>(
      grid.nghostzones, [=](const PointDesc &p) Z4C_INLINE Z4C_GPU {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
        const GF3D2index index1(layout1, p.I);

        // Load
        const mat3<vreal, DN, DN> gammat = gf_gammat1(mask, index1, 1);

        // Calculate Z4c variables (only Gamt)
        const mat3<vreal, UP, UP> gammatu = gammat.inv(1);

        const mat3<vec3<vreal, DN>, DN, DN> dgammat([&](int a, int b) {
          return deriv(mask, gf_gammat1(a, b), p.I, dx);
        });

        const vec3<mat3<vreal, DN, DN>, DN> Gammatl = calc_gammal(dgammat);
        const vec3<mat3<vreal, DN, DN>, UP> Gammat =
            calc_gamma(gammatu, Gammatl);
        const vec3<vreal, UP> Gamt([&](int a) Z4C_INLINE {
          return sum2sym([&](int x, int y) Z4C_INLINE {
            return gammatu(x, y) * Gammat(a)(x, y);
          });
        });

        // Store
        gf_Gamt1.store(mask, index1, Gamt);
      });
}

} // namespace Z4c
