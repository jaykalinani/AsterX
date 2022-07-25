#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "utils.hxx"
#include <array>

namespace AsterX {
using namespace Loop;

template <int dir> void ComputeStaggeredB(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_ComputeBfromA;
  DECLARE_CCTK_PARAMETERS;

  const std::array dx{CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1),
                      CCTK_DELTA_SPACE(2)};
  const std::array idx{1 / dx[0], 1 / dx[1], 1 / dx[2]};

  constexpr auto DI = PointDesc::DI;

  static_assert(dir >= 0 && dir < 3, "");

  constexpr array<int, dim> face_centred = {!(dir == 0), !(dir == 1),
                                            !(dir == 2)};
  grid.loop_int_device<face_centred[0], face_centred[1], face_centred[2]>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

        // Neighbouring "plus" and "minus" cell indices
        const auto imjk = p.I - DI[0];
        const auto ijmk = p.I - DI[1];
        const auto ijkm = p.I - DI[2];

        /* dBx is curl(A) at i+1/2 */
        dBx_stag(p.I) = idx[1]*(Avec_z(p.I) - Avec_z(ijmk)) -
                        idx[2]*(Avec_y(p.I) - Avec_y(ijkm));

        /* dBy is curl(A) at j+1/2 */
        dBy_stag(p.I) = idx[2]*(Avec_x(p.I) - Avec_x(ijkm)) -
                        idx[0]*(Avec_z(p.I) - Avec_z(imjk));

        /* dBy is curl(A) at j+1/2 */
        dBy_stag(p.I) = idx[0]*(Avec_y(p.I) - Avec_y(imjk)) -
                        idx[1]*(Avec_x(p.I) - Avec_x(ijmk));

        //TODO: need to implement copy conditions?
      });
}

extern "C" void AsterX_ComputeBfromA(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_ComputeBfromA;
  DECLARE_CCTK_PARAMETERS;

  ComputeStaggeredB<0>(cctkGH);
  ComputeStaggeredB<1>(cctkGH);
  ComputeStaggeredB<2>(cctkGH);

  constexpr auto DI = PointDesc::DI;
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

        // Neighbouring "plus" and "minus" cell indices
        const auto imjk = p.I - DI[0];
        const auto ijmk = p.I - DI[1];
        const auto ijkm = p.I - DI[2];

        dBx(p.I) = 0.5 * (dBx_stag(p.I) + dBx_stag(imjk));
        dBy(p.I) = 0.5 * (dBy_stag(p.I) + dBy_stag(ijmk));
        dBz(p.I) = 0.5 * (dBz_stag(p.I) + dBz_stag(ijkm));
      });
}

} // namespace AsterX
