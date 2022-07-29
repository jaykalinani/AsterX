#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "utils.hxx"
#include <array>

namespace AsterX
{
  using namespace Loop;

  template <int dir>
  void ComputeStaggeredB(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTSX_AsterX_ComputeBstagFromA;
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
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE
        {
          // Neighbouring "plus" and "minus" cell indices
          const auto ipjk = p.I + DI[0];
          const auto ijpk = p.I + DI[1];
          const auto ijkp = p.I + DI[2];

          if (dir == 0) {
            /* dBx is curl(A) at (i-1/2,j,k) */
            dBx_stag(p.I) = idx[1] * (Avec_z(ijpk) - Avec_z(p.I)) -
                            idx[2] * (Avec_y(ijkp) - Avec_y(p.I));
          } else if (dir == 1) {
            /* dBy is curl(A) at (i,j-1/2,k) */
            dBy_stag(p.I) = idx[2] * (Avec_x(ijkp) - Avec_x(p.I)) -
                            idx[0] * (Avec_z(ipjk) - Avec_z(p.I));
          } else if (dir == 2) {
            /* dBz is curl(A) at (i,j,z-1/2) */
            dBz_stag(p.I) = idx[0] * (Avec_y(ipjk) - Avec_y(p.I)) -
                            idx[1] * (Avec_x(ijpk) - Avec_x(p.I));
          }

          //TODO: need to implement copy conditions?
        });
  }

  extern "C" void AsterX_ComputeBstagFromA(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTSX_AsterX_ComputeBstagFromA;
    DECLARE_CCTK_PARAMETERS;

    ComputeStaggeredB<0>(cctkGH);
    ComputeStaggeredB<1>(cctkGH);
    ComputeStaggeredB<2>(cctkGH);
  }

  extern "C" void AsterX_ComputeBFromBstag(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTSX_AsterX_ComputeBFromBstag;
    DECLARE_CCTK_PARAMETERS;

    constexpr auto DI = PointDesc::DI;
    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE
        {
          // Neighbouring "plus" and "minus" cell indices
          const auto ipjk = p.I + DI[0];
          const auto ijpk = p.I + DI[1];
          const auto ijkp = p.I + DI[2];

          /* Second order interpolation of staggered B components to cell center */
          dBx(p.I) = 0.5 * (dBx_stag(p.I) + dBx_stag(ipjk));
          dBy(p.I) = 0.5 * (dBy_stag(p.I) + dBy_stag(ijpk));
          dBz(p.I) = 0.5 * (dBz_stag(p.I) + dBz_stag(ijkp));
        });
  }

} // namespace AsterX
