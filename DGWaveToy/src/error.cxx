#include <dg.hxx>
#include <loop.hxx>
#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace DGWaveToy {
using namespace DG;
using namespace std;

using Loop::dim;

extern "C" void DGWaveToy_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_DGWaveToy_Error;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> coordx_(cctkGH, coordx);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> coordy_(cctkGH, coordy);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> coordz_(cctkGH, coordz);

  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> u_(cctkGH, u);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> ft_(cctkGH, ft);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> fx_(cctkGH, fx);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> fy_(cctkGH, fy);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> fz_(cctkGH, fz);

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> u_error_(cctkGH, u_error);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> ft_error_(cctkGH, ft_error);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> fx_error_(cctkGH, fx_error);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> fy_error_(cctkGH, fy_error);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> fz_error_(cctkGH, fz_error);

  const CCTK_REAL kx = 2 * M_PI;
  const CCTK_REAL ky = 2 * M_PI;
  const CCTK_REAL kz = 2 * M_PI;
  const CCTK_REAL omega = sqrt(pow(kx, 2) + pow(ky, 2) + pow(kz, 2));

  Loop::loop_all<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL t = cctk_time;
    const CCTK_REAL x = coordx_(p.I);
    const CCTK_REAL y = coordy_(p.I);
    const CCTK_REAL z = coordz_(p.I);

    u_error_(p.I) =
        u_(p.I) - cos(omega * t) * cos(kx * x) * cos(ky * y) * cos(kz * z);

    ft_error_(p.I) = ft_(p.I) - -omega * sin(omega * t) * cos(kx * x) *
                                    cos(ky * y) * cos(kz * z);
    fx_error_(p.I) = fx_(p.I) - -kx * cos(omega * t) * sin(kx * x) *
                                    cos(ky * y) * cos(kz * z);
    fy_error_(p.I) = fy_(p.I) - -ky * cos(omega * t) * cos(kx * x) *
                                    sin(ky * y) * cos(kz * z);
    fz_error_(p.I) = fz_(p.I) - -kz * cos(omega * t) * cos(kx * x) *
                                    cos(ky * y) * sin(kz * z);
  });
}

} // namespace DGWaveToy
