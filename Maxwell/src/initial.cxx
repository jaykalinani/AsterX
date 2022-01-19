#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>

namespace Maxwell {
using namespace Loop;
using namespace std;

namespace {
template <typename T> T pow2(T x) { return x * x; }
} // namespace

////////////////////////////////////////////////////////////////////////////////

extern "C" void Maxwell_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Maxwell_Initial;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;

  const GF3D<CCTK_REAL, 0, 1, 1> dyz_(cctkGH, dyz);
  const GF3D<CCTK_REAL, 1, 0, 1> dzx_(cctkGH, dzx);
  const GF3D<CCTK_REAL, 1, 1, 0> dxy_(cctkGH, dxy);

  const GF3D<CCTK_REAL, 0, 1, 1> byz_(cctkGH, byz);
  const GF3D<CCTK_REAL, 1, 0, 1> bzx_(cctkGH, bzx);
  const GF3D<CCTK_REAL, 1, 1, 0> bxy_(cctkGH, bxy);

  if (CCTK_EQUALS(setup, "plane wave")) {

    // wave number
    const CCTK_REAL kx = CCTK_REAL(M_PI) * spatial_frequency_x;
    const CCTK_REAL ky = CCTK_REAL(M_PI) * spatial_frequency_y;
    const CCTK_REAL kz = CCTK_REAL(M_PI) * spatial_frequency_z;
    const CCTK_REAL omega = sqrt(pow2(kx) + pow2(ky) + pow2(kz));
    // amplitude
    const CCTK_REAL hx = amplitude_x;
    const CCTK_REAL hy = amplitude_y;
    const CCTK_REAL hz = amplitude_z;

    loop_int<0, 1, 1>(cctkGH, [&](const PointDesc &p) {
      dyz_(p.I) = hx * cos(omega * t - kx * p.x - ky * p.y - kz * p.z);
    });
    loop_int<1, 0, 1>(cctkGH, [&](const PointDesc &p) {
      dzx_(p.I) = hy * cos(omega * t - kx * p.x - ky * p.y - kz * p.z);
    });
    loop_int<1, 1, 0>(cctkGH, [&](const PointDesc &p) {
      dxy_(p.I) = hz * cos(omega * t - kx * p.x - ky * p.y - kz * p.z);
    });

    loop_int<0, 1, 1>(cctkGH, [&](const PointDesc &p) {
      byz_(p.I) = hz * cos(omega * t - kx * p.x - ky * p.y - kz * p.z);
    });
    loop_int<1, 0, 1>(cctkGH, [&](const PointDesc &p) {
      bzx_(p.I) = hx * cos(omega * t - kx * p.x - ky * p.y - kz * p.z);
    });
    loop_int<1, 1, 0>(cctkGH, [&](const PointDesc &p) {
      bxy_(p.I) = hy * cos(omega * t - kx * p.x - ky * p.y - kz * p.z);
    });

  } else if (CCTK_EQUALS(setup, "standing wave")) {

    // wave number
    const CCTK_REAL kx = CCTK_REAL(M_PI) * spatial_frequency_x;
    const CCTK_REAL ky = CCTK_REAL(M_PI) * spatial_frequency_y;
    const CCTK_REAL kz = CCTK_REAL(M_PI) * spatial_frequency_z;
    const CCTK_REAL omega = sqrt(pow2(kx) + pow2(ky) + pow2(kz));
    // amplitude
    const CCTK_REAL hx = amplitude_x;
    const CCTK_REAL hy = amplitude_y;
    const CCTK_REAL hz = amplitude_z;

    loop_int<0, 1, 1>(cctkGH, [&](const PointDesc &p) {
      dyz_(p.I) =
          hx * cos(omega * t) * sin(kx * p.x) * sin(ky * p.y) * sin(kz * p.z);
    });
    loop_int<1, 0, 1>(cctkGH, [&](const PointDesc &p) {
      dzx_(p.I) =
          hy * cos(omega * t) * sin(kx * p.x) * sin(ky * p.y) * sin(kz * p.z);
    });
    loop_int<1, 1, 0>(cctkGH, [&](const PointDesc &p) {
      dxy_(p.I) =
          hz * cos(omega * t) * sin(kx * p.x) * sin(ky * p.y) * sin(kz * p.z);
    });

    loop_int<0, 1, 1>(cctkGH, [&](const PointDesc &p) {
      byz_(p.I) =
          hz * cos(omega * t) * cos(kx * p.x) * cos(ky * p.y) * cos(kz * p.z);
    });
    loop_int<1, 0, 1>(cctkGH, [&](const PointDesc &p) {
      bzx_(p.I) =
          hx * cos(omega * t) * cos(kx * p.x) * cos(ky * p.y) * cos(kz * p.z);
    });
    loop_int<1, 1, 0>(cctkGH, [&](const PointDesc &p) {
      bxy_(p.I) =
          hy * cos(omega * t) * cos(kx * p.x) * cos(ky * p.y) * cos(kz * p.z);
    });

  } else {
    assert(0);
  }
}

} // namespace Maxwell
