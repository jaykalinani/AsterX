#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cmath>

namespace ErrorEstimator {
using namespace Loop;
using namespace std;

template <typename T> constexpr T pow2(const T x) { return x * x; }

template <typename T> constexpr T fmax3(const T x, const T y, const T z) {
  return fmax(fmax(x, y), z);
}

enum class shape_t { sphere, cube };

template <typename T>
constexpr T radius_sphere(const T x, const T y, const T z) {
  return sqrt(pow2(x) + pow2(y) + pow2(z));
}

template <typename T> constexpr T radius_cube(const T x, const T y, const T z) {
  return fmax3(fabs(x), fabs(y), fabs(z));
}

template <typename T, int CI, int CJ, int CK>
T lap(const GF3D<const T, CI, CJ, CK> &var, const vect<int, dim> &I) {
  const auto DI = vect<int, dim>::unit(0);
  const auto DJ = vect<int, dim>::unit(1);
  const auto DK = vect<int, dim>::unit(2);
  return fabs(var(I - DI) - 2 * var(I) + var(I + DI)) +
         fabs(var(I - DJ) - 2 * var(I) + var(I + DJ)) +
         fabs(var(I - DK) - 2 * var(I) + var(I + DK));
}

extern "C" void ErrorEstimator_Estimate(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ErrorEstimator_Estimate;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  const CCTK_REAL scalefactor = scale_by_resolution ? cbrt(dx * dy * dz) : 1;

  const GF3D<CCTK_REAL, 1, 1, 1> regrid_error_(cctkGH, regrid_error);

  typedef CCTK_REAL radius_t(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z);
  radius_t *const radius = [](const char *const sh) {
    if (CCTK_EQUALS(sh, "sphere"))
      return (radius_t *)radius_sphere;
    else if (CCTK_EQUALS(sh, "cube"))
      return (radius_t *)radius_cube;
    abort();
  }(region_shape);

  loop_int<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
    CCTK_REAL maxerr = 0;
    // Sample the neighbourhood of this point
    for (int k = -1; k <= +1; ++k) {
      for (int j = -1; j <= +1; ++j) {
        for (int i = -1; i <= +1; ++i) {
          const auto r =
              radius(p.x + i * dx / 2, p.y + j * dy / 2, p.z + k * dz / 2);
          const auto err = scalefactor / fmax(r, epsilon);
          maxerr = fmax(maxerr, err);
        }
      }
    }
    regrid_error_(p.I) = maxerr;
  });
}

} // namespace ErrorEstimator
