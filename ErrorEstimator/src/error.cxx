#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>

namespace ErrorEstimator {

namespace {
template <typename T> T pow2(T x) { return x * x; }

template <typename T>
T lap(const Loop::GF3D<const T, 1, 1, 1> &var, const Loop::PointDesc &p) {
  const auto DI = Loop::vect<int, Loop::dim>::unit(0);
  const auto DJ = Loop::vect<int, Loop::dim>::unit(1);
  const auto DK = Loop::vect<int, Loop::dim>::unit(2);
  return fabs(var(p.I - DI) - 2 * var(p.I) + var(p.I + DI)) +
         fabs(var(p.I - DJ) - 2 * var(p.I) + var(p.I + DJ)) +
         fabs(var(p.I - DK) - 2 * var(p.I) + var(p.I + DK));
}
} // namespace

extern "C" void ErrorEstimator_Estimate(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ErrorEstimator_Estimate;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> regrid_error_(cctkGH, regrid_error);

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    auto r = sqrt(pow2(p.x) + pow2(p.y) + pow2(p.z));
    regrid_error_(p.I) = 1 / fmax(r, CCTK_REAL(1.0e-10));
  });
}

} // namespace ErrorEstimator
