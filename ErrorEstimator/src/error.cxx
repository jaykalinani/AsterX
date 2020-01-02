#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>

namespace ErrorEstimator {
using namespace Loop;
using namespace std;

namespace {
template <typename T> T pow2(T x) { return x * x; }

template <typename T, int CI, int CJ, int CK>
T lap(const GF3D<const T, CI, CJ, CK> &var, const vect<int, dim> &I) {
  const auto DI = vect<int, dim>::unit(0);
  const auto DJ = vect<int, dim>::unit(1);
  const auto DK = vect<int, dim>::unit(2);
  return fabs(var(I - DI) - 2 * var(I) + var(I + DI)) +
         fabs(var(I - DJ) - 2 * var(I) + var(I + DJ)) +
         fabs(var(I - DK) - 2 * var(I) + var(I + DK));
}
} // namespace

extern "C" void ErrorEstimator_Estimate(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ErrorEstimator_Estimate;
  DECLARE_CCTK_PARAMETERS;

  const GF3D<CCTK_REAL, 1, 1, 1> regrid_error_(cctkGH, regrid_error);

  const CCTK_REAL scalefactor = scale_by_resolution ?
    pow(CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2), 1./3.) :
    1.;

  loop_int<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
    auto r = sqrt(pow2(p.x) + pow2(p.y) + pow2(p.z));
    regrid_error_(p.I) = scalefactor / fmax(r, CCTK_REAL(1.0e-10));
  });
}

} // namespace ErrorEstimator
