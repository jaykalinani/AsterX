#include <cctk.h>
#include <cctk_Arguments.h>
#include "reconstruct.hxx"

namespace ReconX {
namespace {
std::array<CCTK_REAL, 2> sample(const int method,
                               const std::array<CCTK_REAL, 6> &q) {
  switch (method) {
  case 0: return minmod_reconstruct(q[1], q[2], q[3], q[4]);
  case 1: return monocentral_reconstruct(q[1], q[2], q[3], q[4]);
  case 2: return wenoz_reconstruct(q[0], q[1], q[2], q[3], q[4], q[5], 1.0e-26);
  case 3: return wenozp_reconstruct(q[0], q[1], q[2], q[3], q[4], q[5],
                                    0.1, 1.0e-26, true);
  case 4: return mp5_reconstruct(q[0], q[1], q[2], q[3], q[4], q[5], 4.0);
  default: {
    reconstruct_params_t params{};
    params.ppm_shock_detection = false;
    params.ppm_zone_flattening = false;
    return ppm_reconstruct(q[0], q[1], q[2], q[3], q[4], q[5],
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, false, params);
  }
  }
}
void check(const CCTK_REAL actual, const CCTK_REAL expected) {
  if (!std::isfinite(actual) || fabs(actual-expected) > 1.0e-11*fmax(1.0, fabs(expected)))
    CCTK_VERROR("ReconX test failed: %.16e != %.16e", actual, expected);
}
} // namespace

extern "C" void ReconX_Test(CCTK_ARGUMENTS) {
  for (int method = 0; method < 6; ++method) {
    for (const CCTK_REAL value : {-0.03, 0.0, 0.5, 1.0e6}) {
      const auto constant = sample(method, {value, value, value, value, value, value});
      check(constant[0], value);
      check(constant[1], value);
    }
    const auto linear = sample(method, {-2.5, -1.5, -0.5, 0.5, 1.5, 2.5});
    check(linear[0], 0.0);
    check(linear[1], 0.0);
    const std::array<CCTK_REAL, 6> q{0.1, 0.2, 0.4, 0.45, 0.46, 0.5};
    const auto forward = sample(method, q);
    const auto reverse = sample(method, {q[5], q[4], q[3], q[2], q[1], q[0]});
    check(forward[0], reverse[1]);
    check(forward[1], reverse[0]);
    const auto jump = sample(method, {0.01, 0.01, 0.01, 1.0, 1.0, 1.0});
    if (!std::isfinite(jump[0]) || !std::isfinite(jump[1]))
      CCTK_ERROR("ReconX test: discontinuity produced a non-finite face");
    if (method < 2 && (jump[0] < 0.01 || jump[0] > 1.0 ||
                       jump[1] < 0.01 || jump[1] > 1.0))
      CCTK_ERROR("ReconX test: TVD reconstruction escaped the jump bounds");
  }
  CCTK_INFO("ReconX constant, linear, reflection and discontinuity tests passed");
}
} // namespace ReconX
