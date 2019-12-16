#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace Maxwell {

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

extern "C" void Maxwell_EstimateError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Maxwell_EstimateError;
  DECLARE_CCTK_PARAMETERS;

  // const Loop::GF3D<const CCTK_REAL, 1, 1, 1> avgphi_(cctkGH, avgphi);

  // const Loop::GF3D<const CCTK_REAL, 1, 1, 1> avgax_(cctkGH, avgax);
  // const Loop::GF3D<const CCTK_REAL, 1, 1, 1> avgay_(cctkGH, avgay);
  // const Loop::GF3D<const CCTK_REAL, 1, 1, 1> avgaz_(cctkGH, avgaz);

  // const Loop::GF3D<const CCTK_REAL, 1, 1, 1> avgex_(cctkGH, avgex);
  // const Loop::GF3D<const CCTK_REAL, 1, 1, 1> avgey_(cctkGH, avgey);
  // const Loop::GF3D<const CCTK_REAL, 1, 1, 1> avgez_(cctkGH, avgez);

  // const Loop::GF3D<const CCTK_REAL, 1, 1, 1> avgbyz_(cctkGH, avgbyz);
  // const Loop::GF3D<const CCTK_REAL, 1, 1, 1> avgbzx_(cctkGH, avgbzx);
  // const Loop::GF3D<const CCTK_REAL, 1, 1, 1> avgbxy_(cctkGH, avgbxy);

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> regrid_error_(cctkGH, regrid_error);

  if (false) {
    // Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    //   regrid_error_(p.I) = lap(avgphi_, p) + lap(avgax_, p) + lap(avgay_, p)
    //   +
    //                        lap(avgaz_, p) + lap(avgex_, p) + lap(avgey_, p) +
    //                        lap(avgez_, p) + lap(avgbyz_, p) + lap(avgbzx_, p)
    //                        + lap(avgbxy_, p);
    // });
  } else if (true) {
    Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      auto r = sqrt(pow2(p.x) + pow2(p.y) + pow2(p.z));
      regrid_error_(p.I) = 1 / fmax(r, CCTK_REAL(1.0e-10));
    });
  } else {
    assert(0);
  }
}

} // namespace Maxwell
