#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace Maxwell {
using namespace Loop;
using namespace std;

extern "C" void Maxwell_Average(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Maxwell_Average;
  DECLARE_CCTK_PARAMETERS;

  const auto DI = vect<int, dim>::unit(0);
  const auto DJ = vect<int, dim>::unit(1);
  const auto DK = vect<int, dim>::unit(2);

  const GF3D<const CCTK_REAL, 0, 1, 1> dyz_(cctkGH, dyz);
  const GF3D<const CCTK_REAL, 1, 0, 1> dzx_(cctkGH, dzx);
  const GF3D<const CCTK_REAL, 1, 1, 0> dxy_(cctkGH, dxy);

  const GF3D<const CCTK_REAL, 0, 1, 1> byz_(cctkGH, byz);
  const GF3D<const CCTK_REAL, 1, 0, 1> bzx_(cctkGH, bzx);
  const GF3D<const CCTK_REAL, 1, 1, 0> bxy_(cctkGH, bxy);

  const GF3D<CCTK_REAL, 1, 1, 1> avgdyz_(cctkGH, avgdyz);
  const GF3D<CCTK_REAL, 1, 1, 1> avgdzx_(cctkGH, avgdzx);
  const GF3D<CCTK_REAL, 1, 1, 1> avgdxy_(cctkGH, avgdxy);

  const GF3D<CCTK_REAL, 1, 1, 1> avgbyz_(cctkGH, avgbyz);
  const GF3D<CCTK_REAL, 1, 1, 1> avgbzx_(cctkGH, avgbzx);
  const GF3D<CCTK_REAL, 1, 1, 1> avgbxy_(cctkGH, avgbxy);

  loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
    avgdyz_(p.I) = (dyz_(p.I) + dyz_(p.I + DI)) / 2;
  });
  loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
    avgdzx_(p.I) = (dzx_(p.I) + dzx_(p.I + DJ)) / 2;
  });
  loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
    avgdxy_(p.I) = (dxy_(p.I) + dxy_(p.I + DK)) / 2;
  });

  loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
    avgbyz_(p.I) = (byz_(p.I) + byz_(p.I + DI)) / 2;
  });
  loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
    avgbzx_(p.I) = (bzx_(p.I) + bzx_(p.I + DJ)) / 2;
  });
  loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
    avgbxy_(p.I) = (bxy_(p.I) + bxy_(p.I + DK)) / 2;
  });
}

} // namespace Maxwell
