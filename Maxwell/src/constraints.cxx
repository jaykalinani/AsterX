#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace Maxwell {
using namespace Loop;
using namespace std;

extern "C" void Maxwell_Constraints(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Maxwell_Constraints;
  DECLARE_CCTK_PARAMETERS;

  const auto DI = vect<int, dim>::unit(0);
  const auto DJ = vect<int, dim>::unit(1);
  const auto DK = vect<int, dim>::unit(2);

  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  const auto dxp = [&](const auto &u, const auto &I) {
    return (u(I + DI) - u(I)) / dx;
  };
  const auto dyp = [&](const auto &u, const auto &I) {
    return (u(I + DJ) - u(I)) / dy;
  };
  const auto dzp = [&](const auto &u, const auto &I) {
    return (u(I + DK) - u(I)) / dz;
  };

  const GF3D<const CCTK_REAL, 0, 1, 1> dyz_(cctkGH, dyz);
  const GF3D<const CCTK_REAL, 1, 0, 1> dzx_(cctkGH, dzx);
  const GF3D<const CCTK_REAL, 1, 1, 0> dxy_(cctkGH, dxy);

  const GF3D<const CCTK_REAL, 0, 1, 1> byz_(cctkGH, byz);
  const GF3D<const CCTK_REAL, 1, 0, 1> bzx_(cctkGH, bzx);
  const GF3D<const CCTK_REAL, 1, 1, 0> bxy_(cctkGH, bxy);

  const GF3D<CCTK_REAL, 1, 1, 1> divd_(cctkGH, divd);

  const GF3D<CCTK_REAL, 1, 1, 1> divb_(cctkGH, divb);

  loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
    divd_(p.I) = dxp(dyz_, p.I) + dyp(dzx_, p.I) + dzp(dxy_, p.I);
  });

  loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
    divb_(p.I) = dxp(byz_, p.I) + dyp(bzx_, p.I) + dzp(bxy_, p.I);
  });
}

} // namespace Maxwell
