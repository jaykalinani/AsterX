#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

namespace Maxwell {

extern "C" void Maxwell_Constraints(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Maxwell_Constraints;
  DECLARE_CCTK_PARAMETERS;

  const auto DI = Loop::vect<int, Loop::dim>::unit(0);
  const auto DJ = Loop::vect<int, Loop::dim>::unit(1);
  const auto DK = Loop::vect<int, Loop::dim>::unit(2);

  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  const auto dxm{[&](const auto &u, const auto &p) {
    return (u(p.I) - u(p.I - DI)) / dx;
  }};
  const auto dym{[&](const auto &u, const auto &p) {
    return (u(p.I) - u(p.I - DJ)) / dy;
  }};
  const auto dzm{[&](const auto &u, const auto &p) {
    return (u(p.I) - u(p.I - DK)) / dz;
  }};

  const auto dxp{[&](const auto &u, const auto &p) {
    return (u(p.I + DI) - u(p.I)) / dx;
  }};
  const auto dyp{[&](const auto &u, const auto &p) {
    return (u(p.I + DJ) - u(p.I)) / dy;
  }};
  const auto dzp{[&](const auto &u, const auto &p) {
    return (u(p.I + DK) - u(p.I)) / dz;
  }};

  const Loop::GF3D<const CCTK_REAL, 1, 0, 0> ax_(cctkGH, ax);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 0> ay_(cctkGH, ay);
  const Loop::GF3D<const CCTK_REAL, 0, 0, 1> az_(cctkGH, az);

  const Loop::GF3D<const CCTK_REAL, 1, 0, 0> ex_(cctkGH, ex);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 0> ey_(cctkGH, ey);
  const Loop::GF3D<const CCTK_REAL, 0, 0, 1> ez_(cctkGH, ez);

  const Loop::GF3D<CCTK_REAL, 0, 1, 1> byz_(cctkGH, byz);
  const Loop::GF3D<CCTK_REAL, 1, 0, 1> bzx_(cctkGH, bzx);
  const Loop::GF3D<CCTK_REAL, 1, 1, 0> bxy_(cctkGH, bxy);

  const Loop::GF3D<CCTK_REAL, 0, 1, 1> curlayz_(cctkGH, curlayz);
  const Loop::GF3D<CCTK_REAL, 1, 0, 1> curlazx_(cctkGH, curlazx);
  const Loop::GF3D<CCTK_REAL, 1, 1, 0> curlaxy_(cctkGH, curlaxy);

  const Loop::GF3D<CCTK_REAL, 0, 0, 0> dive_(cctkGH, dive);

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> divb_(cctkGH, divb);

  Loop::loop_int<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    curlayz_(p.I) = byz_(p.I) - (dyp(az_, p) - dzp(ay_, p));
  });
  Loop::loop_int<1, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    curlazx_(p.I) = bzx_(p.I) - (dzp(ax_, p) - dxp(az_, p));
  });
  Loop::loop_int<1, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    curlaxy_(p.I) = bxy_(p.I) - (dxp(ay_, p) - dyp(ax_, p));
  });

  Loop::loop_int<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    dive_(p.I) = dxm(ex_, p) + dym(ey_, p) + dzm(ez_, p);
  });

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    divb_(p.I) = dxp(byz_, p) + dyp(bzx_, p) + dzp(bxy_, p);
  });
}

} // namespace Maxwell
