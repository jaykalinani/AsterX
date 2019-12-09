#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

namespace Maxwell {

extern "C" void Maxwell_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Maxwell_RHS;
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

  const Loop::GF3D<const CCTK_REAL, 0, 0, 0> phi_(cctkGH, phi);

  const Loop::GF3D<const CCTK_REAL, 1, 0, 0> ax_(cctkGH, ax);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 0> ay_(cctkGH, ay);
  const Loop::GF3D<const CCTK_REAL, 0, 0, 1> az_(cctkGH, az);

  const Loop::GF3D<const CCTK_REAL, 1, 0, 0> ex_(cctkGH, ex);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 0> ey_(cctkGH, ey);
  const Loop::GF3D<const CCTK_REAL, 0, 0, 1> ez_(cctkGH, ez);

  const Loop::GF3D<const CCTK_REAL, 0, 1, 1> byz_(cctkGH, byz);
  const Loop::GF3D<const CCTK_REAL, 1, 0, 1> bzx_(cctkGH, bzx);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 0> bxy_(cctkGH, bxy);

  const Loop::GF3D<CCTK_REAL, 0, 0, 0> dtphi_(cctkGH, dtphi);

  const Loop::GF3D<CCTK_REAL, 1, 0, 0> dtax_(cctkGH, dtax);
  const Loop::GF3D<CCTK_REAL, 0, 1, 0> dtay_(cctkGH, dtay);
  const Loop::GF3D<CCTK_REAL, 0, 0, 1> dtaz_(cctkGH, dtaz);

  const Loop::GF3D<CCTK_REAL, 1, 0, 0> dtex_(cctkGH, dtex);
  const Loop::GF3D<CCTK_REAL, 0, 1, 0> dtey_(cctkGH, dtey);
  const Loop::GF3D<CCTK_REAL, 0, 0, 1> dtez_(cctkGH, dtez);

  const Loop::GF3D<CCTK_REAL, 0, 1, 1> dtbyz_(cctkGH, dtbyz);
  const Loop::GF3D<CCTK_REAL, 1, 0, 1> dtbzx_(cctkGH, dtbzx);
  const Loop::GF3D<CCTK_REAL, 1, 1, 0> dtbxy_(cctkGH, dtbxy);

  Loop::loop_int<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    dtphi_(p.I) = -(dxm(ax_, p) + dym(ay_, p) + dzm(az_, p));
  });

  Loop::loop_int<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    dtax_(p.I) = -dxp(phi_, p) - ex_(p.I);
  });
  Loop::loop_int<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    dtay_(p.I) = -dyp(phi_, p) - ey_(p.I);
  });
  Loop::loop_int<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    dtaz_(p.I) = -dzp(phi_, p) - ez_(p.I);
  });

  Loop::loop_int<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    dtex_(p.I) = -(dzm(bzx_, p) - dym(bxy_, p));
  });
  Loop::loop_int<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    dtey_(p.I) = -(dxm(bxy_, p) - dzm(byz_, p));
  });
  Loop::loop_int<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    dtez_(p.I) = -(dym(byz_, p) - dxm(bzx_, p));
  });

  Loop::loop_int<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    dtbyz_(p.i, p.j, p.k) = dzp(ey_, p) - dyp(ez_, p);
  });
  Loop::loop_int<1, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    dtbzx_(p.i, p.j, p.k) = dxp(ez_, p) - dzp(ex_, p);
  });
  Loop::loop_int<1, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    dtbxy_(p.i, p.j, p.k) = dyp(ex_, p) - dxp(ey_, p);
  });
}

} // namespace Maxwell
