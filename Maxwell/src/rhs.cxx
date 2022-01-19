#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace Maxwell {
using namespace Loop;
using namespace std;

namespace {

template <int CI, int CJ, int CK, typename T>
T star(const GF3D<const T, CI, CJ, CK> &u, const vect<int, dim> &I) {
  constexpr int SI = CI == 0 ? +1 : -1;
  constexpr int SJ = CJ == 0 ? +1 : -1;
  constexpr int SK = CK == 0 ? +1 : -1;
  const auto DI = vect<int, dim>::unit(0);
  const auto DJ = vect<int, dim>::unit(1);
  const auto DK = vect<int, dim>::unit(2);
  const auto SDI = SI * DI;
  const auto SDJ = SJ * DJ;
  const auto SDK = SK * DK;
  return (u(I) + u(I + SDI) + u(I + SDJ) + u(I + SDI + SDJ) + u(I + SDK) +
          u(I + SDI + SDK) + u(I + SDJ + SDK) + u(I + SDI + SDJ + SDK)) /
         8;
}

} // namespace

////////////////////////////////////////////////////////////////////////////////

extern "C" void Maxwell_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Maxwell_RHS;
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

  const GF3D<CCTK_REAL, 0, 1, 1> dtdyz_(cctkGH, dtdyz);
  const GF3D<CCTK_REAL, 1, 0, 1> dtdzx_(cctkGH, dtdzx);
  const GF3D<CCTK_REAL, 1, 1, 0> dtdxy_(cctkGH, dtdxy);

  const GF3D<CCTK_REAL, 0, 1, 1> dtbyz_(cctkGH, dtbyz);
  const GF3D<CCTK_REAL, 1, 0, 1> dtbzx_(cctkGH, dtbzx);
  const GF3D<CCTK_REAL, 1, 1, 0> dtbxy_(cctkGH, dtbxy);

  loop_int<0, 1, 1>(cctkGH, [&](const PointDesc &p) {
    dtdyz_(p.I) = (star(bxy_, p.I + DJ) - star(bxy_, p.I)) / dy -
                  (star(bzx_, p.I + DK) - star(bzx_, p.I)) / dz;
  });
  loop_int<1, 0, 1>(cctkGH, [&](const PointDesc &p) {
    dtdzx_(p.I) = (star(byz_, p.I + DK) - star(byz_, p.I)) / dz -
                  (star(bxy_, p.I + DI) - star(bxy_, p.I)) / dx;
  });
  loop_int<1, 1, 0>(cctkGH, [&](const PointDesc &p) {
    dtdxy_(p.I) = (star(bzx_, p.I + DI) - star(bzx_, p.I)) / dx -
                  (star(byz_, p.I + DJ) - star(byz_, p.I)) / dy;
  });

  loop_int<0, 1, 1>(cctkGH, [&](const PointDesc &p) {
    dtbyz_(p.I) = (star(dzx_, p.I + DK) - star(dzx_, p.I)) / dz -
                  (star(dxy_, p.I + DJ) - star(dxy_, p.I)) / dy;
  });
  loop_int<1, 0, 1>(cctkGH, [&](const PointDesc &p) {
    dtbzx_(p.I) = (star(dxy_, p.I + DI) - star(dxy_, p.I)) / dx -
                  (star(dyz_, p.I + DK) - star(dyz_, p.I)) / dz;
  });
  loop_int<1, 1, 0>(cctkGH, [&](const PointDesc &p) {
    dtbxy_(p.I) = (star(dyz_, p.I + DJ) - star(dyz_, p.I)) / dy -
                  (star(dzx_, p.I + DI) - star(dzx_, p.I)) / dx;
  });
}

} // namespace Maxwell
