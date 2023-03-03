#include "dg.hxx"

#include <div.hxx>
#include <loop.hxx>
#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cstdlib>

namespace DGCoordinates {
using namespace Arith;
using namespace DG;
using namespace std;

using Loop::dim;

extern "C" void DGCoordinates_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_DGCoordinates_Setup;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> coordx_(cctkGH, coordx);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> coordy_(cctkGH, coordy);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> coordz_(cctkGH, coordz);

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> cvol_(cctkGH, cvol);

  const Loop::GridDescBase grid(cctkGH);

  for (int d = 0; d < dim; ++d)
    assert(grid.nghostzones[d] == 1);
  for (int d = 0; d < dim; ++d)
    assert(mod_floor(grid.gsh[d] - 2 * grid.nghostzones[d], dg_npoints) == 0);
  for (int d = 0; d < dim; ++d)
    assert(mod_floor(grid.lbnd[d], dg_npoints) == 0);
  for (int d = 0; d < dim; ++d)
    assert(mod_floor(grid.lsh[d] - 2 * grid.nghostzones[d], dg_npoints) == 0);

  switch (dg_npoints) {

  case 2: {
    constexpr int N = 2;
    grid.loop_all<1, 1, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
      const get_coords<CCTK_REAL, N> coords(grid, p);
      coordx_(p.I) = coords.coord[0];
      coordy_(p.I) = coords.coord[1];
      coordz_(p.I) = coords.coord[2];
      cvol_(p.I) = coords.vol;
    });
    break;
  }

  case 3: {
    constexpr int N = 3;
    grid.loop_all<1, 1, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
      const get_coords<CCTK_REAL, N> coords(grid, p);
      coordx_(p.I) = coords.coord[0];
      coordy_(p.I) = coords.coord[1];
      coordz_(p.I) = coords.coord[2];
      cvol_(p.I) = coords.vol;
    });
    break;
  }

  case 4: {
    constexpr int N = 4;
    grid.loop_all<1, 1, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
      const get_coords<CCTK_REAL, N> coords(grid, p);
      coordx_(p.I) = coords.coord[0];
      coordy_(p.I) = coords.coord[1];
      coordz_(p.I) = coords.coord[2];
      cvol_(p.I) = coords.vol;
    });
    break;
  }

  case 8: {
    constexpr int N = 8;
    grid.loop_all<1, 1, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
      const get_coords<CCTK_REAL, N> coords(grid, p);
      coordx_(p.I) = coords.coord[0];
      coordy_(p.I) = coords.coord[1];
      coordz_(p.I) = coords.coord[2];
      cvol_(p.I) = coords.vol;
    });
    break;
  }

  default:
    assert(0);
  }
}

} // namespace DGCoordinates
