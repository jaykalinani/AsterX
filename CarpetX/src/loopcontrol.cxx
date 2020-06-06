#include "cctk.h"

#include "loop.hxx"
#include "loopcontrol.h"

extern "C" GridDescBase_t LC_CreateGridDesc(const cGH *cctkGH) {
  Loop::GridDescBase grid_(cctkGH);
  GridDescBase_t grid;

  for (int d = 0; d < LC_DIM; ++d) {
    grid.gsh[d] = grid_.gsh[d];
    grid.lbnd[d] = grid_.lbnd[d];
    grid.ubnd[d] = grid_.ubnd[d];
    grid.lsh[d] = grid_.lsh[d];
    grid.ash[d] = grid_.ash[d];
    grid.bbox[2 * d + 0] = grid_.bbox[2 * d + 0];
    grid.bbox[2 * d + 1] = grid_.bbox[2 * d + 1];
    grid.nghostzones[d] = grid_.nghostzones[d];
    grid.tmin[d] = grid_.tmin[d];
    grid.tmax[d] = grid_.tmax[d];
    grid.x0[d] = grid_.x0[d];
    grid.dx[d] = grid_.dx[d];
  }

  return grid;
}
