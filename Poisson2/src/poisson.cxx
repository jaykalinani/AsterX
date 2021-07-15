#include <pdesolvers.hxx>

#include <defs.hxx>
#include <loop.hxx>

#include <fixmath.hxx> // include this before <cctk.h>
#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <array>
#include <cmath>
#include <cstdlib>

namespace Poisson2 {

extern "C" void Poisson2_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Poisson2_Init;

  const int dim = Loop::dim;

  const std::array<int, dim> indextype = {0, 0, 0};
  const std::array<int, dim> nghostzones = {
      cctk_nghostzones[0], cctk_nghostzones[1], cctk_nghostzones[2]};
  Arith::vect<int, dim> imin, imax;
  Loop::GridDescBase(cctkGH).box_int<0, 0, 0>(nghostzones, imin, imax);
  const Loop::GF3D2layout layout1(cctkGH, indextype);

  const Loop::GF3D2<CCTK_REAL> gf_sol(layout1, sol);

  const Loop::GridDescBase grid(cctkGH);
  grid.loop_int<0, 0, 0>(grid.nghostzones,
                         [=](const Loop::PointDesc &p) ARITH_INLINE {
                           const Loop::GF3D2index index1(layout1, p.I);
                           gf_sol(index1) = 0;
                         });
}

extern "C" void Poisson2_Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Poisson2_Boundaries;

  const int dim = Loop::dim;

  const std::array<int, dim> indextype = {0, 0, 0};
  const std::array<int, dim> nghostzones = {
      cctk_nghostzones[0], cctk_nghostzones[1], cctk_nghostzones[2]};
  Arith::vect<int, dim> imin, imax;
  Loop::GridDescBase(cctkGH).box_int<0, 0, 0>(nghostzones, imin, imax);
  const Loop::GF3D2layout layout1(cctkGH, indextype);

  const Loop::GF3D2<CCTK_REAL> gf_sol(layout1, sol);

  const Loop::GridDescBase grid(cctkGH);
  grid.loop_bnd<0, 0, 0>(grid.nghostzones,
                         [=](const Loop::PointDesc &p) ARITH_INLINE {
                           const Loop::GF3D2index index1(layout1, p.I);
                           gf_sol(index1) = 0;
                         });
}

extern "C" void Poisson2_Residual(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Poisson2_Residual;
  DECLARE_CCTK_PARAMETERS;

  const int dim = Loop::dim;

  const std::array<int, dim> indextype = {0, 0, 0};
  const std::array<int, dim> nghostzones = {
      cctk_nghostzones[0], cctk_nghostzones[1], cctk_nghostzones[2]};
  Arith::vect<int, dim> imin, imax;
  Loop::GridDescBase(cctkGH).box_int<0, 0, 0>(nghostzones, imin, imax);
  const Loop::GF3D2layout layout1(cctkGH, indextype);

  for (int d = 0; d < dim; ++d)
    assert(2 * nghostzones[d] >= fdorder);

  const Loop::GF3D2<const CCTK_REAL> gf_sol(layout1, sol);
  const Loop::GF3D2<CCTK_REAL> gf_res(layout1, res);

  const Loop::GridDescBase grid(cctkGH);
  grid.loop_int<0, 0, 0>(
      grid.nghostzones, [=](const Loop::PointDesc &p) ARITH_INLINE {
        const Loop::GF3D2index index1(layout1, p.I);
        using Arith::pow2;
        CCTK_REAL ddsol = 0;
        switch (fdorder) {
        case 2:
          for (int d = 0; d < dim; ++d)
            ddsol += (gf_sol(Loop::GF3D2index(layout1, p.I - p.DI[d])) //
                      - 2 * gf_sol(Loop::GF3D2index(layout1, p.I))     //
                      + gf_sol(Loop::GF3D2index(layout1, p.I + p.DI[d]))) /
                     pow2(p.DX[d]);
          break;
        case 4:
          for (int d = 0; d < dim; ++d)
            ddsol +=
                (-1 / 12.0 *
                     gf_sol(Loop::GF3D2index(layout1, p.I - 2 * p.DI[d]))     //
                 + 4 / 3.0 * gf_sol(Loop::GF3D2index(layout1, p.I - p.DI[d])) //
                 - 5 / 2.0 * gf_sol(Loop::GF3D2index(layout1, p.I))           //
                 + 4 / 3.0 * gf_sol(Loop::GF3D2index(layout1, p.I + p.DI[d])) //
                 - 1 / 12.0 *
                       gf_sol(Loop::GF3D2index(layout1, p.I + 2 * p.DI[d]))) /
                pow2(p.DX[d]);
          break;
        default:
          abort();
        }
        gf_res(index1) = ddsol - 1;
      });
}

extern "C" void Poisson2_Jacobian(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Poisson2_Jacobian;
  DECLARE_CCTK_PARAMETERS;

  const int dim = Loop::dim;

  const std::array<int, dim> indextype = {0, 0, 0};
  const std::array<int, dim> nghostzones = {
      cctk_nghostzones[0], cctk_nghostzones[1], cctk_nghostzones[2]};
  Arith::vect<int, dim> imin, imax;
  Loop::GridDescBase(cctkGH).box_int<0, 0, 0>(nghostzones, imin, imax);
  const Loop::GF3D2layout layout1(cctkGH, indextype);

  const Loop::GF3D2<const CCTK_REAL> gf_idx(layout1, idx);
  const Loop::GF3D2<const CCTK_REAL> gf_sol(layout1, sol);

  assert(PDESolvers::jacobian.has_value());
  const PDESolvers::jacobian_t &J = PDESolvers::jacobian.value();

  const Loop::GridDescBase grid(cctkGH);
  grid.loop_int<0, 0, 0>(
      grid.nghostzones, [=](const Loop::PointDesc &p) ARITH_INLINE {
        const Loop::GF3D2index index1(layout1, p.I);
        using std::lrint;
        using Arith::pow2;
        const int idx = lrint(gf_idx(index1));
        switch (fdorder) {
        case 2:
          for (int d = 0; d < dim; ++d) {
            J.add_value(idx,
                        lrint(gf_idx(Loop::GF3D2index(layout1, p.I - p.DI[d]))),
                        1 / pow2(p.DX[d]));
            J.add_value(idx, lrint(gf_idx(Loop::GF3D2index(layout1, p.I))),
                        -2 / pow2(p.DX[d]));
            J.add_value(idx,
                        lrint(gf_idx(Loop::GF3D2index(layout1, p.I + p.DI[d]))),
                        1 / pow2(p.DX[d]));
          }
          break;
        case 4:
          for (int d = 0; d < dim; ++d) {
            J.add_value(
                idx,
                lrint(gf_idx(Loop::GF3D2index(layout1, p.I - 2 * p.DI[d]))),
                -1 / 12.0 / pow2(p.DX[d]));
            J.add_value(idx,
                        lrint(gf_idx(Loop::GF3D2index(layout1, p.I - p.DI[d]))),
                        4 / 3.0 / pow2(p.DX[d]));
            J.add_value(idx, lrint(gf_idx(Loop::GF3D2index(layout1, p.I))),
                        -5 / 2.0 / pow2(p.DX[d]));
            J.add_value(idx,
                        lrint(gf_idx(Loop::GF3D2index(layout1, p.I + p.DI[d]))),
                        4 / 3.0 / pow2(p.DX[d]));
            J.add_value(
                idx,
                lrint(gf_idx(Loop::GF3D2index(layout1, p.I + 2 * p.DI[d]))),
                -1 / 12.0 / pow2(p.DX[d]));
          }
          break;
        default:
          abort();
        }
      });
}

} // namespace Poisson2
