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

extern "C" void Poisson2_Source(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Poisson2_Source;
  DECLARE_CCTK_PARAMETERS;

  const int dim = Loop::dim;

  const int npoints = 27; // 3 levels

  const std::array<int, dim> indextype = {0, 0, 0};
  const std::array<int, dim> nghostzones = {
      cctk_nghostzones[0], cctk_nghostzones[1], cctk_nghostzones[2]};
  Arith::vect<int, dim> imin, imax;
  Loop::GridDescBase(cctkGH).box_int<0, 0, 0>(nghostzones, imin, imax);
  const Loop::GF3D2layout layout1(cctkGH, indextype);

  const Loop::GF3D2<CCTK_REAL> gf_src(layout1, src);

  const Loop::GridDescBase grid(cctkGH);
  if (CCTK_EQUALS(source, "constant")) {
    grid.loop_all<0, 0, 0>(grid.nghostzones,
                           [=](const Loop::PointDesc &p) ARITH_INLINE {
                             const Loop::GF3D2index index1(layout1, p.I);
                             gf_src(index1) = 1;
                           });
  } else if (CCTK_EQUALS(source, "logo")) {
    grid.loop_all<0, 0, 0>(
        grid.nghostzones, [=](const Loop::PointDesc &p) ARITH_INLINE {
          const Loop::GF3D2index index1(layout1, p.I);
          using std::lrint;
          Arith::vect<int, dim> i;
          for (int d = 0; d < dim; ++d)
            i[d] = lrint((p.X[d] / logo_width + 0.5) * npoints - 0.5);
          if (all(i >= Arith::vect<int, dim>::pure(0) &&
                  i < Arith::vect<int, dim>::pure(npoints))) {
            if (all(i / 9 == Arith::vect<int, dim>::pure(1))) {
              gf_src(index1) = 3; // red box
            } else if (all(i / 3 % 3 == Arith::vect<int, dim>::pure(1))) {
              gf_src(index1) = 2; // green box
            } else if (all(i / 1 % 3 == Arith::vect<int, dim>::pure(1))) {
              gf_src(index1) = 1; // blue box
            } else {
              gf_src(index1) = 0;
            }
          } else {
            gf_src(index1) = 0;
          }
        });
  } else {
    CCTK_ERROR("Unknown value for parameter \"source\"");
  }
}

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
  const Loop::GF3D2<const CCTK_REAL> gf_src(layout1, src);
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
          assert(0);
          abort();
        }
        gf_res(index1) = ddsol - gf_src(index1);
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
  const Loop::GF3D2<const CCTK_REAL> gf_src(layout1, sol);

  assert(PDESolvers::jacobians.has_value());
  PDESolvers::jacobian_t &J = PDESolvers::jacobians.value().get_local();

  const int fdorder1 = fdorder_jac < 0 ? fdorder : fdorder_jac;

  const Loop::GridDescBase grid(cctkGH);
  grid.loop_int<0, 0, 0>(
      grid.nghostzones, [=, &J](const Loop::PointDesc &p) ARITH_INLINE {
        const Loop::GF3D2index index1(layout1, p.I);
        using std::lrint;
        using Arith::pow2;
        const int idx = lrint(gf_idx(index1));
        switch (fdorder1) {
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
          assert(0);
          abort();
        }
      });
}

} // namespace Poisson2
