#include <loop.hxx>
#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>

namespace Poisson {
using namespace Loop;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

constexpr CCTK_REAL rmin = 0.1;

constexpr CCTK_REAL A = 1.0;
constexpr CCTK_REAL W = 0.1;

// Gaussian RHS
template <typename T> constexpr T frhs(T x, T y, T z) {
  T r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
  return 4 * M_PI * A * exp(-pow(r / (sqrt(2) * W), 2));
}

// Solution for the above RHS, assuming homogeneous Dirichlet boundary
// conditions at infinity
template <typename T> constexpr T fsol(T x, T y, T z) {
  T r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
  T C = -A * pow(sqrt(2 * M_PI) * W, 3);
  if (r < 1.0e-12)
    return C * 2 / (sqrt(2 * M_PI) * W);
  return C * erf(r / (sqrt(2) * W)) / r;
}

template <typename T> constexpr T fbnd(T x, T y, T z) {
  T r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
  T alpha = -A * pow(sqrt(2 * M_PI) * W, 3);
  return alpha / fmax(rmin, r);
}

// Robin boundary conditions:
//   f(r) = \alpha / r + C
//   f'(r) = - \alpha / r^2
//
//   r f(r) - r C = \alpha
//   -r^2 f'(r) = \alpha
//
//   r f(r) - r C = -r^2 f'(r)
//
//   f(r) + r f'(r) = C

// AMReX uses a weird Laplace stencil. It takes a second derivative in
// one direction, and smoothes with weights [1, 4, 1] / 6 in the other
// directions.
template <typename T>
T laplace(const GF3D<const T, 0, 0, 0> &u, const PointDesc &p) {
  T r = 0;
  for (int dir = 0; dir < dim; ++dir) {
    vect<T, dim> lap{1, -2, 1};
    lap /= pow(p.DX[dir], 2);
    vect<T, dim> smo{1, 4, 1};
    smo /= 6;
    for (int k = -1; k <= 1; ++k) {
      for (int j = -1; j <= 1; ++j) {
        for (int i = -1; i <= 1; ++i) {
          const vect<int, dim> di{i, j, k};
          T w = 1;
          for (int d = 0; d < dim; ++d)
            if (d == dir)
              w *= lap[1 + di[d]];
            else
              w *= smo[1 + di[d]];
          auto I = p.I;
          for (int d = 0; d < dim; ++d)
            I += di[d] * p.DI[d];
          r += w * u(I);
        }
      }
    }
  }
  return r;
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void Poisson_setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Poisson_setup;
  DECLARE_CCTK_PARAMETERS;

  const GF3D<CCTK_REAL, 0, 0, 0> phi_(cctkGH, phi);
  const GF3D<CCTK_REAL, 0, 0, 0> rhs_(cctkGH, rhs);

  // Note: The boundary conditions are applied on the outermost layer
  // of the interior points, not on the boundary points.

  // Initialize all points (including ghosts) and set boundary conditions
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { phi_(p.I) = 0.0; });

  // Set boundary conditions
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    if (cctk_lbnd[0] + p.i == 1 || cctk_lbnd[0] + p.i == cctk_gsh[0] - 1 ||
        cctk_lbnd[1] + p.j == 1 || cctk_lbnd[1] + p.j == cctk_gsh[1] - 1 ||
        cctk_lbnd[2] + p.k == 1 || cctk_lbnd[2] + p.k == cctk_gsh[2] - 1)
      phi_(p.I) = fbnd(p.x, p.y, p.z);
  });

  // Set RHS
  loop_all<0, 0, 0>(
      cctkGH, [&](const PointDesc &p) { rhs_(p.I) = frhs(p.x, p.y, p.z); });
}

extern "C" void Poisson_solve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Poisson_solve;
  DECLARE_CCTK_PARAMETERS;

  const int gi_sol = CCTK_GroupIndex("Poisson::phi");
  assert(gi_sol >= 0);
  const int gi_rhs = CCTK_GroupIndex("Poisson::rhs");
  assert(gi_rhs >= 0);
  const int gi_res = CCTK_GroupIndex("Poisson::res");
  assert(gi_res >= 0);

  const CCTK_REAL reltol = 0.0;
  const CCTK_REAL abstol = 1.0e-12;
  CCTK_REAL res_initial, res_final;
  SolvePoisson(gi_sol, gi_rhs, gi_res, reltol, abstol, &res_initial,
               &res_final);
  CCTK_VINFO("Residual before solve: %g", double(res_initial));
  CCTK_VINFO("Residual after solve:  %g", double(res_final));
}

extern "C" void Poisson_residual(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Poisson_residual;
  DECLARE_CCTK_PARAMETERS;

  const GF3D<const CCTK_REAL, 0, 0, 0> phi_(cctkGH, phi);

  const GF3D<CCTK_REAL, 0, 0, 0> ires_(cctkGH, ires);

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    if (cctk_lbnd[0] + p.i == 1 || cctk_lbnd[0] + p.i == cctk_gsh[0] - 1 ||
        cctk_lbnd[1] + p.j == 1 || cctk_lbnd[1] + p.j == cctk_gsh[1] - 1 ||
        cctk_lbnd[2] + p.k == 1 || cctk_lbnd[2] + p.k == cctk_gsh[2] - 1)
      ires_(p.I) = phi_(p.I) - fbnd(p.x, p.y, p.z);
    else
      ires_(p.I) = laplace(phi_, p) - frhs(p.x, p.y, p.z);
  });
}

extern "C" void Poisson_residual_boundary(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Poisson_residual_boundary;
  DECLARE_CCTK_PARAMETERS;

  const GF3D<const CCTK_REAL, 0, 0, 0> phi_(cctkGH, phi);

  const GF3D<CCTK_REAL, 0, 0, 0> ires_(cctkGH, ires);

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { ires_(p.I) = 0.0; });
}

extern "C" void Poisson_error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Poisson_error;
  DECLARE_CCTK_PARAMETERS;

  const GF3D<const CCTK_REAL, 0, 0, 0> phi_(cctkGH, phi);

  const GF3D<CCTK_REAL, 0, 0, 0> asol_(cctkGH, asol);
  const GF3D<CCTK_REAL, 0, 0, 0> aerr_(cctkGH, aerr);

  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    if (cctk_lbnd[0] + p.i == 1 || cctk_lbnd[0] + p.i == cctk_gsh[0] - 1 ||
        cctk_lbnd[1] + p.j == 1 || cctk_lbnd[1] + p.j == cctk_gsh[1] - 1 ||
        cctk_lbnd[2] + p.k == 1 || cctk_lbnd[2] + p.k == cctk_gsh[2] - 1)
      asol_(p.I) = 0.0;
    else
      asol_(p.I) = fsol(p.x, p.y, p.z);
    aerr_(p.I) = phi_(p.I) - asol_(p.I);
  });
}

} // namespace Poisson
