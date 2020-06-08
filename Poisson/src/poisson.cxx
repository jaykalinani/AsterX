#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
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

  SolvePoisson();
}

// extern "C" void Poisson_fixup(CCTK_ARGUMENTS) {
//   DECLARE_CCTK_ARGUMENTS_Poisson_fixup;
//   DECLARE_CCTK_PARAMETERS;
//
//   const GF3D<CCTK_REAL, 0, 0, 0> phi_(cctkGH, phi);
//
//   // Set boundary conditions
//   loop_bnd<0, 0, 0>(
//       cctkGH, [&](const PointDesc &p) { phi_(p.I) = fbnd(p.x, p.y, p.z); });
// }

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
      ires_(p.I) = (phi_(p.I - p.DI(0)) - 2 * phi_(p.I) + phi_(p.I + p.DI(0))) /
                       pow(p.DX[0], 2) +
                   (phi_(p.I - p.DI(1)) - 2 * phi_(p.I) + phi_(p.I + p.DI(1))) /
                       pow(p.DX[1], 2) +
                   (phi_(p.I - p.DI(2)) - 2 * phi_(p.I) + phi_(p.I + p.DI(2))) /
                       pow(p.DX[2], 2) -
                   frhs(p.x, p.y, p.z);
  });
}

extern "C" void Poisson_residual_boundary(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Poisson_residual_boundary;
  DECLARE_CCTK_PARAMETERS;

  const GF3D<const CCTK_REAL, 0, 0, 0> phi_(cctkGH, phi);

  const GF3D<CCTK_REAL, 0, 0, 0> ires_(cctkGH, ires);

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    // ires_(p.I) = phi_(p.I) - fbnd(p.x, p.y, p.z);
    ires_(p.I) = 0.0;
  });
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
