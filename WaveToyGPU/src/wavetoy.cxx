#include <fixmath.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

namespace WaveToyGPU {
using namespace std;

constexpr int dim = 3;

// Linear interpolation between (i0, x0) and (i1, x1)
template <typename Y, typename X> Y linterp(Y y0, Y y1, X x0, X x1, X x) {
  return Y(x - x0) / Y(x1 - x0) * y0 + Y(x - x1) / Y(x0 - x1) * y1;
}

// Spline with compact support of radius 1 and volume 1
template <typename T> T spline(T r) {
  if (r >= 1.0)
    return 0.0;
  constexpr CCTK_REAL f = dim == 1   ? 1.0
                          : dim == 2 ? 24.0 / 7.0 / M_PI
                          : dim == 3 ? 4.0 / M_PI
                                     : -1;
  const T r2 = pow(r, 2);
  return f * (r <= 0.5 ? 1 - 2 * r2 : 2 + r * (-4 + 2 * r));
}

// The potential for the spline
template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T
spline_potential(T r) {
  // \Laplace u = 4 \pi \rho
  if (r >= 1.0)
    return -1 / r;
  static_assert(dim == 3, "");
  const T r2 = pow(r, 2);
  return r <= 0.5
             ? -7 / T(3) + r2 * (8 / T(3) - 8 / T(5) * r2)
             : (1 / T(15) +
                r * (-40 / T(15) +
                     r2 * (80 / T(15) + r * (-80 / T(15) + r * 24 / T(15))))) /
                   r;
}

// Time derivative
// TODO: Use dual numbers for derivative
template <typename F, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST auto
timederiv(const F &f, T dt) {
  return [=](T t, T x, T y, T z) {
    return (f(t, x, y, z) - f(t - dt, x, y, z)) / dt;
  };
}

// Gradient
template <typename T> auto xderiv(T f(T t, T x, T y, T z), T dx) {
  return [=](T t, T x, T y, T z) {
    return (f(t, x + dx, y, z) - f(t, x - dx, y, z)) / (2 * dx);
  };
}
template <typename T> auto yderiv(T f(T t, T x, T y, T z), T dy) {
  return [=](T t, T x, T y, T z) {
    return (f(t, x, y + dy, z) - f(t, x, y - dy, z)) / (2 * dy);
  };
}
template <typename T> auto zderiv(T f(T t, T x, T y, T z), T dz) {
  return [=](T t, T x, T y, T z) {
    return (f(t, x, y, z + dz) - f(t, x, y, z - dz)) / (2 * dz);
  };
}

////////////////////////////////////////////////////////////////////////////////

// Standing wave
template <typename T> T standing(T t, T x, T y, T z) {
  DECLARE_CCTK_PARAMETERS;
  T kx = 2 * M_PI * spatial_frequency_x;
  T ky = 2 * M_PI * spatial_frequency_y;
  T kz = 2 * M_PI * spatial_frequency_z;
  T omega = sqrt(pow(kx, 2) + pow(ky, 2) + pow(kz, 2));
  return cos(omega * t) * cos(kx * x) * cos(ky * y) * cos(kz * z);
}

// Periodic Gaussian
template <typename T> T periodic_gaussian(T t, T x, T y, T z) {
  DECLARE_CCTK_PARAMETERS;
  T kx = M_PI * spatial_frequency_x;
  T ky = M_PI * spatial_frequency_y;
  T kz = M_PI * spatial_frequency_z;
  T omega = sqrt(pow(kx, 2) + pow(ky, 2) + pow(kz, 2));
  return exp(-0.5 * pow(sin(kx * x + ky * y + kz * z - omega * t) / width, 2));
}

// Gaussian
template <typename T> T gaussian(T t, T x, T y, T z) {
  DECLARE_CCTK_PARAMETERS;
  // u(t,r) = (f(r-t) - f(r+t)) / r
  T r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
  auto f = [&](T x) { return exp(-0.5 * pow(x / width, 2)); };
  auto fx = [&](T x) { return -x / pow(width, 2) * f(x); };
  if (r < 1.0e-8)
    // Use L'Hôpital's rule for small r
    return fx(r - t) - fx(r + t);
  else
    return (f(r - t) - f(r + t)) / r;
}

// Central potential
template <typename T> T central_potential(T t, T x, T y, T z) {
  DECLARE_CCTK_PARAMETERS;
  T r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
  return -central_point_charge / pow(central_point_radius, dim) *
         spline_potential(r / central_point_radius);
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void WaveToyGPU_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyGPU_Initialize;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;
  const CCTK_REAL dt = CCTK_DELTA_TIME;

  const auto standing =
      [=] CCTK_DEVICE CCTK_HOST(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                                CCTK_REAL z) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL kx = 2 * M_PI * spatial_frequency_x;
        CCTK_REAL ky = 2 * M_PI * spatial_frequency_y;
        CCTK_REAL kz = 2 * M_PI * spatial_frequency_z;
        CCTK_REAL omega = sqrt(pow(kx, 2) + pow(ky, 2) + pow(kz, 2));
        return cos(omega * t) * cos(kx * x) * cos(ky * y) * cos(kz * z);
      };

  const auto periodic_gaussian =
      [=] CCTK_DEVICE CCTK_HOST(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                                CCTK_REAL z) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL kx = M_PI * spatial_frequency_x;
        CCTK_REAL ky = M_PI * spatial_frequency_y;
        CCTK_REAL kz = M_PI * spatial_frequency_z;
        CCTK_REAL omega = sqrt(pow(kx, 2) + pow(ky, 2) + pow(kz, 2));
        return exp(-0.5 *
                   pow(sin(kx * x + ky * y + kz * z - omega * t) / width, 2));
      };

  const auto gaussian =
      [=] CCTK_DEVICE CCTK_HOST(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                                CCTK_REAL z) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // u(t,r) = (f(r-t) - f(r+t)) / r
        CCTK_REAL r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
        auto f = [=](CCTK_REAL x) { return exp(-0.5 * pow(x / width, 2)); };
        auto fx = [=](CCTK_REAL x) { return -x / pow(width, 2) * f(x); };
        if (r < 1.0e-8)
          // Use L'Hôpital's rule for small r
          return fx(r - t) - fx(r + t);
        else
          return (f(r - t) - f(r + t)) / r;
      };

  const auto central_potential =
      [=] CCTK_DEVICE CCTK_HOST(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                                CCTK_REAL z) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
        return -central_point_charge / pow(central_point_radius, dim) *
               spline_potential(r / central_point_radius);
      };

  // Determine loop extent
  const array<int, dim> nghostzones{cctkGH->cctk_nghostzones[0],
                                    cctkGH->cctk_nghostzones[1],
                                    cctkGH->cctk_nghostzones[2]};
  const Loop::GridDescBaseDevice griddesc(cctkGH);

  if (CCTK_EQUALS(initial_condition, "standing wave")) {

    griddesc.loop_int_device<1, 1, 1>(
        nghostzones, [=] CCTK_DEVICE CCTK_HOST(const Loop::PointDesc &p)
                         CCTK_ATTRIBUTE_ALWAYS_INLINE {
                           phi[p.idx] = standing(t, p.x, p.y, p.z);
                           psi[p.idx] =
                               timederiv(standing, dt)(t, p.x, p.y, p.z);
                         });

  } else if (CCTK_EQUALS(initial_condition, "periodic Gaussian")) {

    griddesc.loop_int_device<1, 1, 1>(
        nghostzones,
        [=] CCTK_DEVICE CCTK_HOST(const Loop::PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
              phi[p.idx] = periodic_gaussian(t, p.x, p.y, p.z);
              psi[p.idx] = timederiv(periodic_gaussian, dt)(t, p.x, p.y, p.z);
            });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {

    griddesc.loop_int_device<1, 1, 1>(
        nghostzones, [=] CCTK_DEVICE CCTK_HOST(const Loop::PointDesc &p)
                         CCTK_ATTRIBUTE_ALWAYS_INLINE {
                           phi[p.idx] = gaussian(t, p.x, p.y, p.z);
                           psi[p.idx] =
                               timederiv(gaussian, dt)(t, p.x, p.y, p.z);
                         });

  } else if (CCTK_EQUALS(initial_condition, "central potential")) {

    griddesc.loop_int_device<1, 1, 1>(
        nghostzones,
        [=] CCTK_DEVICE CCTK_HOST(const Loop::PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
              phi[p.idx] = central_potential(t, p.x, p.y, p.z);
              psi[p.idx] = timederiv(central_potential, dt)(t, p.x, p.y, p.z);
            });

  } else {
    assert(0);
  }
}

extern "C" void WaveToyGPU_Sync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyGPU_Sync;
  DECLARE_CCTK_PARAMETERS;

  // Do nothing
}

extern "C" void WaveToyGPU_Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyGPU_Boundaries;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;
  const CCTK_REAL dt = CCTK_DELTA_TIME;

  const auto central_potential =
      [=] CCTK_DEVICE CCTK_HOST(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                                CCTK_REAL z) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
        return -central_point_charge / pow(central_point_radius, dim) *
               spline_potential(r / central_point_radius);
      };

  // Determine loop extent
  const array<int, dim> nghostzones{cctkGH->cctk_nghostzones[0],
                                    cctkGH->cctk_nghostzones[1],
                                    cctkGH->cctk_nghostzones[2]};
  const Loop::GridDescBaseDevice griddesc(cctkGH);

  if (CCTK_EQUALS(boundary_condition, "none")) {

    // do nothing

  } else if (CCTK_EQUALS(boundary_condition, "initial")) {

    if (CCTK_EQUALS(initial_condition, "central potential")) {

      griddesc.loop_bnd_device<1, 1, 1>(
          nghostzones,
          [=] CCTK_DEVICE CCTK_HOST(const Loop::PointDesc &p)
              CCTK_ATTRIBUTE_ALWAYS_INLINE {
                phi[p.idx] = central_potential(t, p.x, p.y, p.z);
                psi[p.idx] = timederiv(central_potential, dt)(t, p.x, p.y, p.z);
              });

    } else {
      assert(0);
    }

  } else {
    assert(0);
  }
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void WaveToyGPU_EstimateError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyGPU_EstimateError;
  DECLARE_CCTK_PARAMETERS;

  // Determine loop extent
  const array<int, dim> nghostzones{cctkGH->cctk_nghostzones[0],
                                    cctkGH->cctk_nghostzones[1],
                                    cctkGH->cctk_nghostzones[2]};
  const Loop::GridDescBaseDevice griddesc(cctkGH);

  griddesc.loop_int_device<1, 1, 1>(
      nghostzones, [=] CCTK_DEVICE CCTK_HOST(
                       const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL base_phi = fabs(phi[p.idx]) + fabs(phi_abs);
        CCTK_REAL errx_phi =
            fabs(phi[p.idx - p.di] - 2 * phi[p.idx] + phi[p.idx + p.di]) /
            base_phi;
        CCTK_REAL erry_phi =
            fabs(phi[p.idx - p.dj] - 2 * phi[p.idx] + phi[p.idx + p.dj]) /
            base_phi;
        CCTK_REAL errz_phi =
            fabs(phi[p.idx - p.dk] - 2 * phi[p.idx] + phi[p.idx + p.dk]) /
            base_phi;
        CCTK_REAL base_psi = fabs(psi[p.idx]) + fabs(psi_abs);
        CCTK_REAL errx_psi =
            fabs(psi[p.idx - p.di] - 2 * psi[p.idx] + psi[p.idx + p.di]) /
            base_psi;
        CCTK_REAL erry_psi =
            fabs(psi[p.idx - p.dj] - 2 * psi[p.idx] + psi[p.idx + p.dj]) /
            base_psi;
        CCTK_REAL errz_psi =
            fabs(psi[p.idx - p.dk] - 2 * psi[p.idx] + psi[p.idx + p.dk]) /
            base_psi;
        regrid_error[p.idx] =
            errx_phi + erry_phi + errz_phi + errx_psi + erry_psi + errz_psi;
      });
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void WaveToyGPU_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyGPU_RHS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;
  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  const auto central_potential =
      [=] CCTK_DEVICE CCTK_HOST(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                                CCTK_REAL z) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
        return -central_point_charge / pow(central_point_radius, dim) *
               spline_potential(r / central_point_radius);
      };

  // Determine loop extent
  const array<int, dim> nghostzones{cctkGH->cctk_nghostzones[0],
                                    cctkGH->cctk_nghostzones[1],
                                    cctkGH->cctk_nghostzones[2]};
  const Loop::GridDescBaseDevice griddesc(cctkGH);

  griddesc.loop_int_device<1, 1, 1>(
      nghostzones, [=] CCTK_DEVICE CCTK_HOST(
                       const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL ddx_phi =
            (phi[p.idx - p.di] - 2 * phi[p.idx] + phi[p.idx + p.di]) /
            pow(dx, 2);
        CCTK_REAL ddy_phi =
            (phi[p.idx - p.dj] - 2 * phi[p.idx] + phi[p.idx + p.dj]) /
            pow(dy, 2);
        CCTK_REAL ddz_phi =
            (phi[p.idx - p.dk] - 2 * phi[p.idx] + phi[p.idx + p.dk]) /
            pow(dz, 2);
        psirhs[p.idx] = ddx_phi + ddy_phi + ddz_phi -
                        pow(mass, 2) * phi[p.idx] +
                        4 * M_PI * central_potential(t, p.x, p.y, p.z);
        phirhs[p.idx] = psi[p.idx];
      });
}

extern "C" void WaveToyGPU_RHSSync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyGPU_RHSSync;
  DECLARE_CCTK_PARAMETERS;

  // Do nothing
}

extern "C" void WaveToyGPU_RHSBoundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyGPU_RHSBoundaries;
  DECLARE_CCTK_PARAMETERS;

  // Determine loop extent
  const array<int, dim> nghostzones{cctkGH->cctk_nghostzones[0],
                                    cctkGH->cctk_nghostzones[1],
                                    cctkGH->cctk_nghostzones[2]};
  const Loop::GridDescBaseDevice griddesc(cctkGH);

  // Dirichlet boundary conditions
  griddesc.loop_bnd_device<1, 1, 1>(
      nghostzones, [=] CCTK_DEVICE CCTK_HOST(const Loop::PointDesc &p)
                       CCTK_ATTRIBUTE_ALWAYS_INLINE {
                         psirhs[p.idx] = 0;
                         phirhs[p.idx] = 0;
                       });
}

////////////////////////////////////////////////////////////////////////////////

// Miguel Alcubierre, Gabrielle Allen, Gerd Lanfermann
//
// Routines for applying radiation boundary conditions
//
// Taken from
// <Cactus/arrangements/CactusBase/Boundary/src/RadiationBoundary.c>
//
// The radiative boundary condition that is implemented is:
//
//   f  =  f0  +  u(r - v*t) / r  +  h(r + v*t) / r
//
// That is, I assume outgoing radial waves with a 1/r fall off, and
// the correct asymptotic value f0, plus I include the possibility of
// incoming waves as well (these incoming waves should be modeled
// somehow).
//
// The condition above leads to the differential equation:
//
//   (x / r) d f  +  v d f  + v x (f - f0) / r^2  =  v x H / r^2
//     i      t         i        i                      i
//
// where x_i is the normal direction to the given boundaries, and H =
// 2 dh(s)/ds.
//
// So at a given boundary I only worry about derivatives in the normal
// direction. Notice that u(r-v*t) has dissapeared, but we still do
// not know the value of H.
//
// To get H what I do is the following: I evaluate the expression one
// point in from the boundary and solve for H there. We now need a way
// of extrapolation H to the boundary. For this I assume that H falls
// off as a power law:
//
//   H = k/r**n  =>  d H  =  - n H/r
//                    i
//
// The value of n is is defined by the parameter "radpower". If this
// parameter is negative, H is forced to be zero (this corresponds to
// pure outgoing waves and is the default).
//
// The behaviour I have observed is the following: Using H=0 is very
// stable, but has a very bad initial transient. Taking n to be 0 or
// positive improves the initial behaviour considerably, but
// introduces a drift that can kill the evolution at very late times.
// Empirically, the best value I have found is n=2, for which the
// initial behaviour is very nice, and the late time drift is quite
// small.
//
// Another problem with this condition is that it does not use the
// physical characteristic speed, but rather it assumes a wave speed
// of v, so the boundaries should be out in the region where the
// characteristic speed is constant. Notice that this speed does not
// have to be 1. For gauge quantities {alpha, phi, trK} we can have a
// different asymptotic speed, which is why the value of v is passed
// as a parameter.
#if 0
extern "C" void WaveToyGPU_RadiativeBoundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  constexpr CCTK_REAL phi_inf = 0; // value at infinity
  constexpr CCTK_REAL phi_vel = 1; // propagation speed
  // constexpr CCTK_REAL psi_inf = 0;
  // constexpr CCTK_REAL psi_vel = 1;

  const CCTK_REAL dt = CCTK_DELTA_TIME;
  CCTK_REAL x0[dim], dx[dim];
  for (int d = 0; d < dim; ++d) {
    x0[d] = CCTK_ORIGIN_SPACE(d);
    dx[d] = CCTK_DELTA_SPACE(d);
  }

  int imin[dim], imax[dim];
  for (int d = 0; d < dim; ++d) {
    imin[d] = cctk_bbox[2 * d] ? 0 : cctk_nghostzones[d];
    imax[d] = cctk_lsh[d] - (cctk_bbox[2 * d + 1] ? 0 : cctk_nghostzones[d]);
  }
  // CCTK_VINFO("imin[%d,%d,%d]", imin[0], imin[1], imin[2]);
  // CCTK_VINFO("imax[%d,%d,%d]", imax[0], imax[1], imax[2]);

  int gbox[2 * dim]; // Are the boundaries (if any) ghosts zones?
  for (int d = 0; d < dim; ++d) {
    // #error "NEED lbnd>=0 ALWAYS"
    gbox[2 * d] = cctk_lbnd[d] > 0;
    gbox[2 * d + 1] = cctk_lbnd[d] + cctk_lsh[d] < cctk_gsh[d];
  }
  // CCTK_VINFO("lbnd[%d,%d,%d]", cctk_lbnd[0], cctk_lbnd[1], cctk_lbnd[2]);
  // CCTK_VINFO("lsh[%d,%d,%d]", cctk_lsh[0], cctk_lsh[1], cctk_lsh[2]);
  // CCTK_VINFO("gsh[%d,%d,%d]", cctk_gsh[0], cctk_gsh[1], cctk_gsh[2]);
  // CCTK_VINFO("gbox[%d,%d,%d,%d,%d,%d]", gbox[0], gbox[1], gbox[2], gbox[3],
  //            gbox[4], gbox[5]);

  constexpr int di = 1;
  const int dj = di * cctk_ash[0];
  const int dk = dj * cctk_ash[1];

  for (int dir = 0; dir < dim; ++dir) {
    for (int face = 0; face < 2; ++face) {
      if (!gbox[2 * dir + face]) {

        int bmin[dim], bmax[dim];
        for (int d = 0; d < dim; ++d) {
          // Skip edges and corners if d < dir
          bmin[d] = d < dir ? imin[d] : 0;
          bmax[d] = d < dir ? imax[d] : cctk_lsh[d];
        }
        if (face == 0)
          bmax[dir] = imin[dir];
        else
          bmin[dir] = imax[dir];
        // CCTK_VINFO("dir=%d face=%d", dir, face);
        // CCTK_VINFO("bmin[%d,%d,%d]", bmin[0], bmin[1], bmin[2]);
        // CCTK_VINFO("bmax[%d,%d,%d]", bmax[0], bmax[1], bmax[2]);

        for (int k = bmin[2]; k < bmax[2]; ++k) {
          for (int j = bmin[1]; j < bmax[1]; ++j) {
#pragma omp simd
            for (int i = bmin[0]; i < bmax[0]; ++i) {
              const int idx = CCTK_GFINDEX3D(cctkGH, i, j, k);
              CCTK_REAL x = x0[0] + (cctk_lbnd[0] + i) * dx[0];
              CCTK_REAL y = x0[1] + (cctk_lbnd[1] + j) * p.dy;
              CCTK_REAL z = x0[2] + (cctk_lbnd[2] + k) * p.dz;
              CCTK_REAL r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

              int ddir = dir == 0 ? di : dir == 1 ? dj : dk;
              int ndir = 2 * face - 1;

              CCTK_REAL gradphi =
                  (phi_p[p.idx] - phi_p[p.idx - ndir * ddir]) / dx[dir];
              psi[p.idx] = -phi_vel * (r * gradphi + (phi_p[p.idx] - phi_inf) / r);
              phi[p.idx] = phi_p[p.idx] + dt * psi[p.idx];
            }
          }
        }
      }
    }
  }
}
#endif

extern "C" void WaveToyGPU_Energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyGPU_Energy;

  // const CCTK_REAL dt = CCTK_DELTA_TIME;
  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  // Determine loop extent
  const Loop::GridDescBaseDevice griddesc(cctkGH);
  griddesc.loop_int_device<1, 1, 1>(
      griddesc.nghostzones,
      [=] CCTK_DEVICE CCTK_HOST(
          const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL dt_phi = psi[p.idx];
        CCTK_REAL dx_phi = (phi[p.idx + p.di] - phi[p.idx - p.di]) / (2 * dx);
        CCTK_REAL dy_phi = (phi[p.idx + p.dj] - phi[p.idx - p.dj]) / (2 * dy);
        CCTK_REAL dz_phi = (phi[p.idx + p.dk] - phi[p.idx - p.dk]) / (2 * dz);
        eps[p.idx] = (pow(dt_phi, 2) + pow(dx_phi, 2) + pow(dy_phi, 2) +
                      pow(dz_phi, 2)) /
                     2;
      });
}

extern "C" void WaveToyGPU_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyGPU_Error;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;
  const CCTK_REAL dt = CCTK_DELTA_TIME;

  const auto standing =
      [=] CCTK_DEVICE CCTK_HOST(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                                CCTK_REAL z) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL kx = 2 * M_PI * spatial_frequency_x;
        CCTK_REAL ky = 2 * M_PI * spatial_frequency_y;
        CCTK_REAL kz = 2 * M_PI * spatial_frequency_z;
        CCTK_REAL omega = sqrt(pow(kx, 2) + pow(ky, 2) + pow(kz, 2));
        return cos(omega * t) * cos(kx * x) * cos(ky * y) * cos(kz * z);
      };

  const auto periodic_gaussian =
      [=] CCTK_DEVICE CCTK_HOST(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                                CCTK_REAL z) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL kx = M_PI * spatial_frequency_x;
        CCTK_REAL ky = M_PI * spatial_frequency_y;
        CCTK_REAL kz = M_PI * spatial_frequency_z;
        CCTK_REAL omega = sqrt(pow(kx, 2) + pow(ky, 2) + pow(kz, 2));
        return exp(-0.5 *
                   pow(sin(kx * x + ky * y + kz * z - omega * t) / width, 2));
      };

  const auto gaussian =
      [=] CCTK_DEVICE CCTK_HOST(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                                CCTK_REAL z) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // u(t,r) = (f(r-t) - f(r+t)) / r
        CCTK_REAL r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
        auto f = [=](CCTK_REAL x) { return exp(-0.5 * pow(x / width, 2)); };
        auto fx = [=](CCTK_REAL x) { return -x / pow(width, 2) * f(x); };
        if (r < 1.0e-8)
          // Use L'Hôpital's rule for small r
          return fx(r - t) - fx(r + t);
        else
          return (f(r - t) - f(r + t)) / r;
      };

  const auto central_potential =
      [=] CCTK_DEVICE CCTK_HOST(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                                CCTK_REAL z) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
        return -central_point_charge / pow(central_point_radius, dim) *
               spline_potential(r / central_point_radius);
      };

  // Determine loop extent
  const array<int, dim> nghostzones{cctkGH->cctk_nghostzones[0],
                                    cctkGH->cctk_nghostzones[1],
                                    cctkGH->cctk_nghostzones[2]};
  const Loop::GridDescBaseDevice griddesc(cctkGH);

  if (CCTK_EQUALS(initial_condition, "standing wave")) {

    griddesc.loop_all_device<1, 1, 1>(
        nghostzones,
        [=] CCTK_DEVICE CCTK_HOST(const Loop::PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
              phierr[p.idx] = phi[p.idx] - standing(t, p.x, p.y, p.z);
              psierr[p.idx] =
                  psi[p.idx] - timederiv(standing, dt)(t, p.x, p.y, p.z);
            });

  } else if (CCTK_EQUALS(initial_condition, "periodic Gaussian")) {

    griddesc.loop_all_device<1, 1, 1>(
        nghostzones,
        [=] CCTK_DEVICE CCTK_HOST(const Loop::PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
              phierr[p.idx] = phi[p.idx] - periodic_gaussian(t, p.x, p.y, p.z);
              psierr[p.idx] = psi[p.idx] - timederiv(periodic_gaussian,
                                                     dt)(t, p.x, p.y, p.z);
            });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {

    griddesc.loop_all_device<1, 1, 1>(
        nghostzones,
        [=] CCTK_DEVICE CCTK_HOST(const Loop::PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
              phierr[p.idx] = phi[p.idx] - gaussian(t, p.x, p.y, p.z);
              psierr[p.idx] =
                  psi[p.idx] - timederiv(gaussian, dt)(t, p.x, p.y, p.z);
            });

  } else if (CCTK_EQUALS(initial_condition, "central potential")) {

    griddesc.loop_all_device<1, 1, 1>(
        nghostzones,
        [=] CCTK_DEVICE CCTK_HOST(const Loop::PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
              phierr[p.idx] = phi[p.idx] - central_potential(t, p.x, p.y, p.z);
              psierr[p.idx] = psi[p.idx] - timederiv(central_potential,
                                                     dt)(t, p.x, p.y, p.z);
            });

  } else {
    assert(0);
  }
}

} // namespace WaveToyGPU
