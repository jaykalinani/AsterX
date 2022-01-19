#include <fixmath.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop.hxx>

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

namespace WaveToyCPU {
using namespace std;
using namespace Loop;

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
template <typename T> T spline_potential(T r) {
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
template <typename F, typename T> auto timederiv(const F &f, T dt) {
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
CCTK_REAL standing(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z) {
  DECLARE_CCTK_PARAMETERS;
  auto kx = 2 * M_PI * spatial_frequency_x;
  auto ky = 2 * M_PI * spatial_frequency_y;
  auto kz = 2 * M_PI * spatial_frequency_z;
  auto omega = sqrt(pow(kx, 2) + pow(ky, 2) + pow(kz, 2));
  return cos(omega * t) * cos(kx * x) * cos(ky * y) * cos(kz * z);
}

// Periodic Gaussian
CCTK_REAL periodic_gaussian(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                            CCTK_REAL z) {
  DECLARE_CCTK_PARAMETERS;
  auto kx = M_PI * spatial_frequency_x;
  auto ky = M_PI * spatial_frequency_y;
  auto kz = M_PI * spatial_frequency_z;
  auto omega = sqrt(pow(kx, 2) + pow(ky, 2) + pow(kz, 2));
  return exp(-0.5 * pow(sin(kx * x + ky * y + kz * z - omega * t) / width, 2));
}

// Gaussian
CCTK_REAL gaussian(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z) {
  DECLARE_CCTK_PARAMETERS;
  // u(t,r) = (f(r-t) - f(r+t)) / r
  auto r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
  auto f = [&](auto x) { return exp(-0.5 * pow(x / width, 2)); };
  auto fx = [&](auto x) { return -x / pow(width, 2) * f(x); };
  if (r < 1.0e-8)
    // Use L'Hôpital's rule for small r
    return fx(r - t) - fx(r + t);
  else
    return (f(r - t) - f(r + t)) / r;
}

// Central potential
CCTK_REAL central_potential(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                            CCTK_REAL z) {
  DECLARE_CCTK_PARAMETERS;
  auto r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
  return -central_point_charge / pow(central_point_radius, dim) *
         spline_potential(r / central_point_radius);
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void WaveToyCPU_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyCPU_Initialize;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;
  const CCTK_REAL dt = CCTK_DELTA_TIME;

  const array<int, dim> indextype = {1, 1, 1};
  const GF3D2layout layout(cctkGH, indextype);
  const GF3D2<CCTK_REAL> gf_phi(layout, phi);
  const GF3D2<CCTK_REAL> gf_psi(layout, psi);

  if (CCTK_EQUALS(initial_condition, "standing wave")) {

    loop_int<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
      gf_phi(p.I) = standing(t, p.x, p.y, p.z);
      gf_psi(p.I) = timederiv(standing, dt)(t, p.x, p.y, p.z);
    });

  } else if (CCTK_EQUALS(initial_condition, "periodic Gaussian")) {

    loop_int<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
      gf_phi(p.I) = periodic_gaussian(t, p.x, p.y, p.z);
      gf_psi(p.I) = timederiv(periodic_gaussian, dt)(t, p.x, p.y, p.z);
    });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {

    loop_int<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
      gf_phi(p.I) = gaussian(t, p.x, p.y, p.z);
      gf_psi(p.I) = timederiv(gaussian, dt)(t, p.x, p.y, p.z);
    });

  } else if (CCTK_EQUALS(initial_condition, "central potential")) {

    loop_int<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
      gf_phi(p.I) = central_potential(t, p.x, p.y, p.z);
      gf_psi(p.I) = timederiv(central_potential, dt)(t, p.x, p.y, p.z);
    });

  } else {
    assert(0);
  }
}

extern "C" void WaveToyCPU_Sync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyCPU_Sync;
  DECLARE_CCTK_PARAMETERS;

  // Do nothing
}

extern "C" void WaveToyCPU_Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyCPU_Boundaries;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;
  const CCTK_REAL dt = CCTK_DELTA_TIME;

  const array<int, dim> indextype = {1, 1, 1};
  const GF3D2layout layout(cctkGH, indextype);
  const GF3D2<CCTK_REAL> gf_phi(layout, phi);
  const GF3D2<CCTK_REAL> gf_psi(layout, psi);

  if (CCTK_EQUALS(boundary_condition, "none")) {

    // do nothing

  } else if (CCTK_EQUALS(boundary_condition, "initial")) {

    if (CCTK_EQUALS(initial_condition, "central potential")) {

      loop_bnd<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
        gf_phi(p.I) = central_potential(t, p.x, p.y, p.z);
        gf_psi(p.I) = timederiv(central_potential, dt)(t, p.x, p.y, p.z);
      });

    } else {
      assert(0);
    }

  } else {
    assert(0);
  }
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void WaveToyCPU_EstimateError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyCPU_EstimateError;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {1, 1, 1};
  const GF3D2layout layout(cctkGH, indextype);
  const GF3D2<const CCTK_REAL> gf_phi(layout, phi);
  const GF3D2<const CCTK_REAL> gf_psi(layout, psi);
  const GF3D2<CCTK_REAL> gf_regrid_error(layout, regrid_error);

  loop_int<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
    CCTK_REAL base_phi = fabs(gf_phi(p.I)) + fabs(phi_abs);
    CCTK_REAL errx_phi =
        fabs(gf_phi(p.I - p.DI[0]) - 2 * gf_phi(p.I) + gf_phi(p.I + p.DI[0])) /
        base_phi;
    CCTK_REAL erry_phi =
        fabs(gf_phi(p.I - p.DI[1]) - 2 * gf_phi(p.I) + gf_phi(p.I + p.DI[1])) /
        base_phi;
    CCTK_REAL errz_phi =
        fabs(gf_phi(p.I - p.DI[2]) - 2 * gf_phi(p.I) + gf_phi(p.I + p.DI[2])) /
        base_phi;
    CCTK_REAL base_psi = fabs(gf_psi(p.I)) + fabs(psi_abs);
    CCTK_REAL errx_psi =
        fabs(gf_psi(p.I - p.DI[0]) - 2 * gf_psi(p.I) + gf_psi(p.I + p.DI[0])) /
        base_psi;
    CCTK_REAL erry_psi =
        fabs(gf_psi(p.I - p.DI[1]) - 2 * gf_psi(p.I) + gf_psi(p.I + p.DI[1])) /
        base_psi;
    CCTK_REAL errz_psi =
        fabs(gf_psi(p.I - p.DI[2]) - 2 * gf_psi(p.I) + gf_psi(p.I + p.DI[2])) /
        base_psi;
    gf_regrid_error(p.I) =
        errx_phi + erry_phi + errz_phi + errx_psi + erry_psi + errz_psi;
  });
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void WaveToyCPU_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyCPU_RHS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;

  const array<int, dim> indextype = {1, 1, 1};
  const GF3D2layout layout(cctkGH, indextype);
  const GF3D2<const CCTK_REAL> gf_phi(layout, phi);
  const GF3D2<const CCTK_REAL> gf_psi(layout, psi);
  const GF3D2<CCTK_REAL> gf_phirhs(layout, phirhs);
  const GF3D2<CCTK_REAL> gf_psirhs(layout, psirhs);

  loop_int<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
    CCTK_REAL ddx_phi =
        (gf_phi(p.I - p.DI[0]) - 2 * gf_phi(p.I) + gf_phi(p.I + p.DI[0])) /
        pow(p.dx, 2);
    CCTK_REAL ddy_phi =
        (gf_phi(p.I - p.DI[1]) - 2 * gf_phi(p.I) + gf_phi(p.I + p.DI[1])) /
        pow(p.dy, 2);
    CCTK_REAL ddz_phi =
        (gf_phi(p.I - p.DI[2]) - 2 * gf_phi(p.I) + gf_phi(p.I + p.DI[2])) /
        pow(p.dz, 2);
    gf_psirhs(p.I) = ddx_phi + ddy_phi + ddz_phi - pow(mass, 2) * gf_phi(p.I) +
                     4 * M_PI * central_potential(t, p.x, p.y, p.z);
    gf_phirhs(p.I) = gf_psi(p.I);
  });
}

extern "C" void WaveToyCPU_RHSSync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyCPU_RHSSync;
  DECLARE_CCTK_PARAMETERS;

  // Do nothing
}

extern "C" void WaveToyCPU_RHSBoundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyCPU_RHSBoundaries;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {1, 1, 1};
  const GF3D2layout layout(cctkGH, indextype);
  const GF3D2<CCTK_REAL> gf_phirhs(layout, phirhs);
  const GF3D2<CCTK_REAL> gf_psirhs(layout, psirhs);

  // Dirichlet boundary conditions
  loop_bnd<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
    gf_psirhs(p.I) = 0;
    gf_phirhs(p.I) = 0;
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
extern "C" void WaveToyCPU_RadiativeBoundaries(CCTK_ARGUMENTS) {
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
              gf_psi(p.I) = -phi_vel * (r * gradphi + (phi_p[p.idx] - phi_inf) / r);
              gf_phi(p.I) = phi_p[p.idx] + dt * gf_psi(p.I);
            }
          }
        }
      }
    }
  }
}
#endif

extern "C" void WaveToyCPU_Energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyCPU_Energy;

  const array<int, dim> indextype = {1, 1, 1};
  const GF3D2layout layout(cctkGH, indextype);
  const GF3D2<const CCTK_REAL> gf_phi(layout, phi);
  const GF3D2<const CCTK_REAL> gf_psi(layout, psi);
  const GF3D2<CCTK_REAL> gf_eps(layout, eps);

  loop_int<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
    CCTK_REAL dt_phi = gf_psi(p.I);
    CCTK_REAL dx_phi =
        (gf_phi(p.I + p.DI[0]) - gf_phi(p.I - p.DI[0])) / (2 * p.dx);
    CCTK_REAL dy_phi =
        (gf_phi(p.I + p.DI[1]) - gf_phi(p.I - p.DI[1])) / (2 * p.dy);
    CCTK_REAL dz_phi =
        (gf_phi(p.I + p.DI[2]) - gf_phi(p.I - p.DI[2])) / (2 * p.dz);
    gf_eps(p.I) =
        (pow(dt_phi, 2) + pow(dx_phi, 2) + pow(dy_phi, 2) + pow(dz_phi, 2)) / 2;
  });
}

extern "C" void WaveToyCPU_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyCPU_Error;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;
  const CCTK_REAL dt = CCTK_DELTA_TIME;

  const array<int, dim> indextype = {1, 1, 1};
  const GF3D2layout layout(cctkGH, indextype);
  const GF3D2<const CCTK_REAL> gf_phi(layout, phi);
  const GF3D2<const CCTK_REAL> gf_psi(layout, psi);
  const GF3D2<CCTK_REAL> gf_phierr(layout, phierr);
  const GF3D2<CCTK_REAL> gf_psierr(layout, psierr);

  if (CCTK_EQUALS(initial_condition, "standing wave")) {

    loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
      gf_phierr(p.I) = gf_phi(p.I) - standing(t, p.x, p.y, p.z);
      gf_psierr(p.I) = gf_psi(p.I) - timederiv(standing, dt)(t, p.x, p.y, p.z);
    });

  } else if (CCTK_EQUALS(initial_condition, "periodic Gaussian")) {

    loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
      gf_phierr(p.I) = gf_phi(p.I) - periodic_gaussian(t, p.x, p.y, p.z);
      gf_psierr(p.I) =
          gf_psi(p.I) - timederiv(periodic_gaussian, dt)(t, p.x, p.y, p.z);
    });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {

    loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
      gf_phierr(p.I) = gf_phi(p.I) - gaussian(t, p.x, p.y, p.z);
      gf_psierr(p.I) = gf_psi(p.I) - timederiv(gaussian, dt)(t, p.x, p.y, p.z);
    });

  } else if (CCTK_EQUALS(initial_condition, "central potential")) {

    loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
      gf_phierr(p.I) = gf_phi(p.I) - central_potential(t, p.x, p.y, p.z);
      gf_psierr(p.I) =
          gf_psi(p.I) - timederiv(central_potential, dt)(t, p.x, p.y, p.z);
    });

  } else {
    assert(0);
  }
}

} // namespace WaveToyCPU
