#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop.hxx>

#include <cassert>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <iostream>

namespace MaxwellToyAMReX {
using namespace std;

constexpr int dim = 3;

////////////////////////////////////////////////////////////////////////////////

template <typename T> struct potential {
  T phi, ax, ay, az;

  static potential pure(T a) { return {.phi = a, .ax = a, .ay = a, .az = a}; }
  template <typename F> potential map1(const F &f) const {
    const auto &x = *this;
    return {.phi = f(x.phi), .ax = f(x.ax), .ay = f(x.ay), .az = f(x.az)};
  }
  template <typename F> potential map2(const F &f, const potential &y) const {
    const auto &x = *this;
    return {.phi = f(x.phi, y.phi),
            .ax = f(x.ax, y.ax),
            .ay = f(x.ay, y.ay),
            .az = f(x.az, y.az)};
  }

  potential operator+() const { return *this; }
  potential operator-() const { return map1(std::negate<T>()); }
  potential operator+(const potential &y) const {
    return map2(std::plus<T>(), y);
  }
  potential operator-(const potential &y) const {
    return map2(std::minus<T>(), y);
  }
  potential operator*(const potential &y) const {
    return map2(std::multiplies<T>(), y);
  }
  potential operator/(const potential &y) const {
    return map2(std::divides<T>(), y);
  }
  potential operator+(T a) const { return *this + pure(a); }
  potential operator-(T a) const { return *this - pure(a); }
  potential operator*(T a) const { return *this * pure(a); }
  potential operator/(T a) const { return *this / pure(a); }
};

template <typename T> potential<T> plane_wave(T t, T x, T y, T z) {
  DECLARE_CCTK_PARAMETERS;
  T kx = 2 * M_PI * spatial_frequency_x;
  T ky = 2 * M_PI * spatial_frequency_y;
  T kz = 2 * M_PI * spatial_frequency_z;
  T omega = sqrt(pow(kx, 2) + pow(ky, 2) + pow(kz, 2));
  return {.phi = cos(omega * t - kx * x - ky * y - kz * z),
          .ax = omega / kx * cos(omega * t - kx * x - ky * y - kz * z),
          .ay = 0,
          .az = 0};
}

////////////////////////////////////////////////////////////////////////////////

// Derivatives
template <typename T, typename R> auto tderiv(R f(T t, T x, T y, T z), T dt) {
  return [=](T t, T x, T y, T z) {
    return (f(t + dt / 2, x, y, z) - f(t - dt / 2, x, y, z)) / dt;
  };
}
template <typename T, typename R> auto xderiv(R f(T t, T x, T y, T z), T dx) {
  return [=](T t, T x, T y, T z) {
    return (f(t, x + dx / 2, y, z) - f(t, x - dx / 2, y, z)) / dx;
  };
}
template <typename T, typename R> auto yderiv(R f(T t, T x, T y, T z), T dy) {
  return [=](T t, T x, T y, T z) {
    return (f(t, x, y + dy / 2, z) - f(t, x, y - dy / 2, z)) / dy;
  };
}
template <typename T, typename R> auto zderiv(R f(T t, T x, T y, T z), T dz) {
  return [=](T t, T x, T y, T z) {
    return (f(t, x, y, z + dz / 2) - f(t, x, y, z - dz / 2)) / dz;
  };
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void MaxwellToyAMReX_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;
  const CCTK_REAL dt = CCTK_DELTA_TIME;

  const Loop::GF3D<CCTK_REAL, 0, 0, 0> phi_(cctkGH, phi);

  const Loop::GF3D<CCTK_REAL, 1, 0, 0> ax_(cctkGH, ax);
  const Loop::GF3D<CCTK_REAL, 0, 1, 0> ay_(cctkGH, ay);
  const Loop::GF3D<CCTK_REAL, 0, 0, 1> az_(cctkGH, az);

  const Loop::GF3D<CCTK_REAL, 1, 0, 0> ex_(cctkGH, ex);
  const Loop::GF3D<CCTK_REAL, 0, 1, 0> ey_(cctkGH, ey);
  const Loop::GF3D<CCTK_REAL, 0, 0, 1> ez_(cctkGH, ez);

  const Loop::GF3D<CCTK_REAL, 0, 1, 1> bx_(cctkGH, bx);
  const Loop::GF3D<CCTK_REAL, 1, 0, 1> by_(cctkGH, by);
  const Loop::GF3D<CCTK_REAL, 1, 1, 0> bz_(cctkGH, bz);

  const Loop::GF3D<CCTK_REAL, 0, 0, 0> rho_(cctkGH, rho);

  const Loop::GF3D<CCTK_REAL, 1, 0, 0> jx_(cctkGH, jx);
  const Loop::GF3D<CCTK_REAL, 0, 1, 0> jy_(cctkGH, jy);
  const Loop::GF3D<CCTK_REAL, 0, 0, 1> jz_(cctkGH, jz);

  if (CCTK_EQUALS(initial_condition, "plane wave")) {

    Loop::loop_all<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      phi_(p.i, p.j, p.k) = plane_wave(t + dt / 2, p.x, p.y, p.z).phi;
    });

    Loop::loop_all<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      ax_(p.i, p.j, p.k) = plane_wave(t, p.x, p.y, p.z).ax;
    });
    Loop::loop_all<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      ay_(p.i, p.j, p.k) = plane_wave(t, p.x, p.y, p.z).ay;
    });
    Loop::loop_all<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      az_(p.i, p.j, p.k) = plane_wave(t, p.x, p.y, p.z).az;
    });

    Loop::loop_all<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      ex_(p.i, p.j, p.k) =
          -xderiv(plane_wave<CCTK_REAL>, p.dx)(t + dt / 2, p.x, p.y, p.z).phi -
          tderiv(plane_wave<CCTK_REAL>, dt)(t + dt / 2, p.x, p.y, p.z).ax;
    });
    Loop::loop_all<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      ey_(p.i, p.j, p.k) =
          -yderiv(plane_wave<CCTK_REAL>, p.dy)(t + dt / 2, p.x, p.y, p.z).phi -
          tderiv(plane_wave<CCTK_REAL>, dt)(t + dt / 2, p.x, p.y, p.z).ay;
    });
    Loop::loop_all<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      ez_(p.i, p.j, p.k) =
          -zderiv(plane_wave<CCTK_REAL>, p.dz)(t + dt / 2, p.x, p.y, p.z).phi -
          tderiv(plane_wave<CCTK_REAL>, dt)(t + dt / 2, p.x, p.y, p.z).az;
    });

    Loop::loop_all<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      bx_(p.i, p.j, p.k) =
          yderiv(plane_wave<CCTK_REAL>, p.dy)(t, p.x, p.y, p.z).az -
          zderiv(plane_wave<CCTK_REAL>, p.dz)(t, p.x, p.y, p.z).ay;
    });
    Loop::loop_all<1, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      by_(p.i, p.j, p.k) =
          zderiv(plane_wave<CCTK_REAL>, p.dz)(t, p.x, p.y, p.z).ax -
          xderiv(plane_wave<CCTK_REAL>, p.dx)(t, p.x, p.y, p.z).az;
    });
    Loop::loop_all<1, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      bz_(p.i, p.j, p.k) =
          xderiv(plane_wave<CCTK_REAL>, p.dx)(t, p.x, p.y, p.z).ay -
          yderiv(plane_wave<CCTK_REAL>, p.dy)(t, p.x, p.y, p.z).ax;
    });

    Loop::loop_all<0, 0, 0>(
        cctkGH, [&](const Loop::PointDesc &p) { rho_(p.i, p.j, p.k) = 0.0; });

    Loop::loop_all<1, 0, 0>(
        cctkGH, [&](const Loop::PointDesc &p) { jx_(p.i, p.j, p.k) = 0.0; });
    Loop::loop_all<0, 1, 0>(
        cctkGH, [&](const Loop::PointDesc &p) { jy_(p.i, p.j, p.k) = 0.0; });
    Loop::loop_all<0, 0, 1>(
        cctkGH, [&](const Loop::PointDesc &p) { jz_(p.i, p.j, p.k) = 0.0; });

  } else {
    assert(0);
  }
}

extern "C" void MaxwellToyAMReX_Evolve1(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dt = CCTK_DELTA_TIME;

  const Loop::GF3D<const CCTK_REAL, 0, 0, 0> phi_p_(cctkGH, phi_p);

  const Loop::GF3D<const CCTK_REAL, 1, 0, 0> ax_p_(cctkGH, ax_p);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 0> ay_p_(cctkGH, ay_p);
  const Loop::GF3D<const CCTK_REAL, 0, 0, 1> az_p_(cctkGH, az_p);

  const Loop::GF3D<const CCTK_REAL, 1, 0, 0> ex_p_(cctkGH, ex_p);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 0> ey_p_(cctkGH, ey_p);
  const Loop::GF3D<const CCTK_REAL, 0, 0, 1> ez_p_(cctkGH, ez_p);

  const Loop::GF3D<const CCTK_REAL, 0, 1, 1> bx_p_(cctkGH, bx_p);
  const Loop::GF3D<const CCTK_REAL, 1, 0, 1> by_p_(cctkGH, by_p);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 0> bz_p_(cctkGH, bz_p);

  const Loop::GF3D<const CCTK_REAL, 0, 0, 0> rho_p_(cctkGH, rho_p);

  const Loop::GF3D<const CCTK_REAL, 1, 0, 0> jx_p_(cctkGH, jx_p);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 0> jy_p_(cctkGH, jy_p);
  const Loop::GF3D<const CCTK_REAL, 0, 0, 1> jz_p_(cctkGH, jz_p);

  const Loop::GF3D<CCTK_REAL, 1, 0, 0> ax_(cctkGH, ax);
  const Loop::GF3D<CCTK_REAL, 0, 1, 0> ay_(cctkGH, ay);
  const Loop::GF3D<CCTK_REAL, 0, 0, 1> az_(cctkGH, az);

  const Loop::GF3D<CCTK_REAL, 0, 1, 1> bx_(cctkGH, bx);
  const Loop::GF3D<CCTK_REAL, 1, 0, 1> by_(cctkGH, by);
  const Loop::GF3D<CCTK_REAL, 1, 1, 0> bz_(cctkGH, bz);

  const Loop::GF3D<CCTK_REAL, 1, 0, 0> jx_(cctkGH, jx);
  const Loop::GF3D<CCTK_REAL, 0, 1, 0> jy_(cctkGH, jy);
  const Loop::GF3D<CCTK_REAL, 0, 0, 1> jz_(cctkGH, jz);

  Loop::loop_int<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    jx_(p.i, p.j, p.k) = jx_p_(p.i, p.j, p.k);
  });
  Loop::loop_int<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    jy_(p.i, p.j, p.k) = jy_p_(p.i, p.j, p.k);
  });
  Loop::loop_int<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    jz_(p.i, p.j, p.k) = jz_p_(p.i, p.j, p.k);
  });

  Loop::loop_int<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    bx_(p.i, p.j, p.k) =
        bx_p_(p.i, p.j, p.k) +
        dt * ((ey_p_(p.i, p.j, p.k + 1) - ey_p_(p.i, p.j, p.k)) / p.dz -
              (ez_p_(p.i, p.j + 1, p.k) - ez_p_(p.i, p.j, p.k)) / p.dy);
  });
  Loop::loop_int<1, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    by_(p.i, p.j, p.k) =
        by_p_(p.i, p.j, p.k) +
        dt * ((ez_p_(p.i + 1, p.j, p.k) - ez_p_(p.i, p.j, p.k)) / p.dx -
              (ex_p_(p.i, p.j, p.k + 1) - ex_p_(p.i, p.j, p.k)) / p.dz);
  });
  Loop::loop_int<1, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    bz_(p.i, p.j, p.k) =
        bz_p_(p.i, p.j, p.k) +
        dt * ((ex_p_(p.i, p.j + 1, p.k) - ex_p_(p.i, p.j, p.k)) / p.dy -
              (ey_p_(p.i + 1, p.j, p.k) - ey_p_(p.i, p.j, p.k)) / p.dx);
  });

  Loop::loop_int<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    ax_(p.i, p.j, p.k) =
        ax_p_(p.i, p.j, p.k) -
        dt * ((phi_p_(p.i + 1, p.j, p.k) - phi_p_(p.i, p.j, p.k)) / p.dx +
              ex_p_(p.i, p.j, p.k));
  });
  Loop::loop_int<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    ay_(p.i, p.j, p.k) =
        ay_p_(p.i, p.j, p.k) -
        dt * ((phi_p_(p.i, p.j + 1, p.k) - phi_p_(p.i, p.j, p.k)) / p.dy +
              ey_p_(p.i, p.j, p.k));
  });
  Loop::loop_int<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    az_(p.i, p.j, p.k) =
        az_p_(p.i, p.j, p.k) -
        dt * ((phi_p_(p.i, p.j, p.k + 1) - phi_p_(p.i, p.j, p.k)) / p.dz +
              ez_p_(p.i, p.j, p.k));
  });
}

extern "C" void MaxwellToyAMReX_Evolve2(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dt = CCTK_DELTA_TIME;

  const Loop::GF3D<const CCTK_REAL, 0, 0, 0> phi_p_(cctkGH, phi_p);

  const Loop::GF3D<const CCTK_REAL, 1, 0, 0> ex_p_(cctkGH, ex_p);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 0> ey_p_(cctkGH, ey_p);
  const Loop::GF3D<const CCTK_REAL, 0, 0, 1> ez_p_(cctkGH, ez_p);

  const Loop::GF3D<const CCTK_REAL, 0, 0, 0> rho_p_(cctkGH, rho_p);

  const Loop::GF3D<const CCTK_REAL, 1, 0, 0> ax_(cctkGH, ax);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 0> ay_(cctkGH, ay);
  const Loop::GF3D<const CCTK_REAL, 0, 0, 1> az_(cctkGH, az);

  const Loop::GF3D<const CCTK_REAL, 0, 1, 1> bx_(cctkGH, bx);
  const Loop::GF3D<const CCTK_REAL, 1, 0, 1> by_(cctkGH, by);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 0> bz_(cctkGH, bz);

  const Loop::GF3D<const CCTK_REAL, 1, 0, 0> jx_(cctkGH, jx);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 0> jy_(cctkGH, jy);
  const Loop::GF3D<const CCTK_REAL, 0, 0, 1> jz_(cctkGH, jz);

  const Loop::GF3D<CCTK_REAL, 0, 0, 0> phi_(cctkGH, phi);

  const Loop::GF3D<CCTK_REAL, 1, 0, 0> ex_(cctkGH, ex);
  const Loop::GF3D<CCTK_REAL, 0, 1, 0> ey_(cctkGH, ey);
  const Loop::GF3D<CCTK_REAL, 0, 0, 1> ez_(cctkGH, ez);

  const Loop::GF3D<CCTK_REAL, 0, 0, 0> rho_(cctkGH, rho);

  Loop::loop_int<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    rho_(p.i, p.j, p.k) =
        rho_p_(p.i, p.j, p.k) -
        dt * ((jx_(p.i, p.j, p.k) - jx_(p.i - 1, p.j, p.k)) / p.dx +
              (jy_(p.i, p.j, p.k) - jy_(p.i, p.j - 1, p.k)) / p.dy +
              (jz_(p.i, p.j, p.k) - jz_(p.i, p.j, p.k - 1)) / p.dz);
  });

  Loop::loop_int<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    ex_(p.i, p.j, p.k) =
        ex_p_(p.i, p.j, p.k) -
        dt * ((by_(p.i, p.j, p.k) - by_(p.i, p.j, p.k - 1)) / p.dz -
              (bz_(p.i, p.j, p.k) - bz_(p.i, p.j - 1, p.k)) / p.dy +
              jx_(p.i, p.j, p.k));
  });
  Loop::loop_int<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    ey_(p.i, p.j, p.k) =
        ey_p_(p.i, p.j, p.k) -
        dt * ((bz_(p.i, p.j, p.k) - bz_(p.i - 1, p.j, p.k)) / p.dx -
              (bx_(p.i, p.j, p.k) - bx_(p.i, p.j, p.k - 1)) / p.dz +
              jy_(p.i, p.j, p.k));
  });
  Loop::loop_int<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    ez_(p.i, p.j, p.k) =
        ez_p_(p.i, p.j, p.k) -
        dt * ((bx_(p.i, p.j, p.k) - bx_(p.i, p.j - 1, p.k)) / p.dy -
              (by_(p.i, p.j, p.k) - by_(p.i - 1, p.j, p.k)) / p.dx +
              jz_(p.i, p.j, p.k));
  });

  Loop::loop_int<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    phi_(p.i, p.j, p.k) =
        phi_p_(p.i, p.j, p.k) -
        dt * ((ax_(p.i, p.j, p.k) - ax_(p.i - 1, p.j, p.k)) / p.dx +
              (ay_(p.i, p.j, p.k) - ay_(p.i, p.j - 1, p.k)) / p.dy +
              (ay_(p.i, p.j, p.k) - az_(p.i, p.j, p.k - 1)) / p.dz);
  });
}

extern "C" void MaxwellToyAMReX_EstimateError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  Loop::loop_int<1, 1, 1>(
      cctkGH, [&](const Loop::PointDesc &p) { regrid_error[p.idx] = 0; });
}

extern "C" void MaxwellToyAMReX_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;
  const CCTK_REAL dt = CCTK_DELTA_TIME;

  const Loop::GF3D<CCTK_REAL, 0, 0, 0> phi_(cctkGH, phi);

  const Loop::GF3D<const CCTK_REAL, 1, 0, 0> ax_(cctkGH, ax);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 0> ay_(cctkGH, ay);
  const Loop::GF3D<const CCTK_REAL, 0, 0, 1> az_(cctkGH, az);

  const Loop::GF3D<const CCTK_REAL, 1, 0, 0> ex_(cctkGH, ex);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 0> ey_(cctkGH, ey);
  const Loop::GF3D<const CCTK_REAL, 0, 0, 1> ez_(cctkGH, ez);

  const Loop::GF3D<const CCTK_REAL, 0, 1, 1> bx_(cctkGH, bx);
  const Loop::GF3D<const CCTK_REAL, 1, 0, 1> by_(cctkGH, by);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 0> bz_(cctkGH, bz);

  const Loop::GF3D<const CCTK_REAL, 0, 0, 0> rho_(cctkGH, rho);

  const Loop::GF3D<const CCTK_REAL, 1, 0, 0> jx_(cctkGH, jx);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 0> jy_(cctkGH, jy);
  const Loop::GF3D<const CCTK_REAL, 0, 0, 1> jz_(cctkGH, jz);

  const Loop::GF3D<CCTK_REAL, 0, 0, 0> phierr_(cctkGH, phierr);

  const Loop::GF3D<CCTK_REAL, 1, 0, 0> axerr_(cctkGH, axerr);
  const Loop::GF3D<CCTK_REAL, 0, 1, 0> ayerr_(cctkGH, ayerr);
  const Loop::GF3D<CCTK_REAL, 0, 0, 1> azerr_(cctkGH, azerr);

  const Loop::GF3D<CCTK_REAL, 1, 0, 0> exerr_(cctkGH, exerr);
  const Loop::GF3D<CCTK_REAL, 0, 1, 0> eyerr_(cctkGH, eyerr);
  const Loop::GF3D<CCTK_REAL, 0, 0, 1> ezerr_(cctkGH, ezerr);

  const Loop::GF3D<CCTK_REAL, 0, 1, 1> bxerr_(cctkGH, bxerr);
  const Loop::GF3D<CCTK_REAL, 1, 0, 1> byerr_(cctkGH, byerr);
  const Loop::GF3D<CCTK_REAL, 1, 1, 0> bzerr_(cctkGH, bzerr);

  const Loop::GF3D<CCTK_REAL, 0, 0, 0> rhoerr_(cctkGH, rhoerr);

  const Loop::GF3D<CCTK_REAL, 1, 0, 0> jxerr_(cctkGH, jxerr);
  const Loop::GF3D<CCTK_REAL, 0, 1, 0> jyerr_(cctkGH, jyerr);
  const Loop::GF3D<CCTK_REAL, 0, 0, 1> jzerr_(cctkGH, jzerr);

  if (CCTK_EQUALS(initial_condition, "plane wave")) {

    Loop::loop_all<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      phierr_(p.i, p.j, p.k) =
          phi_(p.i, p.j, p.k) - plane_wave(t + dt / 2, p.x, p.y, p.z).phi;
    });

    Loop::loop_all<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      axerr_(p.i, p.j, p.k) =
          ax_(p.i, p.j, p.k) - plane_wave(t, p.x, p.y, p.z).ax;
    });
    Loop::loop_all<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      ayerr_(p.i, p.j, p.k) =
          ay_(p.i, p.j, p.k) - plane_wave(t, p.x, p.y, p.z).ay;
    });
    Loop::loop_all<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      azerr_(p.i, p.j, p.k) =
          az_(p.i, p.j, p.k) - plane_wave(t, p.x, p.y, p.z).az;
    });

    Loop::loop_all<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      exerr_(p.i, p.j, p.k) =
          ex_(p.i, p.j, p.k) -
          (-xderiv(plane_wave<CCTK_REAL>, p.dx)(t + dt / 2, p.x, p.y, p.z).phi -
           tderiv(plane_wave<CCTK_REAL>, dt)(t + dt / 2, p.x, p.y, p.z).ax);
    });
    Loop::loop_all<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      eyerr_(p.i, p.j, p.k) =
          ey_(p.i, p.j, p.k) -
          (-yderiv(plane_wave<CCTK_REAL>, p.dy)(t + dt / 2, p.x, p.y, p.z).phi -
           tderiv(plane_wave<CCTK_REAL>, dt)(t + dt / 2, p.x, p.y, p.z).ay);
    });
    Loop::loop_all<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      ezerr_(p.i, p.j, p.k) =
          ez_(p.i, p.j, p.k) -
          (-zderiv(plane_wave<CCTK_REAL>, p.dz)(t + dt / 2, p.x, p.y, p.z).phi -
           tderiv(plane_wave<CCTK_REAL>, dt)(t + dt / 2, p.x, p.y, p.z).az);
    });

    Loop::loop_all<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      bxerr_(p.i, p.j, p.k) =
          bx_(p.i, p.j, p.k) -
          (yderiv(plane_wave<CCTK_REAL>, p.dy)(t, p.x, p.y, p.z).az -
           zderiv(plane_wave<CCTK_REAL>, p.dz)(t, p.x, p.y, p.z).ay);
    });
    Loop::loop_all<1, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      byerr_(p.i, p.j, p.k) =
          by_(p.i, p.j, p.k) -
          (zderiv(plane_wave<CCTK_REAL>, p.dz)(t, p.x, p.y, p.z).ax -
           xderiv(plane_wave<CCTK_REAL>, p.dx)(t, p.x, p.y, p.z).az);
    });
    Loop::loop_all<1, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      bzerr_(p.i, p.j, p.k) =
          bz_(p.i, p.j, p.k) -
          (xderiv(plane_wave<CCTK_REAL>, p.dx)(t, p.x, p.y, p.z).ay -
           yderiv(plane_wave<CCTK_REAL>, p.dy)(t, p.x, p.y, p.z).ax);
    });

    Loop::loop_all<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      rhoerr_(p.i, p.j, p.k) = 0.0;
    });

    Loop::loop_all<1, 0, 0>(
        cctkGH, [&](const Loop::PointDesc &p) { jxerr_(p.i, p.j, p.k) = 0.0; });
    Loop::loop_all<0, 1, 0>(
        cctkGH, [&](const Loop::PointDesc &p) { jyerr_(p.i, p.j, p.k) = 0.0; });
    Loop::loop_all<0, 0, 1>(
        cctkGH, [&](const Loop::PointDesc &p) { jzerr_(p.i, p.j, p.k) = 0.0; });

  } else {
    assert(0);
  }
}

} // namespace MaxwellToyAMReX
