#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <iostream>

namespace MaxwellToyCarpetX {
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

////////////////////////////////////////////////////////////////////////////////

// Linear
template <typename T> potential<T> linear(T t, T x, T y, T z) {
  DECLARE_CCTK_PARAMETERS;
  T phi = t - x;
  // T phi = x;
  // T phi = pow(x, 2);
  // T phi = y;
  // T phi = pow(y, 2);
  return {.phi = phi, .ax = -phi, .ay = 0, .az = 0};
}

// Plane wave
template <typename T> potential<T> plane_wave(T t, T x, T y, T z) {
  DECLARE_CCTK_PARAMETERS;
  T kx = 2 * M_PI * spatial_frequency_x;
  T ky = 2 * M_PI * spatial_frequency_y;
  T kz = 2 * M_PI * spatial_frequency_z;
  T omega = sqrt(pow(kx, 2) + pow(ky, 2) + pow(kz, 2));
  T phi = cos(omega * t - kx * x - ky * y - kz * z);
  return {.phi = phi, .ax = omega / kx * phi, .ay = 0, .az = 0};
}

// Plane wave with Gaussian profile
template <typename T> potential<T> gaussian_wave(T t, T x, T y, T z) {
  DECLARE_CCTK_PARAMETERS;
  T kx = M_PI * spatial_frequency_x;
  T ky = M_PI * spatial_frequency_y;
  T kz = M_PI * spatial_frequency_z;
  T omega = sqrt(pow(kx, 2) + pow(ky, 2) + pow(kz, 2));
  T phi = exp(-pow(sin(omega * t - kx * x - ky * y - kz * z) / width, 2) / 2);
  return {.phi = phi, .ax = omega / kx * phi, .ay = 0, .az = 0};
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

extern "C" void MaxwellToyCarpetX_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MaxwellToyCarpetX_Initialize;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;
  const CCTK_REAL dt = CCTK_DELTA_TIME;

  const Loop::GF3D<CCTK_REAL, 0, 0, 0> phi_(cctkGH, phi);

  const Loop::GF3D<CCTK_REAL, 1, 0, 0> ax_(cctkGH, ax);
  const Loop::GF3D<CCTK_REAL, 0, 1, 0> ay_(cctkGH, ay);
  const Loop::GF3D<CCTK_REAL, 0, 0, 1> az_(cctkGH, az);

  const Loop::GF3D<CCTK_REAL, 0, 0, 0> rho_(cctkGH, rho);

  const Loop::GF3D<CCTK_REAL, 1, 0, 0> jx_(cctkGH, jx);
  const Loop::GF3D<CCTK_REAL, 0, 1, 0> jy_(cctkGH, jy);
  const Loop::GF3D<CCTK_REAL, 0, 0, 1> jz_(cctkGH, jz);

  if (CCTK_EQUALS(initial_condition, "linear")) {

    Loop::loop_all<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      phi_(p.i, p.j, p.k) = linear(t + dt / 2, p.x, p.y, p.z).phi;
    });

    Loop::loop_all<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      ax_(p.i, p.j, p.k) = linear(t, p.x, p.y, p.z).ax;
    });
    Loop::loop_all<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      ay_(p.i, p.j, p.k) = linear(t, p.x, p.y, p.z).ay;
    });
    Loop::loop_all<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      az_(p.i, p.j, p.k) = linear(t, p.x, p.y, p.z).az;
    });

    Loop::loop_all<0, 0, 0>(
        cctkGH, [&](const Loop::PointDesc &p) { rho_(p.i, p.j, p.k) = 0.0; });

    Loop::loop_all<1, 0, 0>(
        cctkGH, [&](const Loop::PointDesc &p) { jx_(p.i, p.j, p.k) = 0.0; });
    Loop::loop_all<0, 1, 0>(
        cctkGH, [&](const Loop::PointDesc &p) { jy_(p.i, p.j, p.k) = 0.0; });
    Loop::loop_all<0, 0, 1>(
        cctkGH, [&](const Loop::PointDesc &p) { jz_(p.i, p.j, p.k) = 0.0; });

  } else if (CCTK_EQUALS(initial_condition, "plane wave")) {

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

    Loop::loop_all<0, 0, 0>(
        cctkGH, [&](const Loop::PointDesc &p) { rho_(p.i, p.j, p.k) = 0.0; });

    Loop::loop_all<1, 0, 0>(
        cctkGH, [&](const Loop::PointDesc &p) { jx_(p.i, p.j, p.k) = 0.0; });
    Loop::loop_all<0, 1, 0>(
        cctkGH, [&](const Loop::PointDesc &p) { jy_(p.i, p.j, p.k) = 0.0; });
    Loop::loop_all<0, 0, 1>(
        cctkGH, [&](const Loop::PointDesc &p) { jz_(p.i, p.j, p.k) = 0.0; });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian wave")) {

    Loop::loop_all<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      phi_(p.i, p.j, p.k) = gaussian_wave(t + dt / 2, p.x, p.y, p.z).phi;
    });

    Loop::loop_all<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      ax_(p.i, p.j, p.k) = gaussian_wave(t, p.x, p.y, p.z).ax;
    });
    Loop::loop_all<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      ay_(p.i, p.j, p.k) = gaussian_wave(t, p.x, p.y, p.z).ay;
    });
    Loop::loop_all<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      az_(p.i, p.j, p.k) = gaussian_wave(t, p.x, p.y, p.z).az;
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

extern "C" void MaxwellToyCarpetX_InitializeFields(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MaxwellToyCarpetX_InitializeFields;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  const Loop::GF3D<const CCTK_REAL, 0, 0, 0> phi_(cctkGH, phi);

  const Loop::GF3D<const CCTK_REAL, 1, 0, 0> ax_(cctkGH, ax);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 0> ay_(cctkGH, ay);
  const Loop::GF3D<const CCTK_REAL, 0, 0, 1> az_(cctkGH, az);

  const Loop::GF3D<CCTK_REAL, 1, 0, 0> ex_(cctkGH, ex);
  const Loop::GF3D<CCTK_REAL, 0, 1, 0> ey_(cctkGH, ey);
  const Loop::GF3D<CCTK_REAL, 0, 0, 1> ez_(cctkGH, ez);

  const Loop::GF3D<CCTK_REAL, 0, 1, 1> bx_(cctkGH, bx);
  const Loop::GF3D<CCTK_REAL, 1, 0, 1> by_(cctkGH, by);
  const Loop::GF3D<CCTK_REAL, 1, 1, 0> bz_(cctkGH, bz);

  Loop::loop_int<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    ex_(p.i, p.j, p.k) = (phi_(p.i + 1, p.j, p.k) - phi_(p.i, p.j, p.k)) / dx;
  });
  Loop::loop_int<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    ey_(p.i, p.j, p.k) = (phi_(p.i, p.j + 1, p.k) - phi_(p.i, p.j, p.k)) / dy;
  });
  Loop::loop_int<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    ez_(p.i, p.j, p.k) = (phi_(p.i, p.j, p.k + 1) - phi_(p.i, p.j, p.k)) / dz;
  });

  Loop::loop_int<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    bx_(p.i, p.j, p.k) = (az_(p.i, p.j + 1, p.k) - az_(p.i, p.j, p.k)) / dy -
                         (ay_(p.i, p.j, p.k + 1) - ay_(p.i, p.j, p.k)) / dz;
  });
  Loop::loop_int<1, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    by_(p.i, p.j, p.k) = (ax_(p.i, p.j, p.k + 1) - ax_(p.i, p.j, p.k)) / dz -
                         (az_(p.i + 1, p.j, p.k) - az_(p.i, p.j, p.k)) / dx;
  });
  Loop::loop_int<1, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    bz_(p.i, p.j, p.k) = (ay_(p.i + 1, p.j, p.k) - ay_(p.i, p.j, p.k)) / dx -
                         (ax_(p.i, p.j + 1, p.k) - ax_(p.i, p.j, p.k)) / dy;
  });
}

extern "C" void MaxwellToyCarpetX_EstimateError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MaxwellToyCarpetX_EstimateError;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<const CCTK_REAL, 0, 0, 0> phi_(cctkGH, phi);

  const Loop::GF3D<const CCTK_REAL, 1, 0, 0> ax_(cctkGH, ax);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 0> ay_(cctkGH, ay);
  const Loop::GF3D<const CCTK_REAL, 0, 0, 1> az_(cctkGH, az);

  const Loop::GF3D<const CCTK_REAL, 1, 0, 0> ex_(cctkGH, ex);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 0> ey_(cctkGH, ey);
  const Loop::GF3D<const CCTK_REAL, 0, 0, 1> ez_(cctkGH, ez);

  const Loop::GF3D<const CCTK_REAL, 0, 1, 1> bx_(cctkGH, bx);
  const Loop::GF3D<const CCTK_REAL, 1, 0, 1> by_(cctkGH, by);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 0> bz_(cctkGH, bz);

#if 0
  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    regrid_error[p.idx] = fabs(p.x) < 0.5 && fabs(p.y) < 0.5 && fabs(p.z) <
    0.5;
  });
#endif

#if 0
  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    CCTK_REAL err = 0;
    for (int b = 0; b < 2; ++b)
      for (int a = 0; a < 2; ++a)
        err += fabs(ax_(p.i - 1, p.j + a, p.k + b) -
                    2 * ax_(p.i, p.j + a, p.k + b) +
                    ax_(p.i + 1, p.j + a, p.k + b)) +
               fabs(ay_(p.i + a, p.j - 1, p.k + b) -
                    2 * ay_(p.i + a, p.j, p.k + b) +
                    ay_(p.i + a, p.j + 1, p.k + b)) +
               fabs(az_(p.i + a, p.j + b, p.k - 1) -
                    2 * az_(p.i + a, p.j + b, p.k) +
                    az_(p.i + a, p.j + b, p.k + 1));
    regrid_error[p.idx] = err;
  });
#endif

#if 0
  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    auto closeto = [&](CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL r) {
      return fabs(p.x - x) <= r && fabs(p.y - y) <= r && fabs(p.z - z) <= r;
    };
    constexpr CCTK_REAL x = 0.25;
    constexpr CCTK_REAL y = 0.25;
    constexpr CCTK_REAL z = 0.25;
    constexpr CCTK_REAL r = 0.2;
    regrid_error[p.idx] = // closeto(-x, -y, -z, r) ||
        closeto(+x, -y, -z, r) || closeto(-x, +y, -z, r) ||
        closeto(+x, +y, -z, r) || closeto(-x, -y, +z, r) ||
        closeto(+x, -y, +z, r) || closeto(-x, +y, +z, r) ||
        closeto(+x, +y, +z, r);
  });
#endif

#if 1
  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    CCTK_REAL err = 0;
    const auto accum = [&](CCTK_REAL x) { err += fabs(x); };
    const auto diffx = [&](const auto &var_, int i, int j, int k) {
      for (int b = 0; b < 2; ++b)
        for (int a = 0; a < 2; ++a)
          accum(var_(i + 1, j + a, k + b) - var_(i, j + a, k + b));
    };
    const auto diffy = [&](const auto &var_, int i, int j, int k) {
      for (int b = 0; b < 2; ++b)
        for (int a = 0; a < 2; ++a)
          accum(var_(i + a, j + 1, k + b) - var_(i + a, j, k + b));
    };
    const auto diffz = [&](const auto &var_, int i, int j, int k) {
      for (int b = 0; b < 2; ++b)
        for (int a = 0; a < 2; ++a)
          accum(var_(i + a, j + b, k + 1) - var_(i + a, j + b, k));
    };
    const auto diff000 = [&](const auto &var_) {
      diffx(var_, p.i, p.j, p.k);
      diffy(var_, p.i, p.j, p.k);
      diffz(var_, p.i, p.j, p.k);
    };
    const auto diff100 = [&](const auto &var_) {
      diffx(var_, p.i - 1, p.j, p.k);
      diffx(var_, p.i, p.j, p.k);
      diffy(var_, p.i, p.j, p.k);
      diffz(var_, p.i, p.j, p.k);
    };
    const auto diff010 = [&](const auto &var_) {
      diffx(var_, p.i, p.j, p.k);
      diffy(var_, p.i, p.j - 1, p.k);
      diffy(var_, p.i, p.j, p.k);
      diffz(var_, p.i, p.j, p.k);
    };
    const auto diff001 = [&](const auto &var_) {
      diffx(var_, p.i, p.j, p.k);
      diffy(var_, p.i, p.j, p.k);
      diffz(var_, p.i, p.j, p.k - 1);
      diffz(var_, p.i, p.j, p.k);
    };
    const auto diff011 = [&](const auto &var_) {
      diffx(var_, p.i, p.j, p.k);
      diffy(var_, p.i, p.j - 1, p.k);
      diffy(var_, p.i, p.j, p.k);
      diffz(var_, p.i, p.j, p.k - 1);
      diffz(var_, p.i, p.j, p.k);
    };
    const auto diff101 = [&](const auto &var_) {
      diffx(var_, p.i - 1, p.j, p.k);
      diffx(var_, p.i, p.j, p.k);
      diffy(var_, p.i, p.j, p.k);
      diffz(var_, p.i, p.j, p.k - 1);
      diffz(var_, p.i, p.j, p.k);
    };
    const auto diff110 = [&](const auto &var_) {
      diffx(var_, p.i - 1, p.j, p.k);
      diffx(var_, p.i, p.j, p.k);
      diffy(var_, p.i, p.j - 1, p.k);
      diffy(var_, p.i, p.j, p.k);
      diffz(var_, p.i, p.j, p.k);
    };
    diff000(phi_);
    diff100(ax_);
    diff010(ay_);
    diff001(az_);
    diff100(ex_);
    diff010(ey_);
    diff001(ez_);
    diff011(bx_);
    diff101(by_);
    diff110(bz_);
    regrid_error[p.idx] = err;
  });
#endif
}

extern "C" void MaxwellToyCarpetX_Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MaxwellToyCarpetX_Boundaries;
  DECLARE_CCTK_PARAMETERS;

  // Do nothing; boundary conditions consist of synchronization only
}

extern "C" void MaxwellToyCarpetX_Evolve1(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MaxwellToyCarpetX_Evolve1;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

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

#if 1

  Loop::loop_int<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    jx_(p.i, p.j, p.k) = jx_p_(p.i, p.j, p.k);
  });
  Loop::loop_int<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    jy_(p.i, p.j, p.k) = jy_p_(p.i, p.j, p.k);
  });
  Loop::loop_int<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    jz_(p.i, p.j, p.k) = jz_p_(p.i, p.j, p.k);
  });

  if (evolve_b) {
    Loop::loop_int<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      bx_(p.i, p.j, p.k) =
          bx_p_(p.i, p.j, p.k) +
          dt * ((ey_p_(p.i, p.j, p.k + 1) - ey_p_(p.i, p.j, p.k)) / dz -
                (ez_p_(p.i, p.j + 1, p.k) - ez_p_(p.i, p.j, p.k)) / dy);
    });
    Loop::loop_int<1, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      by_(p.i, p.j, p.k) =
          by_p_(p.i, p.j, p.k) +
          dt * ((ez_p_(p.i + 1, p.j, p.k) - ez_p_(p.i, p.j, p.k)) / dx -
                (ex_p_(p.i, p.j, p.k + 1) - ex_p_(p.i, p.j, p.k)) / dz);
    });
    Loop::loop_int<1, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      bz_(p.i, p.j, p.k) =
          bz_p_(p.i, p.j, p.k) +
          dt * ((ex_p_(p.i, p.j + 1, p.k) - ex_p_(p.i, p.j, p.k)) / dy -
                (ey_p_(p.i + 1, p.j, p.k) - ey_p_(p.i, p.j, p.k)) / dx);
    });
  }

  Loop::loop_int<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    ax_(p.i, p.j, p.k) =
        ax_p_(p.i, p.j, p.k) -
        dt * ((phi_p_(p.i + 1, p.j, p.k) - phi_p_(p.i, p.j, p.k)) / dx +
              ex_p_(p.i, p.j, p.k));
  });
  Loop::loop_int<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    ay_(p.i, p.j, p.k) =
        ay_p_(p.i, p.j, p.k) -
        dt * ((phi_p_(p.i, p.j + 1, p.k) - phi_p_(p.i, p.j, p.k)) / dy +
              ey_p_(p.i, p.j, p.k));
  });
  Loop::loop_int<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    az_(p.i, p.j, p.k) =
        az_p_(p.i, p.j, p.k) -
        dt * ((phi_p_(p.i, p.j, p.k + 1) - phi_p_(p.i, p.j, p.k)) / dz +
              ez_p_(p.i, p.j, p.k));
  });

#else

  Loop::loop_int<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    jx_(p.i, p.j, p.k) = jx_p_(p.i, p.j, p.k);
  });
  Loop::loop_int<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    jy_(p.i, p.j, p.k) = jy_p_(p.i, p.j, p.k);
  });
  Loop::loop_int<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    jz_(p.i, p.j, p.k) = jz_p_(p.i, p.j, p.k);
  });

  if (evolve_b) {
    Loop::loop_int<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      bx_(p.i, p.j, p.k) = bx_p_(p.i, p.j, p.k);
    });
    Loop::loop_int<1, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      by_(p.i, p.j, p.k) = by_p_(p.i, p.j, p.k);
    });
    Loop::loop_int<1, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      bz_(p.i, p.j, p.k) = bz_p_(p.i, p.j, p.k);
    });
  }

  Loop::loop_int<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    ax_(p.i, p.j, p.k) = ax_p_(p.i, p.j, p.k);
  });
  Loop::loop_int<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    ay_(p.i, p.j, p.k) = ay_p_(p.i, p.j, p.k);
  });
  Loop::loop_int<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    az_(p.i, p.j, p.k) = az_p_(p.i, p.j, p.k);
  });

#endif
}

extern "C" void MaxwellToyCarpetX_Evolve2(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MaxwellToyCarpetX_Evolve2;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

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

#if 1

  Loop::loop_int<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    rho_(p.i, p.j, p.k) =
        rho_p_(p.i, p.j, p.k) -
        dt * ((jx_(p.i, p.j, p.k) - jx_(p.i - 1, p.j, p.k)) / dx +
              (jy_(p.i, p.j, p.k) - jy_(p.i, p.j - 1, p.k)) / dy +
              (jz_(p.i, p.j, p.k) - jz_(p.i, p.j, p.k - 1)) / dz);
  });

  Loop::loop_int<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    ex_(p.i, p.j, p.k) =
        ex_p_(p.i, p.j, p.k) -
        dt * ((by_(p.i, p.j, p.k) - by_(p.i, p.j, p.k - 1)) / dz -
              (bz_(p.i, p.j, p.k) - bz_(p.i, p.j - 1, p.k)) / dy +
              jx_(p.i, p.j, p.k));
  });
  Loop::loop_int<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    ey_(p.i, p.j, p.k) =
        ey_p_(p.i, p.j, p.k) -
        dt * ((bz_(p.i, p.j, p.k) - bz_(p.i - 1, p.j, p.k)) / dx -
              (bx_(p.i, p.j, p.k) - bx_(p.i, p.j, p.k - 1)) / dz +
              jy_(p.i, p.j, p.k));
  });
  Loop::loop_int<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    ez_(p.i, p.j, p.k) =
        ez_p_(p.i, p.j, p.k) -
        dt * ((bx_(p.i, p.j, p.k) - bx_(p.i, p.j - 1, p.k)) / dy -
              (by_(p.i, p.j, p.k) - by_(p.i - 1, p.j, p.k)) / dx +
              jz_(p.i, p.j, p.k));
  });

  Loop::loop_int<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    phi_(p.i, p.j, p.k) =
        phi_p_(p.i, p.j, p.k) -
        dt * ((ax_(p.i, p.j, p.k) - ax_(p.i - 1, p.j, p.k)) / dx +
              (ay_(p.i, p.j, p.k) - ay_(p.i, p.j - 1, p.k)) / dy +
              (az_(p.i, p.j, p.k) - az_(p.i, p.j, p.k - 1)) / dz);
  });

#else

  Loop::loop_int<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    rho_(p.i, p.j, p.k) = rho_p_(p.i, p.j, p.k);
  });

  Loop::loop_int<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    ex_(p.i, p.j, p.k) = ex_p_(p.i, p.j, p.k);
  });
  Loop::loop_int<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    ey_(p.i, p.j, p.k) = ey_p_(p.i, p.j, p.k);
  });
  Loop::loop_int<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    ez_(p.i, p.j, p.k) = ez_p_(p.i, p.j, p.k);
  });

  Loop::loop_int<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    phi_(p.i, p.j, p.k) = phi_p_(p.i, p.j, p.k);
  });

#endif
}

extern "C" void MaxwellToyCarpetX_Dependents1(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MaxwellToyCarpetX_Dependents1;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  const Loop::GF3D<const CCTK_REAL, 1, 0, 0> ax_(cctkGH, ax);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 0> ay_(cctkGH, ay);
  const Loop::GF3D<const CCTK_REAL, 0, 0, 1> az_(cctkGH, az);

  const Loop::GF3D<CCTK_REAL, 0, 1, 1> bx_(cctkGH, bx);
  const Loop::GF3D<CCTK_REAL, 1, 0, 1> by_(cctkGH, by);
  const Loop::GF3D<CCTK_REAL, 1, 1, 0> bz_(cctkGH, bz);

  Loop::loop_int<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    bx_(p.i, p.j, p.k) = (az_(p.i, p.j + 1, p.k) - az_(p.i, p.j, p.k)) / dy -
                         (ay_(p.i, p.j, p.k + 1) - ay_(p.i, p.j, p.k)) / dz;
  });
  Loop::loop_int<1, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    by_(p.i, p.j, p.k) = (ax_(p.i, p.j, p.k + 1) - ax_(p.i, p.j, p.k)) / dz -
                         (az_(p.i + 1, p.j, p.k) - az_(p.i, p.j, p.k)) / dx;
  });
  Loop::loop_int<1, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    bz_(p.i, p.j, p.k) = (ay_(p.i + 1, p.j, p.k) - ay_(p.i, p.j, p.k)) / dx -
                         (ax_(p.i, p.j + 1, p.k) - ax_(p.i, p.j, p.k)) / dy;
  });
}

extern "C" void MaxwellToyCarpetX_Constraints(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MaxwellToyCarpetX_Constraints;
  DECLARE_CCTK_PARAMETERS;

  if (analyse_every <= 0 || cctk_iteration % analyse_every != 0)
    return;

  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  const Loop::GF3D<const CCTK_REAL, 1, 0, 0> ex_(cctkGH, ex);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 0> ey_(cctkGH, ey);
  const Loop::GF3D<const CCTK_REAL, 0, 0, 1> ez_(cctkGH, ez);

  const Loop::GF3D<const CCTK_REAL, 0, 1, 1> bx_(cctkGH, bx);
  const Loop::GF3D<const CCTK_REAL, 1, 0, 1> by_(cctkGH, by);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 0> bz_(cctkGH, bz);

  const Loop::GF3D<const CCTK_REAL, 0, 0, 0> rho_(cctkGH, rho);

  const Loop::GF3D<CCTK_REAL, 0, 0, 0> dive_(cctkGH, dive);

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> divb_(cctkGH, divb);

  Loop::loop_int<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    dive_(p.i, p.j, p.k) = (ex_(p.i, p.j, p.k) - ex_(p.i - 1, p.j, p.k)) / dx +
                           (ey_(p.i, p.j, p.k) - ey_(p.i, p.j - 1, p.k)) / dy +
                           (ez_(p.i, p.j, p.k) - ez_(p.i, p.j, p.k - 1)) / dz -
                           rho_(p.i, p.j, p.k);
  });

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    divb_(p.i, p.j, p.k) = (bx_(p.i + 1, p.j, p.k) - bx_(p.i, p.j, p.k)) / dx +
                           (by_(p.i, p.j + 1, p.k) - by_(p.i, p.j, p.k)) / dy +
                           (bz_(p.i, p.j, p.k + 1) - bz_(p.i, p.j, p.k)) / dz;
  });
}

extern "C" void MaxwellToyCarpetX_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MaxwellToyCarpetX_Error;
  DECLARE_CCTK_PARAMETERS;

  if (analyse_every <= 0 || cctk_iteration % analyse_every != 0)
    return;

  const CCTK_REAL t = cctk_time;
  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  const Loop::GF3D<const CCTK_REAL, 0, 0, 0> phi_(cctkGH, phi);

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

  if (CCTK_EQUALS(initial_condition, "linear")) {

    Loop::loop_all<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      phierr_(p.i, p.j, p.k) =
          phi_(p.i, p.j, p.k) - linear(t + dt / 2, p.x, p.y, p.z).phi;
    });

    Loop::loop_all<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      axerr_(p.i, p.j, p.k) = ax_(p.i, p.j, p.k) - linear(t, p.x, p.y, p.z).ax;
    });
    Loop::loop_all<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      ayerr_(p.i, p.j, p.k) = ay_(p.i, p.j, p.k) - linear(t, p.x, p.y, p.z).ay;
    });
    Loop::loop_all<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      azerr_(p.i, p.j, p.k) = az_(p.i, p.j, p.k) - linear(t, p.x, p.y, p.z).az;
    });

    Loop::loop_all<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      exerr_(p.i, p.j, p.k) =
          ex_(p.i, p.j, p.k) -
          (-xderiv(linear<CCTK_REAL>, dx)(t + dt / 2, p.x, p.y, p.z).phi -
           tderiv(linear<CCTK_REAL>, dt)(t + dt / 2, p.x, p.y, p.z).ax);
    });
    Loop::loop_all<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      eyerr_(p.i, p.j, p.k) =
          ey_(p.i, p.j, p.k) -
          (-yderiv(linear<CCTK_REAL>, dy)(t + dt / 2, p.x, p.y, p.z).phi -
           tderiv(linear<CCTK_REAL>, dt)(t + dt / 2, p.x, p.y, p.z).ay);
    });
    Loop::loop_all<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      ezerr_(p.i, p.j, p.k) =
          ez_(p.i, p.j, p.k) -
          (-zderiv(linear<CCTK_REAL>, dz)(t + dt / 2, p.x, p.y, p.z).phi -
           tderiv(linear<CCTK_REAL>, dt)(t + dt / 2, p.x, p.y, p.z).az);
    });

    Loop::loop_all<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      bxerr_(p.i, p.j, p.k) =
          bx_(p.i, p.j, p.k) -
          (yderiv(linear<CCTK_REAL>, dy)(t, p.x, p.y, p.z).az -
           zderiv(linear<CCTK_REAL>, dz)(t, p.x, p.y, p.z).ay);
    });
    Loop::loop_all<1, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      byerr_(p.i, p.j, p.k) =
          by_(p.i, p.j, p.k) -
          (zderiv(linear<CCTK_REAL>, dz)(t, p.x, p.y, p.z).ax -
           xderiv(linear<CCTK_REAL>, dx)(t, p.x, p.y, p.z).az);
    });
    Loop::loop_all<1, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      bzerr_(p.i, p.j, p.k) =
          bz_(p.i, p.j, p.k) -
          (xderiv(linear<CCTK_REAL>, dx)(t, p.x, p.y, p.z).ay -
           yderiv(linear<CCTK_REAL>, dy)(t, p.x, p.y, p.z).ax);
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

  } else if (CCTK_EQUALS(initial_condition, "plane wave")) {

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
          (-xderiv(plane_wave<CCTK_REAL>, dx)(t + dt / 2, p.x, p.y, p.z).phi -
           tderiv(plane_wave<CCTK_REAL>, dt)(t + dt / 2, p.x, p.y, p.z).ax);
    });
    Loop::loop_all<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      eyerr_(p.i, p.j, p.k) =
          ey_(p.i, p.j, p.k) -
          (-yderiv(plane_wave<CCTK_REAL>, dy)(t + dt / 2, p.x, p.y, p.z).phi -
           tderiv(plane_wave<CCTK_REAL>, dt)(t + dt / 2, p.x, p.y, p.z).ay);
    });
    Loop::loop_all<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      ezerr_(p.i, p.j, p.k) =
          ez_(p.i, p.j, p.k) -
          (-zderiv(plane_wave<CCTK_REAL>, dz)(t + dt / 2, p.x, p.y, p.z).phi -
           tderiv(plane_wave<CCTK_REAL>, dt)(t + dt / 2, p.x, p.y, p.z).az);
    });

    Loop::loop_all<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      bxerr_(p.i, p.j, p.k) =
          bx_(p.i, p.j, p.k) -
          (yderiv(plane_wave<CCTK_REAL>, dy)(t, p.x, p.y, p.z).az -
           zderiv(plane_wave<CCTK_REAL>, dz)(t, p.x, p.y, p.z).ay);
    });
    Loop::loop_all<1, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      byerr_(p.i, p.j, p.k) =
          by_(p.i, p.j, p.k) -
          (zderiv(plane_wave<CCTK_REAL>, dz)(t, p.x, p.y, p.z).ax -
           xderiv(plane_wave<CCTK_REAL>, dx)(t, p.x, p.y, p.z).az);
    });
    Loop::loop_all<1, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      bzerr_(p.i, p.j, p.k) =
          bz_(p.i, p.j, p.k) -
          (xderiv(plane_wave<CCTK_REAL>, dx)(t, p.x, p.y, p.z).ay -
           yderiv(plane_wave<CCTK_REAL>, dy)(t, p.x, p.y, p.z).ax);
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

  } else if (CCTK_EQUALS(initial_condition, "Gaussian wave")) {

    Loop::loop_all<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      phierr_(p.i, p.j, p.k) =
          phi_(p.i, p.j, p.k) - gaussian_wave(t + dt / 2, p.x, p.y, p.z).phi;
    });

    Loop::loop_all<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      axerr_(p.i, p.j, p.k) =
          ax_(p.i, p.j, p.k) - gaussian_wave(t, p.x, p.y, p.z).ax;
    });
    Loop::loop_all<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      ayerr_(p.i, p.j, p.k) =
          ay_(p.i, p.j, p.k) - gaussian_wave(t, p.x, p.y, p.z).ay;
    });
    Loop::loop_all<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      azerr_(p.i, p.j, p.k) =
          az_(p.i, p.j, p.k) - gaussian_wave(t, p.x, p.y, p.z).az;
    });

    Loop::loop_all<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      exerr_(p.i, p.j, p.k) =
          ex_(p.i, p.j, p.k) -
          (-xderiv(gaussian_wave<CCTK_REAL>, dx)(t + dt / 2, p.x, p.y, p.z)
                .phi -
           tderiv(gaussian_wave<CCTK_REAL>, dt)(t + dt / 2, p.x, p.y, p.z).ax);
    });
    Loop::loop_all<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      eyerr_(p.i, p.j, p.k) =
          ey_(p.i, p.j, p.k) -
          (-yderiv(gaussian_wave<CCTK_REAL>, dy)(t + dt / 2, p.x, p.y, p.z)
                .phi -
           tderiv(gaussian_wave<CCTK_REAL>, dt)(t + dt / 2, p.x, p.y, p.z).ay);
    });
    Loop::loop_all<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      ezerr_(p.i, p.j, p.k) =
          ez_(p.i, p.j, p.k) -
          (-zderiv(gaussian_wave<CCTK_REAL>, dz)(t + dt / 2, p.x, p.y, p.z)
                .phi -
           tderiv(gaussian_wave<CCTK_REAL>, dt)(t + dt / 2, p.x, p.y, p.z).az);
    });

    Loop::loop_all<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      bxerr_(p.i, p.j, p.k) =
          bx_(p.i, p.j, p.k) -
          (yderiv(gaussian_wave<CCTK_REAL>, dy)(t, p.x, p.y, p.z).az -
           zderiv(gaussian_wave<CCTK_REAL>, dz)(t, p.x, p.y, p.z).ay);
    });
    Loop::loop_all<1, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      byerr_(p.i, p.j, p.k) =
          by_(p.i, p.j, p.k) -
          (zderiv(gaussian_wave<CCTK_REAL>, dz)(t, p.x, p.y, p.z).ax -
           xderiv(gaussian_wave<CCTK_REAL>, dx)(t, p.x, p.y, p.z).az);
    });
    Loop::loop_all<1, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      bzerr_(p.i, p.j, p.k) =
          bz_(p.i, p.j, p.k) -
          (xderiv(gaussian_wave<CCTK_REAL>, dx)(t, p.x, p.y, p.z).ay -
           yderiv(gaussian_wave<CCTK_REAL>, dy)(t, p.x, p.y, p.z).ax);
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

#if 0
extern "C" void MaxwellToyCarpetX_NaNCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MaxwellToyCarpetX_NaNCheck;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<const CCTK_REAL, 0, 0, 0> phi_(cctkGH, phi);

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

  const Loop::GF3D<const CCTK_REAL, 0, 0, 0> dive_(cctkGH, dive);

  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> divb_(cctkGH, divb);

  const auto nancheck = [&](const CCTK_REAL *restrict var,
                            const Loop::PointDesc &p) {
    if (isnan(var[p.idx]))
      CCTK_VERROR("Found nan at [%d,%d,%d]", p.i, p.j, p.k);
  };

  Loop::loop_all<0, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { nancheck(phi, p); });

  Loop::loop_all<1, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { nancheck(ax, p); });
  Loop::loop_all<0, 1, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { nancheck(ay, p); });
  Loop::loop_all<0, 0, 1>(cctkGH,
                          [&](const Loop::PointDesc &p) { nancheck(az, p); });

  Loop::loop_all<1, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { nancheck(ex, p); });
  Loop::loop_all<0, 1, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { nancheck(ey, p); });
  Loop::loop_all<0, 0, 1>(cctkGH,
                          [&](const Loop::PointDesc &p) { nancheck(ez, p); });

  Loop::loop_all<0, 1, 1>(cctkGH,
                          [&](const Loop::PointDesc &p) { nancheck(bx, p); });
  Loop::loop_all<1, 0, 1>(cctkGH,
                          [&](const Loop::PointDesc &p) { nancheck(by, p); });
  Loop::loop_all<1, 1, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { nancheck(bz, p); });

  Loop::loop_all<0, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { nancheck(rho, p); });

  Loop::loop_all<1, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { nancheck(jx, p); });
  Loop::loop_all<0, 1, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { nancheck(jy, p); });
  Loop::loop_all<0, 0, 1>(cctkGH,
                          [&](const Loop::PointDesc &p) { nancheck(jz, p); });

  // Loop::loop_all<0, 0, 0>(cctkGH,
  //                         [&](const Loop::PointDesc &p) { nancheck(dive, p);
  //                         });

  // Loop::loop_all<1, 1, 1>(cctkGH,
  //                         [&](const Loop::PointDesc &p) { nancheck(divb, p);
  //                         });
}
#endif

} // namespace MaxwellToyCarpetX
