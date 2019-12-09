#include "dual.hxx"

#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>
#include <complex>
#include <ostream>

namespace Maxwell {
using namespace std;

template <typename T> T pow2(T x) { return x * x; }

////////////////////////////////////////////////////////////////////////////////

template <typename T> struct potential { T phi, ax, ay, az; };

// template <typename F, typename T>
// potential<T> calc_dtp(const F &f, T t, T x, T y, T z, T dt) {
//   auto fm = f(t, x, y, z);
//   auto fp = f(t + dt, x, y, z);
//   return {
//       .phi = (fp.phi - fm.phi) / dt,
//       .ax = (fp.ax - fm.ax) / dt,
//       .ay = (fp.ay - fm.ay) / dt,
//       .az = (fp.az - fm.az) / dt,
//   };
// }

template <typename F, typename T>
potential<T> calc_dxp(const F &f, T t, T x, T y, T z, T dx) {
  auto fm = f(t, x, y, z);
  auto fp = f(t, x + dx, y, z);
  return {
      .phi = (fp.phi - fm.phi) / dx,
      .ax = (fp.ax - fm.ax) / dx,
      .ay = (fp.ay - fm.ay) / dx,
      .az = (fp.az - fm.az) / dx,
  };
}

template <typename F, typename T>
potential<T> calc_dyp(const F &f, T t, T x, T y, T z, T dy) {
  auto fm = f(t, x, y, z);
  auto fp = f(t, x, y + dy, z);
  return {
      .phi = (fp.phi - fm.phi) / dy,
      .ax = (fp.ax - fm.ax) / dy,
      .ay = (fp.ay - fm.ay) / dy,
      .az = (fp.az - fm.az) / dy,
  };
}

template <typename F, typename T>
potential<T> calc_dzp(const F &f, T t, T x, T y, T z, T dz) {
  auto fm = f(t, x, y, z);
  auto fp = f(t, x, y, z + dz);
  return {
      .phi = (fp.phi - fm.phi) / dz,
      .ax = (fp.ax - fm.ax) / dz,
      .ay = (fp.ay - fm.ay) / dz,
      .az = (fp.az - fm.az) / dz,
  };
}

// template <typename F, typename T>
// potential<T> calc_dtm(const F &f, T t, T x, T y, T z, T dt) {
//   auto fm = f(t - dt, x, y, z);
//   auto fp = f(t, x, y, z);
//   return {
//       .phi = (fp.phi - fm.phi) / dt,
//       .ax = (fp.ax - fm.ax) / dt,
//       .ay = (fp.ay - fm.ay) / dt,
//       .az = (fp.az - fm.az) / dt,
//   };
// }
//
// template <typename F, typename T>
// potential<T> calc_dxm(const F &f, T t, T x, T y, T z, T dx) {
//   auto fm = f(t, x - dx, y, z);
//   auto fp = f(t, x, y, z);
//   return {
//       .phi = (fp.phi - fm.phi) / dx,
//       .ax = (fp.ax - fm.ax) / dx,
//       .ay = (fp.ay - fm.ay) / dx,
//       .az = (fp.az - fm.az) / dx,
//   };
// }
//
// template <typename F, typename T>
// potential<T> calc_dym(const F &f, T t, T x, T y, T z, T dy) {
//   auto fm = f(t, x, y - dy, z);
//   auto fp = f(t, x, y, z);
//   return {
//       .phi = (fp.phi - fm.phi) / dy,
//       .ax = (fp.ax - fm.ax) / dy,
//       .ay = (fp.ay - fm.ay) / dy,
//       .az = (fp.az - fm.az) / dy,
//   };
// }
//
// template <typename F, typename T>
// potential<T> calc_dzm(const F &f, T t, T x, T y, T z, T dz) {
//   auto fm = f(t, x, y, z - dz);
//   auto fp = f(t, x, y, z);
//   return {
//       .phi = (fp.phi - fm.phi) / dz,
//       .ax = (fp.ax - fm.ax) / dz,
//       .ay = (fp.ay - fm.ay) / dz,
//       .az = (fp.az - fm.az) / dz,
//   };
// }

template <typename F, typename T>
potential<T> calc_dt(const F &f, T t, T x, T y, T z) {
  auto ff = f(dual<T>(t, 1), dual<T>(x), dual<T>(y), dual<T>(z));
  return {
      ff.phi.eps,
      ff.ax.eps,
      ff.ay.eps,
      ff.az.eps,
  };
}

// template <typename F, typename T>
// potential<T> calc_dx(const F &f, T t, T x, T y, T z) {
//   auto ff = f(dual<T>(t), dual<T>(x, 1), dual<T>(y), dual<T>(z));
//   return {
//       ff.phi.eps,
//       ff.ax.eps,
//       ff.ay.eps,
//       ff.az.eps,
//   };
// }
//
// template <typename F, typename T>
// potential<T> calc_dy(const F &f, T t, T x, T y, T z) {
//   auto ff = f(dual<T>(t), dual<T>(x), dual<T>(y, 1), dual<T>(z));
//   return {
//       ff.phi.eps,
//       ff.ax.eps,
//       ff.ay.eps,
//       ff.az.eps,
//   };
// }
//
// template <typename F, typename T>
// potential<T> calc_dz(const F &f, T t, T x, T y, T z) {
//   auto ff = f(dual<T>(t), dual<T>(x), dual<T>(y), dual<T>(z, 1));
//   return {
//       ff.phi.eps,
//       ff.ax.eps,
//       ff.ay.eps,
//       ff.az.eps,
//   };
// }

////////////////////////////////////////////////////////////////////////////////

// Plane wave
template <typename T>
potential<T> plane_wave(T t, T x, T y, T z, T dx, T dy, T dz) {
  DECLARE_CCTK_PARAMETERS;
  // wave number
  T kx = M_PI * spatial_frequency_x;
  T ky = M_PI * spatial_frequency_y;
  T kz = M_PI * spatial_frequency_z;
  // choose omega to ensure Lorenz gauge
  // continuum
  T omega = sqrt(pow2(kx) + pow2(ky) + pow2(kz));
  // note: we don't have a discrete choice for omega
  T hx = amplitude_x;
  T hy = amplitude_y;
  T hz = amplitude_z;
  // choose ht to ensure div E = 0
  // continuum
  // T ht =
  //     omega * (hx * kx + hy * ky + hz * kz) / (pow2(kx) + pow2(ky) +
  //     pow2(kz));
  // discrete
  typedef complex<T> CT;
  CT ht =
      omega *
      (hx / dx * CT(sin(dx * kx), 2 * pow2(sin(dx * kx / 2))) +
       hy / dy * CT(sin(dy * ky), 2 * pow2(sin(dy * ky / 2))) +
       hz / dz * CT(sin(dz * kz), 2 * pow2(sin(dz * kz / 2)))) /
      CT(4 * ((pow2(sin(dx * kx / 2) / dx)) + (pow2(sin(dy * ky / 2) / dy)) +
              (pow2(sin(dz * kz / 2) / dz))));
  CT u = CT(cos(omega * t - kx * x - ky * y - kz * z),
            sin(omega * t - kx * x - ky * y - kz * z));
  return {
      .phi = real(ht * u),
      .ax = real(hx * u),
      .ay = real(hy * u),
      .az = real(hz * u),
  };
}

// Plane wave with Gaussian profile
template <typename T> potential<T> gaussian_wave(T t, T x, T y, T z) {
  DECLARE_CCTK_PARAMETERS;
  T kx = M_PI * spatial_frequency_x;
  T ky = M_PI * spatial_frequency_y;
  T kz = M_PI * spatial_frequency_z;
  T omega = sqrt(pow2(kx) + pow2(ky) + pow2(kz));
  T hx = amplitude_x;
  T hy = amplitude_y;
  T hz = amplitude_z;
  T ht =
      omega * (hx * kx + hy * ky + hz * kz) / (pow2(kx) + pow2(ky) + pow2(kz));
  T u = exp(-pow2(sin(omega * t - kx * x - ky * y - kz * z) / width) / 2);
  return {
      .phi = ht * u,
      .ax = hx * u,
      .ay = hy * u,
      .az = hz * u,
  };
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void Maxwell_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Maxwell_Initial;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;

  // const CCTK_REAL dt = CCTK_DELTA_TIME;
  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  const Loop::GF3D<CCTK_REAL, 0, 0, 0> phi_(cctkGH, phi);

  const Loop::GF3D<CCTK_REAL, 1, 0, 0> ax_(cctkGH, ax);
  const Loop::GF3D<CCTK_REAL, 0, 1, 0> ay_(cctkGH, ay);
  const Loop::GF3D<CCTK_REAL, 0, 0, 1> az_(cctkGH, az);

  const Loop::GF3D<CCTK_REAL, 1, 0, 0> ex_(cctkGH, ex);
  const Loop::GF3D<CCTK_REAL, 0, 1, 0> ey_(cctkGH, ey);
  const Loop::GF3D<CCTK_REAL, 0, 0, 1> ez_(cctkGH, ez);

  const Loop::GF3D<CCTK_REAL, 0, 1, 1> byz_(cctkGH, byz);
  const Loop::GF3D<CCTK_REAL, 1, 0, 1> bzx_(cctkGH, bzx);
  const Loop::GF3D<CCTK_REAL, 1, 1, 0> bxy_(cctkGH, bxy);

  const auto loop_setup{[&](const auto &f4) {
    const auto f{[&](const auto &p) { return f4(t, p.x, p.y, p.z); }};
    const auto dtf{[&](const auto &p) {
      return calc_dt(
          [&](auto t, auto x, auto y, auto z) { return f4(t, x, y, z); }, t,
          p.x, p.y, p.z);
    }};
    const auto dxpf{[&](const auto &p) {
      return calc_dxp(
          [&](auto t, auto x, auto y, auto z) { return f4(t, x, y, z); }, t,
          p.x, p.y, p.z, dx);
    }};
    const auto dypf{[&](const auto &p) {
      return calc_dyp(
          [&](auto t, auto x, auto y, auto z) { return f4(t, x, y, z); }, t,
          p.x, p.y, p.z, dy);
    }};
    const auto dzpf{[&](const auto &p) {
      return calc_dzp(
          [&](auto t, auto x, auto y, auto z) { return f4(t, x, y, z); }, t,
          p.x, p.y, p.z, dz);
    }};

    Loop::loop_int<0, 0, 0>(
        cctkGH, [&](const Loop::PointDesc &p) { phi_(p.I) = f(p).phi; });

    Loop::loop_int<1, 0, 0>(
        cctkGH, [&](const Loop::PointDesc &p) { ax_(p.I) = f(p).ax; });
    Loop::loop_int<0, 1, 0>(
        cctkGH, [&](const Loop::PointDesc &p) { ay_(p.I) = f(p).ay; });
    Loop::loop_int<0, 0, 1>(
        cctkGH, [&](const Loop::PointDesc &p) { az_(p.I) = f(p).az; });

    Loop::loop_int<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      ex_(p.I) = -dxpf(p).phi - dtf(p).ax;
    });
    Loop::loop_int<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      ey_(p.I) = -dypf(p).phi - dtf(p).ay;
    });
    Loop::loop_int<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      ez_(p.I) = -dzpf(p).phi - dtf(p).az;
    });

    Loop::loop_int<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      byz_(p.I) = dypf(p).az - dzpf(p).ay;
    });
    Loop::loop_int<1, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      bzx_(p.I) = dzpf(p).az - dxpf(p).az;
    });
    Loop::loop_int<1, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      bxy_(p.I) = dxpf(p).ay - dypf(p).ax;
    });
  }};

  if (CCTK_EQUALS(setup, "plane wave")) {
    loop_setup([&](auto t, auto x, auto y, auto z) {
      typedef decltype(t) T;
      return plane_wave(t, x, y, z, T(dx), T(dy), T(dz));
    });
    // } else if (CCTK_EQUALS(setup, "Gaussian wave")) {
    //   loop_setup([&](auto t, auto x, auto y, auto z) {
    //     return gaussian_wave(t, x, y, z);
    //   });
  } else {
    assert(0);
  }
}

} // namespace Maxwell
