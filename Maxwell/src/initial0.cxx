#include "dual.hxx"

#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>
#include <complex>

namespace Maxwell {
using namespace std;

namespace {
int bitsign(bool s) { return s ? -1 : +1; }

template <typename T> T pow2(T x) { return x * x; }
template <typename T> T sinc(T x) { return x == T(0) ? T(1) : sin(x) / x; }
template <typename T> complex<T> cis(T x) { return {cos(x), sin(x)}; }
} // namespace

////////////////////////////////////////////////////////////////////////////////

template <typename T> struct potential { T phi, ax, ay, az; };

// Continuous derivative
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

// Discrete derivative
template <typename F, typename T>
potential<T> calc_dxc(const F &f, T t, T x, T y, T z, T dx) {
  auto fm = f(t, x - dx / 2, y, z);
  auto fp = f(t, x + dx / 2, y, z);
  return {
      .phi = (fp.phi - fm.phi) / dx,
      .ax = (fp.ax - fm.ax) / dx,
      .ay = (fp.ay - fm.ay) / dx,
      .az = (fp.az - fm.az) / dx,
  };
}

template <typename F, typename T>
potential<T> calc_dyc(const F &f, T t, T x, T y, T z, T dy) {
  auto fm = f(t, x, y - dy / 2, z);
  auto fp = f(t, x, y + dy / 2, z);
  return {
      .phi = (fp.phi - fm.phi) / dy,
      .ax = (fp.ax - fm.ax) / dy,
      .ay = (fp.ay - fm.ay) / dy,
      .az = (fp.az - fm.az) / dy,
  };
}

template <typename F, typename T>
potential<T> calc_dzc(const F &f, T t, T x, T y, T z, T dz) {
  auto fm = f(t, x, y, z - dz / 2);
  auto fp = f(t, x, y, z + dz / 2);
  return {
      .phi = (fp.phi - fm.phi) / dz,
      .ax = (fp.ax - fm.ax) / dz,
      .ay = (fp.ay - fm.ay) / dz,
      .az = (fp.az - fm.az) / dz,
  };
}

////////////////////////////////////////////////////////////////////////////////

// Plane wave implementation
template <typename T>
potential<complex<T> > plane_wave_impl(T t, T x, T y, T z, T dx, T dy, T dz,
                                       T kx, T ky, T kz, T hx, T hy, T hz) {
  DECLARE_CCTK_PARAMETERS;
  typedef complex<T> CT;
  // choose frequency to ensure div E = 0
  T omega = sqrt(pow2(sinc(kx * dx / 2) * kx) + pow2(sinc(ky * dy / 2) * ky) +
                 pow2(sinc(kz * dz / 2) * kz));
  // choose amplitude to ensure Lorenz gauge
  CT ht = (CT(hx * kx * sinc(kx * dx / 2)) * cis(kx * dx / 2) +
           CT(hy * ky * sinc(ky * dy / 2)) * cis(ky * dy / 2) +
           CT(hz * kz * sinc(kz * dz / 2)) * cis(kz * dz / 2)) /
          CT(omega);
  CT u = cis(omega * t - kx * x - ky * y - kz * z);
  return {
      .phi = ht * u,
      .ax = hx * u,
      .ay = hy * u,
      .az = hz * u,
  };
}

// Plane wave
template <typename T>
potential<T> plane_wave(T t, T x, T y, T z, T dx, T dy, T dz) {
  DECLARE_CCTK_PARAMETERS;
  // wave number
  T kx = M_PI * spatial_frequency_x;
  T ky = M_PI * spatial_frequency_y;
  T kz = M_PI * spatial_frequency_z;
  // amplitude
  T hx = amplitude_x;
  T hy = amplitude_y;
  T hz = amplitude_z;
  //
  auto p = plane_wave_impl(t, x, y, z, dx, dy, dz, kx, ky, kz, hx, hy, hz);
  return {
      .phi = real(p.phi),
      .ax = real(p.ax),
      .ay = real(p.ay),
      .az = real(p.az),
  };
}

// Plane wave with a triangle profile
template <typename T>
potential<T> triangle_wave(T t, T x, T y, T z, T dx, T dy, T dz) {
  DECLARE_CCTK_PARAMETERS;
  // wave number
  T kx = M_PI * spatial_frequency_x;
  T ky = M_PI * spatial_frequency_y;
  T kz = M_PI * spatial_frequency_z;
  // amplitude
  T hx = amplitude_x;
  T hy = amplitude_y;
  T hz = amplitude_z;
  //
  potential<T> p{0, 0, 0, 0};
  for (int i = 0; i < num_coefficients; ++i) {
    const int k = 2 * i + 1;
    const T kf = k;
    const T hf = bitsign(i & 1) / pow2(kf);
    const auto pk = plane_wave_impl(t, x, y, z, dx, dy, dz, kf * kx, kf * ky,
                                    kf * kz, hf * hx, hf * hy, hf * hz);
    p.phi += imag(pk.phi);
    p.ax += imag(pk.ax);
    p.ay += imag(pk.ay);
    p.az += imag(pk.az);
  }
  return p;
}

// Plane wave with Gaussian profile (NOT WORKING)
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
    const auto dxf{[&](const auto &p) {
      return calc_dxc(
          [&](auto t, auto x, auto y, auto z) { return f4(t, x, y, z); }, t,
          p.x, p.y, p.z, dx);
    }};
    const auto dyf{[&](const auto &p) {
      return calc_dyc(
          [&](auto t, auto x, auto y, auto z) { return f4(t, x, y, z); }, t,
          p.x, p.y, p.z, dy);
    }};
    const auto dzf{[&](const auto &p) {
      return calc_dzc(
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
      ex_(p.I) = -dxf(p).phi - dtf(p).ax;
    });
    Loop::loop_int<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      ey_(p.I) = -dyf(p).phi - dtf(p).ay;
    });
    Loop::loop_int<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      ez_(p.I) = -dzf(p).phi - dtf(p).az;
    });

    Loop::loop_int<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      byz_(p.I) = dyf(p).az - dzf(p).ay;
    });
    Loop::loop_int<1, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      bzx_(p.I) = dzf(p).ax - dxf(p).az;
    });
    Loop::loop_int<1, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      bxy_(p.I) = dxf(p).ay - dyf(p).ax;
    });
  }};

  if (CCTK_EQUALS(setup, "plane wave")) {
    loop_setup([&](auto t, auto x, auto y, auto z) {
      typedef decltype(t) T;
      return plane_wave(t, x, y, z, T(dx), T(dy), T(dz));
    });
  } else if (CCTK_EQUALS(setup, "triangle wave")) {
    loop_setup([&](auto t, auto x, auto y, auto z) {
      typedef decltype(t) T;
      return triangle_wave(t, x, y, z, T(dx), T(dy), T(dz));
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
