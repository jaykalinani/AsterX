#ifndef UTILS_HXX
#define UTILS_HXX

#include <fixmath.hxx>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <mat.hxx>
#include <vec.hxx>
#include <sum.hxx>
#include <simd.hxx>

#include <algorithm>
#include <array>
#include <cmath>

namespace AsterX {
using namespace std;
using namespace Loop;
using namespace Arith;

// Computes the contraction of smat and vec
template <typename T, int D>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<T, D>
calc_contraction(const smat<T, D> &g, const vec<T, D> &v) {
  return [&](int i) ARITH_INLINE {
    return sum<D>([&](int j) ARITH_INLINE { return g(i, j) * v(j); });
  };
}

template <typename T, int D>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
calc_contraction(const vec<T, D> &v_up, const vec<T, D> &v_dn) {
  return sum<D>([&](int i) ARITH_INLINE { return v_up(i) * v_dn(i); });
}

template <typename T, int D>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
calc_contraction(const smat<T, D> &g_dn, const smat<T, D> &T_up) {
  return sum_symm<D>([&](int i, int j)
                         ARITH_INLINE { return g_dn(i, j) * T_up(i, j); });
}

// Computes the Lorentz factor
template <typename T, int D>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_wlorentz(const vec<T, D> &v_up, const vec<T, D> &v_dn) {
  return 1.0 / sqrt(1.0 - calc_contraction(v_up, v_dn));
}

template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T pow2(T x) {
  return x * x;
}

// FD2: vertex centered input, vertex centered output, oneside stencil
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_fd2_v2v_oneside(const GF3D2<const T> &gf, const PointDesc &p,
                     const int dir, const int sign) {
  constexpr auto DI = PointDesc::DI;
  return -sign *
         (gf(p.I + 2 * sign * DI[dir]) - 4.0 * gf(p.I + sign * DI[dir]) +
          3.0 * gf(p.I)) *
         (0.5 / p.DX[dir]);
}

// FD2: cell centered input, cell centered output
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_fd2_c2c(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  constexpr auto DI = PointDesc::DI;
  return (0.5 / p.DX[dir]) * (gf(p.I + DI[dir]) - gf(p.I - DI[dir]));
}

// FD2: vertex centered input, edge centered output
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_fd2_v2e(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  constexpr auto DI = PointDesc::DI;
  return (gf(p.I + DI[dir]) - gf(p.I)) / p.DX[dir];
}

// FD2: vertex centered input, cell centered output
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_fd2_v2c(const GF3D2<const T> &gf, const PointDesc &p, int dir) {
  constexpr auto DI = PointDesc::DI;
  T dgf1, dgf2, dgf3, dgf4;
  int dir1, dir2;
  if (dir == 0) {
    dir1 = 1;
    dir2 = 2;
  } else if (dir == 1) {
    dir1 = 0;
    dir2 = 2;
  } else {
    dir1 = 0;
    dir2 = 1;
  }

  dgf1 = (gf(p.I + DI[dir]) - gf(p.I)) / p.DX[dir];
  dgf2 = (gf(p.I + DI[dir1] + DI[dir]) - gf(p.I + DI[dir1])) / p.DX[dir];
  dgf3 = (gf(p.I + DI[dir2] + DI[dir]) - gf(p.I + DI[dir2])) / p.DX[dir];
  dgf4 = (gf(p.I + DI[dir1] + DI[dir2] + DI[dir]) -
          gf(p.I + DI[dir1] + DI[dir2])) /
         p.DX[dir];

  return 0.25 * (dgf1 + dgf2 + dgf3 + dgf4);
}

// FD4: cell centered input, cell centered output
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_fd4_c2c(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  constexpr auto DI = PointDesc::DI;
  return (1.0 / (12.0 * p.DX[dir])) *
         (-gf(p.I + 2 * DI[dir]) + 8.0 * gf(p.I + DI[dir]) -
          8.0 * gf(p.I - DI[dir]) + gf(p.I - 2 * DI[dir]));
}

// FD4: vertex centered input, cell centered output
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_fd4_v2c(const GF3D2<const T> &gf, const PointDesc &p, int dir) {
  constexpr auto DI = PointDesc::DI;
  T dgf1, dgf2, dgf3, dgf4;
  T dgf5, dgf6, dgf7, dgf8;
  int dir1, dir2;

  if (dir == 0) {
    dir1 = 1;
    dir2 = 2;
  } else if (dir == 1) {
    dir1 = 0;
    dir2 = 2;
  } else {
    dir1 = 0;
    dir2 = 1;
  }

  dgf1 =
      ((1.0 / 24.0) * gf(p.I + 2 * DI[dir]) - (9.0 / 8.0) * gf(p.I + DI[dir]) +
       (9.0 / 8.0) * gf(p.I) - (1.0 / 24.0) * gf(p.I - DI[dir])) /
      p.DX[dir];
  dgf2 = ((1.0 / 24.0) * gf(p.I + DI[dir1] + 2 * DI[dir]) -
          (9.0 / 8.0) * gf(p.I + DI[dir1] + DI[dir]) +
          (9.0 / 8.0) * gf(p.I + DI[dir1]) -
          (1.0 / 24.0) * gf(p.I + DI[dir1] - DI[dir])) /
         p.DX[dir];
  dgf3 = ((1.0 / 24.0) * gf(p.I + DI[dir2] + 2 * DI[dir]) -
          (9.0 / 8.0) * gf(p.I + DI[dir2] + DI[dir]) +
          (9.0 / 8.0) * gf(p.I + DI[dir2]) -
          (1.0 / 24.0) * gf(p.I + DI[dir2] - DI[dir])) /
         p.DX[dir];
  dgf4 = ((1.0 / 24.0) * gf(p.I + DI[dir1] + DI[dir2] + 2 * DI[dir]) -
          (9.0 / 8.0) * gf(p.I + DI[dir1] + DI[dir2] + DI[dir]) +
          (9.0 / 8.0) * gf(p.I + DI[dir1] + DI[dir2]) -
          (1.0 / 24.0) * gf(p.I + DI[dir1] + DI[dir2] - DI[dir])) /
         p.DX[dir];

  dgf5 = ((1.0 / 24.0) * gf(p.I + 2 * DI[dir1] + 2 * DI[dir2] + 2 * DI[dir]) -
          (9.0 / 8.0) * gf(p.I + 2 * DI[dir1] + 2 * DI[dir2] + DI[dir]) +
          (9.0 / 8.0) * gf(p.I + 2 * DI[dir1] + 2 * DI[dir2]) -
          (1.0 / 24.0) * gf(p.I + 2 * DI[dir1] + 2 * DI[dir2] - DI[dir])) /
         p.DX[dir];

  dgf6 = ((1.0 / 24.0) * gf(p.I + 2 * DI[dir1] - DI[dir2] + 2 * DI[dir]) -
          (9.0 / 8.0) * gf(p.I + 2 * DI[dir1] - DI[dir2] + DI[dir]) +
          (9.0 / 8.0) * gf(p.I + 2 * DI[dir1] - DI[dir2]) -
          (1.0 / 24.0) * gf(p.I + 2 * DI[dir1] - DI[dir2] - DI[dir])) /
         p.DX[dir];

  dgf7 = ((1.0 / 24.0) * gf(p.I - DI[dir1] + 2 * DI[dir2] + 2 * DI[dir]) -
          (9.0 / 8.0) * gf(p.I - DI[dir1] + 2 * DI[dir2] + DI[dir]) +
          (9.0 / 8.0) * gf(p.I - DI[dir1] + 2 * DI[dir2]) -
          (1.0 / 24.0) * gf(p.I - DI[dir1] + 2 * DI[dir2] - DI[dir])) /
         p.DX[dir];

  dgf8 = ((1.0 / 24.0) * gf(p.I - DI[dir1] - DI[dir2] + 2 * DI[dir]) -
          (9.0 / 8.0) * gf(p.I - DI[dir1] - DI[dir2] + DI[dir]) +
          (9.0 / 8.0) * gf(p.I - DI[dir1] - DI[dir2]) -
          (1.0 / 24.0) * gf(p.I - DI[dir1] - DI[dir2] - DI[dir])) /
         p.DX[dir];

  return (27.0 / 14.0) * (dgf1 + dgf2 + dgf3 + dgf4) -
         (10.0 / 7.0) * (dgf5 + dgf6 + dgf7 + dgf8);
}

// Second-order average of vertex-centered grid functions to cell center
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_avg_v2c(const GF3D2<const T> &gf, const PointDesc &p) {
  constexpr auto DI = PointDesc::DI;
  T gf_avg = 0.0;

  for (int dk = 0; dk < 2; ++dk) {
    for (int dj = 0; dj < 2; ++dj) {
      for (int di = 0; di < 2; ++di) {
        gf_avg += gf(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
      }
    }
  }
  return gf_avg / 8.0;
}

// Second-order average of edge-centered grid functions to vertex-centered
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_avg_e2v(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  constexpr auto DI = PointDesc::DI;
  T gf_avg = 0.0;

  for (int di = 0; di < 2; ++di) {
    gf_avg += gf(p.I - DI[dir] * di);
  }
  return gf_avg / 2.0;
}

// Second-order average of edge-centered grid functions (along dir) to cell
// center
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_avg_e2c(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  constexpr auto DI = PointDesc::DI;
  T gf_avg = 0.0;

  for (int dk = 0; dk < (dir == 2 ? 1 : 2); ++dk) {
    for (int dj = 0; dj < (dir == 1 ? 1 : 2); ++dj) {
      for (int di = 0; di < (dir == 0 ? 1 : 2); ++di) {
        gf_avg += gf(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
      }
    }
  }
  return gf_avg / 4.0;
}

// Second-order average of vertex-centered grid functionsto face
// center (perp to dir)
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_avg_v2f(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  constexpr auto DI = PointDesc::DI;
  T gf_avg = 0.0;

  for (int dk = 0; dk < (dir == 2 ? 1 : 2); ++dk) {
    for (int dj = 0; dj < (dir == 1 ? 1 : 2); ++dj) {
      for (int di = 0; di < (dir == 0 ? 1 : 2); ++di) {
        gf_avg += gf(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
      }
    }
  }
  return gf_avg / 4.0;
}

// Second-order average of cell-centered grid functions to edge center
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_avg_c2e(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  constexpr auto DI = PointDesc::DI;
  T gf_avg = 0.0;

  for (int dk = 0; dk < (dir == 2 ? 1 : 2); ++dk) {
    for (int dj = 0; dj < (dir == 1 ? 1 : 2); ++dj) {
      for (int di = 0; di < (dir == 0 ? 1 : 2); ++di) {
        gf_avg += gf(p.I - DI[0] * di - DI[1] * dj - DI[2] * dk);
      }
    }
  }
  return gf_avg / 4.0;
}

} // namespace AsterX

#endif // #ifndef UTILS_HXX
