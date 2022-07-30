#include <fixmath.hxx>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cmath>

namespace AsterX {
using namespace std;
using namespace Loop;

// Computes the determinant of spatial metric
template <typename T>
CCTK_DEVICE CCTK_HOST const T calc_detg(const T &gxx, const T &gxy,
                                        const T &gxz, const T &gyy,
                                        const T &gyz, const T &gzz) {
  return -gxz * gxz * gyy + 2.0 * gxy * gxz * gyz - gxx * gyz * gyz -
         gxy * gxy * gzz + gxx * gyy * gzz;
}

// Computes the upper spatial metric components
template <typename T>
CCTK_DEVICE CCTK_HOST array<T, 6>
calc_upperg(const T &gxx, const T &gxy, const T &gxz, const T &gyy,
            const T &gyz, const T &gzz, const T &detg) {
  return {
      (gyy * gzz - gyz * gyz) / detg, // uxx
      (gxz * gyz - gxy * gzz) / detg, // uxy
      (gxy * gyz - gxz * gyy) / detg, // uxz
      (gxx * gzz - gxz * gxz) / detg, // uyy
      (gxy * gxz - gyz * gxx) / detg, // uyz
      (gxx * gyy - gxy * gxy) / detg  // uzz
  };
}

// FD2: vertex centered input, vertex centered output, oneside stencil
template <typename T>
CCTK_DEVICE CCTK_HOST T calc_fd2_v2v_oneside(const GF3D2<const T> &gf,
                                             const PointDesc &p, const int dir,
                                             const int sign) {
  constexpr auto DI = PointDesc::DI;
  return -sign *
         (gf(p.I + 2.0 * sign * DI[dir]) - 4.0 * gf(p.I + sign * DI[dir]) +
          3.0 * gf(p.I)) *
         (0.5 / p.DX[dir]);
}

// FD2: cell centered input, cell centered output
template <typename T>
CCTK_DEVICE CCTK_HOST T calc_fd2_c2c(const GF3D2<const T> &gf,
                                     const PointDesc &p, const int dir) {
  constexpr auto DI = PointDesc::DI;
  return (0.5 / p.DX[dir]) * (gf(p.I + DI[dir]) - gf(p.I - DI[dir]));
}

// FD2: vertex centered input, edge centered output
template <typename T>
CCTK_DEVICE CCTK_HOST T calc_fd2_v2e(const GF3D2<const T> &gf,
                                     const PointDesc &p, const int dir) {
  constexpr auto DI = PointDesc::DI;
  return (gf(p.I + DI[dir]) - gf(p.I)) / p.DX[dir];
}

// FD2: vertex centered input, cell centered output
template <typename T>
CCTK_DEVICE CCTK_HOST T calc_fd2_v2c(const GF3D2<const T> &gf,
                                     const PointDesc &p, int dir) {
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
CCTK_DEVICE CCTK_HOST T calc_fd4_c2c(const GF3D2<const T> &gf,
                                     const PointDesc &p, const int dir) {
  constexpr auto DI = PointDesc::DI;
  return (1.0 / (12.0 * p.DX[dir])) *
         (-gf(p.I + 2 * DI[dir]) + 8.0 * gf(p.I + DI[dir]) -
          8.0 * gf(p.I - DI[dir]) + gf(p.I - 2 * DI[dir]));
}

// FD4: vertex centered input, cell centered output
template <typename T>
CCTK_DEVICE CCTK_HOST T calc_fd4_v2c(const GF3D2<const T> &gf,
                                     const PointDesc &p, int dir) {
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

// Computes the components of a lower indice quantity
template <typename T>
CCTK_DEVICE CCTK_HOST const array<T, 3>
calc_vlow(const array<T, 3> &v_up, const T &gxx, const T &gxy, const T &gxz,
          const T &gyy, const T &gyz, const T &gzz) {
  return {
      gxx * v_up[0] + gxy * v_up[1] + gxz * v_up[2], // vlowx
      gxy * v_up[0] + gyy * v_up[1] + gyz * v_up[2], // vlowy
      gxz * v_up[0] + gyz * v_up[1] + gzz * v_up[2]  // vlowz
  };
}

// Computes the Lorentz factor
template <typename T>
CCTK_DEVICE CCTK_HOST const T calc_wlor(const array<T, 3> &v_low,
                                        const array<T, 3> &v_up) {
  return 1.0 / sqrt(1.0 - (v_low[0] * v_up[0] + v_low[1] * v_up[1] +
                           v_low[2] * v_up[2]));
}

// Second-order average of vertex-centered grid functions to cell center
template <typename T>
CCTK_DEVICE CCTK_HOST T calc_avg_v2c(const GF3D2<const T> &gf,
                                     const PointDesc &p) {
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
CCTK_DEVICE CCTK_HOST T calc_avg_e2v(const GF3D2<const T> &gf,
                                     const PointDesc &p, const int dir) {
  constexpr auto DI = PointDesc::DI;
  T gf_avg = 0.0;

  for (int di = -1; di < 1; ++di) {
    gf_avg += gf(p.I + DI[dir] * di);
  }
  return gf_avg / 2.0;
}

// Second-order average of edge-centered grid functions (along dir) to cell
// center
template <typename T>
CCTK_DEVICE CCTK_HOST T calc_avg_v2c(const GF3D2<const T> &gf,
                                     const PointDesc &p, const int dir) {
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

} // namespace AsterX
