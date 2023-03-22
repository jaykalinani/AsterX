#ifndef INTERP_HXX
#define INTERP_HXX

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

template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_avg_c2v(const GF3D2<const T> &gf, const PointDesc &p) {
  constexpr auto DI = PointDesc::DI;
  T gf_avg = 0.0;

  for (int dk = 0; dk < 2; ++dk) {
    for (int dj = 0; dj < 2; ++dj) {
      for (int di = 0; di < 2; ++di) {
        gf_avg += gf(p.I - DI[0] * di - DI[1] * dj - DI[2] * dk);
      }
    }
  }
  return gf_avg / 8.0;
}

} // namespace AsterX

#endif // #ifndef INTERP_HXX
