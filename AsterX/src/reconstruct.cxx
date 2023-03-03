#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

#include "utils.hxx"
#include "reconstruct.hxx"

namespace AsterX {
using namespace std;
using namespace Loop;

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T minmod(const T &x,
                                                                   const T &y) {
  if (signbit(x) != signbit(y))
    return T(0);
  if (fabs(x) < fabs(y))
    return x;
  else
    return y;
}

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T
monocentral(const T &x, const T &y) {
  if (sgn(x) != sgn(y))
    return 0;
  else
    return sgn(x) * min(2 * fabs(x), min(2 * fabs(y), fabs(x + y) / 2));
}

CCTK_DEVICE array<CCTK_REAL, 2>
reconstruct(const GF3D2<const CCTK_REAL> &gf_var, const PointDesc &p,
            reconstruction_t reconstruction, int dir) {
  constexpr auto DI = PointDesc::DI;
  // Neighbouring "plus" and "minus" cell indices
  const auto Immm = p.I - 3 * DI[dir];
  const auto Imm = p.I - 2 * DI[dir];
  const auto Im = p.I - DI[dir];
  const auto Ip = p.I;
  const auto Ipp = p.I + DI[dir];
  const auto Ippp = p.I + 2 * DI[dir];

  switch (reconstruction) {

  case reconstruction_t::Godunov: {
    CCTK_REAL var_m = gf_var(Im);
    CCTK_REAL var_p = gf_var(Ip);
    return array<CCTK_REAL, 2>{var_m, var_p};
  }

  case reconstruction_t::minmod: {
    // reconstructs values of Im and Ip at the common face between these
    // two cells
    CCTK_REAL var_slope_p = gf_var(Ipp) - gf_var(Ip);
    CCTK_REAL var_slope_c = gf_var(Ip) - gf_var(Im);
    CCTK_REAL var_slope_m = gf_var(Im) - gf_var(Imm);
    // reconstructed Im on its "plus/right" side
    CCTK_REAL var_m = gf_var(Im) + minmod(var_slope_c, var_slope_m) / 2;
    // reconstructed Ip on its "minus/left" side
    CCTK_REAL var_p = gf_var(Ip) - minmod(var_slope_p, var_slope_c) / 2;
    return array<CCTK_REAL, 2>{var_m, var_p};
  }

  case reconstruction_t::monocentral: {
    // reconstructs values of Im and Ip at the common face between these
    // two cells
    // reconstructed Im on its "plus/right" side
    CCTK_REAL var_slope_p = gf_var(Ip) - gf_var(Im);
    CCTK_REAL var_slope_m = gf_var(Im) - gf_var(Imm);
    CCTK_REAL var_m = gf_var(Im) + monocentral(var_slope_p, var_slope_m) / 2;
    // reconstructed Ip on its "minus/left" side
    var_slope_p = gf_var(Ipp) - gf_var(Ip);
    var_slope_m = gf_var(Ip) - gf_var(Im);
    CCTK_REAL var_p = gf_var(Ip) - monocentral(var_slope_p, var_slope_m) / 2;
    return array<CCTK_REAL, 2>{var_m, var_p};
  }

  case reconstruction_t::ppm: {
    // Usually recon methods return the left and right states of a given
    // cell face, ie (i-1/2-eps) and (i-1/2+eps).
    // But ppm return the states at the left and right
    // faces of a given cell, ie (i-1/2) and (i+1/2).
    // So, to get left and right states from (i-1/2),
    // we first apply ppm to cell i and keep left face (i-1/2+eps)
    // and then apply ppm to cell i-1 and keep right face (i-1/2-eps)

    // Start calculating left (i-1/2) and right (i+1/2) faces unique
    // values: Eq. (A1) in https://arxiv.org/pdf/astro-ph/0503420.pdf with
    // 1/8 --> 1/6 Equiv. to Eq. (1.6) in C&W (1984)
    CCTK_REAL var_slope_p = gf_var(Ip) - gf_var(Im);
    CCTK_REAL var_slope_m = gf_var(Im) - gf_var(Imm);
    CCTK_REAL grad_m = monocentral(var_slope_p, var_slope_m);
    var_slope_p = gf_var(Ipp) - gf_var(Ip);
    var_slope_m = gf_var(Ip) - gf_var(Im);
    CCTK_REAL grad_p = monocentral(var_slope_p, var_slope_m);
    CCTK_REAL left_face =
        (gf_var(Im) + gf_var(Ip)) / 2.0 + (grad_m - grad_p) / 6.0;

    var_slope_p = gf_var(Ipp) - gf_var(Ip);
    var_slope_m = gf_var(Ip) - gf_var(Im);
    grad_m = monocentral(var_slope_p, var_slope_m);
    var_slope_p = gf_var(Ippp) - gf_var(Ipp);
    var_slope_m = gf_var(Ipp) - gf_var(Ip);
    grad_p = monocentral(var_slope_p, var_slope_m);
    CCTK_REAL right_face =
        (gf_var(Ip) + gf_var(Ipp)) / 2.0 + (grad_m - grad_p) / 6.0;

    // Now apply conditions in Eq. (1.11) of C&W (1984)
    CCTK_REAL qa = (right_face - gf_var(Ip)) * (gf_var(Ip) - left_face);
    CCTK_REAL qd = (right_face - left_face);
    CCTK_REAL qe = 6.0 * (gf_var(Ip) - (left_face + right_face) / 2.0);
    if (qa <= 0.) {
      left_face = gf_var(Ip);
      right_face = gf_var(Ip);
    }
    if (qd * (qd - qe) < 0.0)
      left_face = 3.0 * gf_var(Ip) - 2.0 * right_face;
    if (qd * (qd + qe) < 0.0)
      right_face = 3.0 * gf_var(Ip) - 2.0 * left_face;

    // Keep left value of cell i as i-1/2+eps
    CCTK_REAL var_p = left_face;

    // Start calculating left (i-1-1/2) and right (i-1+1/2) faces unique
    // values: Eq. (A1) in https://arxiv.org/pdf/astro-ph/0503420.pdf with
    // 1/8 --> 1/6 Equiv. to Eq. (1.6) in C&W (1984)
    var_slope_p = gf_var(Im) - gf_var(Imm);
    var_slope_m = gf_var(Imm) - gf_var(Immm);
    grad_m = monocentral(var_slope_p, var_slope_m);
    var_slope_p = gf_var(Ip) - gf_var(Im);
    var_slope_m = gf_var(Im) - gf_var(Imm);
    grad_p = monocentral(var_slope_p, var_slope_m);
    left_face = (gf_var(Imm) + gf_var(Im)) / 2.0 + (grad_m - grad_p) / 6.0;

    var_slope_p = gf_var(Ip) - gf_var(Im);
    var_slope_m = gf_var(Im) - gf_var(Imm);
    grad_m = monocentral(var_slope_p, var_slope_m);
    var_slope_p = gf_var(Ipp) - gf_var(Ip);
    var_slope_m = gf_var(Ip) - gf_var(Im);
    grad_p = monocentral(var_slope_p, var_slope_m);
    right_face = (gf_var(Im) + gf_var(Ip)) / 2.0 + (grad_m - grad_p) / 6.0;

    // Now apply conditions in Eq. (1.11) of C&W (1984)
    qa = (right_face - gf_var(Im)) * (gf_var(Im) - left_face);
    qd = (right_face - left_face);
    qe = 6.0 * (gf_var(Im) - (left_face + right_face) / 2.0);
    if (qa <= 0.) {
      left_face = gf_var(Im);
      right_face = gf_var(Im);
    }
    if (qd * (qd - qe) < 0.0)
      left_face = 3.0 * gf_var(Im) - 2.0 * right_face;
    if (qd * (qd + qe) < 0.0)
      right_face = 3.0 * gf_var(Im) - 2.0 * left_face;

    // Keep right value of cell i-1 as i-1/2-eps
    CCTK_REAL var_m = right_face;

    return array<CCTK_REAL, 2>{var_m, var_p};
  }

  default:
    assert(0);
  }
}

} // namespace AsterX
