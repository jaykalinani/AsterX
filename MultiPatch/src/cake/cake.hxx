#ifndef MULTIPATCH_CAKE_HXX
#define MULTIPATCH_CAKE_HXX

#include "multipatch.hxx"
#include "tests.hxx"

#include <string>
#include <type_traits>

namespace MultiPatch {
namespace Cake {

/**
 * Stores a 3-vector with all indices down
 */
using svec_d = vec<CCTK_REAL, dim, DN>;

/**
 * Stores a 3-vector with all indices up
 */
using svec_u = vec<CCTK_REAL, dim, UP>;

/**
 * Stores a 3-matrix with all indices down
 */
using smat_d = smat<CCTK_REAL, dim, DN, DN>;

/**
 * Stores a Jacobian matrix of 3 dimensions
 */
using jac_t = vec<vec<CCTK_REAL, dim, DN>, dim, UP>;

/**
 * Stores a the derivatives of a Jacobian matrix of 3 dimensions
 */
using djac_t = vec<smat<CCTK_REAL, dim, DN, DN>, dim, UP>;

/**
 * Precondition assertion. If the precondition fails, the code is aborted.
 *
 * @param predicate The predicate to test.
 * @param msg The message to display while aborting the code.
 */
CCTK_DEVICE CCTK_HOST inline void expects(bool predicate, const char *msg) {
  if (!predicate) {
#ifndef __CUDACC__
    CCTK_ERROR(msg);
#else
    assert(0);
#endif
  }
}

/**
 * Computes the power of a number at compile time. Comptible
 * with Mathematica's Power function
 *
 * @param n The exponent.
 * @param x The base.
 * @return The n-th power of x
 */
template <typename T>
CCTK_DEVICE CCTK_HOST static inline constexpr T Power(T x, int n) {
  if (n == 0) {
    return T{1};
  } else if (n < 0) {
    return T{1} / Power(x, -n);
  } else {
    return x * Power(x, n - 1);
  }
}

/**
 * Compatibility wrapper for replacing Mathematica's Sqrt function.
 *
 * @param x The radicand.
 */
template <typename T> CCTK_DEVICE CCTK_HOST static inline T Sqrt(T x) {
  using std::sqrt;
  return sqrt(x);
}

/**
 * Tags for each patch piece in the cake
 */
enum class patch_piece : int {
  cartesian,

  plus_x,
  minus_x,

  plus_y,
  minus_y,

  plus_z,
  minus_z,

  inner_boundary,
  outer_boundary,

  exterior = -1
};

/**
 * Gets name from patch piece
 *
 * @tparam p A patch piece type.
 * @return A string representing the name of the piece.
 */
inline const std::string piece_name(const patch_piece &p) {
  if (p == patch_piece::cartesian)
    return "cartesian";
  else if (p == patch_piece::plus_x)
    return "plus x";
  else if (p == patch_piece::minus_x)
    return "minus x";
  else if (p == patch_piece::plus_y)
    return "plus y";
  else if (p == patch_piece::minus_y)
    return "minus y";
  else if (p == patch_piece::plus_z)
    return "plus z";
  else if (p == patch_piece::minus_z)
    return "minus z";
  else if (p == patch_piece::minus_z)
    return "interpatch boundary";
  else if (p == patch_piece::inner_boundary)
    return "inner cube boundary";
  else if (p == patch_piece::outer_boundary)
    return "outer boundary";
  else
    return "exterior";
}

/**
 * Get the patch piece that owns a global coordinate point.
 *
 * @param pt The PatchTransformations structure describing the patch system.
 * @param global_vars The global coordinate triplet to locate the owner for.
 * @return The patch piece owning the global coordinates.
 */
CCTK_DEVICE CCTK_HOST patch_piece
get_owner_patch(const PatchTransformations &pt, const svec_u &global_vars);

} // namespace Cake
} // namespace MultiPatch

#endif // MULTIPATCH_CAKE_HXX
