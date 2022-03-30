#ifndef MULTIPATCH_CAKE_HXX
#define MULTIPATCH_CAKE_HXX

#include <fixmath.hxx>
#include <cctk.h>

#include <loop.hxx>
#include <mat.hxx>
#include <vec.hxx>
#include <vect.hxx>

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
 * Stores a 3-matrx with all indices down
 */
using smat_d = smat<CCTK_REAL, dim, DN, DN>;

/**
 * Stores a Jacobian matrix of 3 dimentions
 */
using jac_t = vec<vec<CCTK_REAL, dim, DN>, dim, UP>;

/**
 * Stores a the derivatives of a Jacobian matrix of 3 dimentions
 */
using djac_t = vec<smat<CCTK_REAL, dim, DN, DN>, dim, UP>;

/**
 * Precondition assertion macro. If the precondition fails, the code is aborted.
 *
 * @param predicate The predicate to test.
 * @param msg The message to display while aborting the code.
 */
inline void expects(bool predicate, const char *msg) {
  if (!predicate) {
    CCTK_ERROR(msg);
  }
}

/**
 * Computes the positive integer power of a number at compile time. Comptible
 * with Mathematica's Power function
 *
 * @param n The exponent.
 * @param x The base.
 * @return The n-th power of x
 */
template <typename T> static inline constexpr T Power(T x, unsigned n) {
  return (n == 0) ? T(1) : x * Power(x, n - 1);
}

/**
 * Computes the integer power of a number at compile time. Comptible with
 * Mathematica's Power function
 *
 * @param n The exponent.
 * @param x The base.
 * @return The n-th power of x
 */
template <typename T> static inline constexpr T Power(T x, int n) {
  return (n < 0) ? T(1) / Power(x, unsigned(-n)) : Power(x, unsigned(n));
}

/**
 * Compatibility wrapper for replacing Mathematica's Sqrt function.
 *
 * @param x The radicand.
 */
template <typename T> static inline T Sqrt(T x) {
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
  else
    return "exterior";
}

/**
 * Determine which patch piece owns a global coordinate triplet
 *
 * @param pt The patch data
 * @param global_vars The values of the local global (x, y, z)
 */
patch_piece get_owner_patch(const PatchTransformations &pt,
                            const svec_u &global_vars);

} // namespace Cake
} // namespace MultiPatch

#endif // MULTIPATCH_CAKE_HXX