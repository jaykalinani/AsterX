#include "multipatch.hxx"
#include "tests.hxx"
#include "cake_jacobians.hxx"
#include "cake.hxx"

namespace MultiPatch {
namespace Cake {

/**
 * Core function of the cake local -> global coordinate transformations.
 *
 * @param pt The cake patch data.
 * @param a The a local coordinate, ranging from (-1, 1)
 * @param b The b local coordiante, ranging from (-1, 1)
 * @param c the c local coordiante, ranging from (-1, 1)
 * @return The value of the core.
 */
inline CCTK_REAL local_to_global_cake_core(const PatchTransformations &pt,
                                           CCTK_REAL a, CCTK_REAL b,
                                           CCTK_REAL c) {
  using MultiPatchTests::at_boundary;
  using MultiPatchTests::within;
  using std::sqrt;

  // Ensure that -1 <= a,b,c <= 1
  expects(
      within(a, 1.0) || at_boundary(a, 1.0),
      "Invoked local_to_global_cake_core with the variable a not in the (-1,1) "
      "range.");
  expects(
      within(b, 1.0) || at_boundary(b, 1.0),
      "Invoked local_to_global_cake_core with the variable b not in the (-1,1) "
      "range.");
  expects(
      within(c, 1.0) || at_boundary(c, 1.0),
      "Invoked local_to_global_cake_core with the variable c not in the (-1,1) "
      "range.");

  const auto r0 = pt.cake_inner_boundary_radius;
  const auto r1 = pt.cake_outer_boundary_radius;

  const auto numerator = r0 * (1 - c) + r1 * (1 + c);
  const auto denominator =
      sqrt(4 + 2 * Power(a, 2) * (1 + c) + 2 * Power(b, 2) * (1 + c));

  return numerator / denominator;
}

/**
 * Derivative of the local-> global core function of the cake coordinate
 * transformation with respect to a.
 *
 * @param pt The cake patch data.
 * @param a The a local coordinate, ranging from (-1, 1)
 * @param b The b local coordiante, ranging from (-1, 1)
 * @param c the c local coordiante, ranging from (-1, 1)
 * @return The value of the core derivative with respect to a.
 */
inline CCTK_REAL local_to_global_cake_core_da(const PatchTransformations &pt,
                                              CCTK_REAL a, CCTK_REAL b,
                                              CCTK_REAL c) {
  using std::sqrt;

  const auto r0 = pt.cake_inner_boundary_radius;
  const auto r1 = pt.cake_outer_boundary_radius;

  const auto factor = 4 + 2 * Power(a, 2) * (1 + c) + 2 * Power(b, 2) * (1 + c);
  const auto denominator = sqrt(factor) * factor;
  const auto numerator = 2 * a * (1 + c) * ((1 - c) * r0 + (1 + c) * r1);

  return -numerator / denominator;
}

/**
 * Derivative of the local -> global core function of the cake coordinate
 * transformation with respect to b.
 *
 * @param pt The cake patch data.
 * @param a The a local coordinate, ranging from (-1, 1)
 * @param b The b local coordiante, ranging from (-1, 1)
 * @param c the c local coordiante, ranging from (-1, 1)
 * @return The value of the core derivative with respect to b.
 */
inline CCTK_REAL local_to_global_cake_core_db(const PatchTransformations &pt,
                                              CCTK_REAL a, CCTK_REAL b,
                                              CCTK_REAL c) {
  using std::sqrt;

  const auto r0 = pt.cake_inner_boundary_radius;
  const auto r1 = pt.cake_outer_boundary_radius;

  const auto factor = 4 + 2 * Power(a, 2) * (1 + c) + 2 * Power(b, 2) * (1 + c);
  const auto denominator = factor * sqrt(factor);
  const auto numerator = 2 * b * (1 + c) * ((1 - c) * r0 + (1 + c) * r1);

  return -numerator / denominator;
}

/**
 * Derivative of the local -> global core function of the cake coordinate
 * transformation with respect to c.
 *
 * @param pt The cake patch data.
 * @param a The a local coordinate, ranging from (-1, 1)
 * @param b The b local coordiante, ranging from (-1, 1)
 * @param c the c local coordiante, ranging from (-1, 1)
 * @return The value of the core derivative with respect to c.
 */
inline CCTK_REAL local_to_global_cake_core_dc(const PatchTransformations &pt,
                                              CCTK_REAL a, CCTK_REAL b,
                                              CCTK_REAL c) {
  using std::sqrt;

  const auto r0 = pt.cake_inner_boundary_radius;
  const auto r1 = pt.cake_outer_boundary_radius;

  const auto factor = 4 + 2 * Power(a, 2) * (1 + c) + 2 * Power(b, 2) * (1 + c);
  const auto denominator = sqrt(factor) * factor;
  const auto numerator =
      (Power(a, 2) + Power(b, 2)) * ((c - 1) * r0 - (c + 1) * r1) +
      (r1 - r0) * factor;

  return numerator / denominator;
}

/**
 * Second derivative of the local -> global core function of the cake coordinate
 * transformation with respect to a and a.
 *
 * @param pt The cake patch data.
 * @param a The a local coordinate, ranging from (-1, 1)
 * @param b The b local coordiante, ranging from (-1, 1)
 * @param c the c local coordiante, ranging from (-1, 1)
 * @return The value of the core derivative with respect to c.
 */
inline CCTK_REAL local_to_global_cake_core_da_da(const PatchTransformations &pt,
                                                 CCTK_REAL a, CCTK_REAL b,
                                                 CCTK_REAL c) {
  using std::sqrt;

  const auto r0 = pt.cake_inner_boundary_radius;
  const auto r1 = pt.cake_outer_boundary_radius;

  const auto factor_1 = 2 + Power(a, 2) * (1 + c) + Power(b, 2) * (1 + c);
  const auto factor_2 = 2 * factor_1;
  const auto denominator = sqrt(factor_1) * Power(factor_2, 3);
  const auto numerator = sqrt(2) * (1 + c) * ((c - 1) * r0 - (c + 1) * r1) *
                         (6 * Power(a, 2) * (1 + c) - Power(factor_2, 2));

  return -numerator / denominator;
}

/**
 * Second derivative of the local -> global core function of the cake coordinate
 * transformation with respect to a and b.
 *
 * @param pt The cake patch data.
 * @param a The a local coordinate, ranging from (-1, 1)
 * @param b The b local coordiante, ranging from (-1, 1)
 * @param c the c local coordiante, ranging from (-1, 1)
 * @return The value of the core derivative with respect to c.
 */
inline CCTK_REAL local_to_global_cake_core_da_db(const PatchTransformations &pt,
                                                 CCTK_REAL a, CCTK_REAL b,
                                                 CCTK_REAL c) {
  using std::sqrt;

  const auto r0 = pt.cake_inner_boundary_radius;
  const auto r1 = pt.cake_outer_boundary_radius;

  const auto factor = 4 + 2 * Power(a, 2) * (1 + c) + 2 * Power(b, 2) * (1 + c);
  const auto denominator = sqrt(factor) * Power(factor, 3);
  const auto numerator =
      12 * a * b * Power(1 + c, 2) * ((1 - c) * r0 + (1 + c) * r1);

  return numerator / denominator;
}

/**
 * Second derivative of the local -> global core function of the cake coordinate
 * transformation with respect to a and c.
 *
 * @param pt The cake patch data.
 * @param a The a local coordinate, ranging from (-1, 1)
 * @param b The b local coordiante, ranging from (-1, 1)
 * @param c the c local coordiante, ranging from (-1, 1)
 * @return The value of the core derivative with respect to c.
 */
inline CCTK_REAL local_to_global_cake_core_da_dc(const PatchTransformations &pt,
                                                 CCTK_REAL a, CCTK_REAL b,
                                                 CCTK_REAL c) {
  using std::sqrt;

  const auto r0 = pt.cake_inner_boundary_radius;
  const auto r1 = pt.cake_outer_boundary_radius;

  const auto factor_1 = 2 + Power(a, 2) * (1 + c) + Power(b, 2) * (1 + c);
  const auto factor_2 = 2 * factor_1;
  const auto denominator = sqrt(factor_1) * Power(factor_2, 3);
  const auto numerator = Sqrt(2) * a *
                         (3 * (Power(a, 2) + Power(b, 2)) * (1 + c) *
                              ((-1 + c) * r0 - (1 + c) * r1) +
                          2 * (r1 + c * (-r0 + r1)) * Power(factor_2, 2));

  return -numerator / denominator;
}

/**
 * Second derivative of the local -> global core function of the cake coordinate
 * transformation with respect to b and b.
 *
 * @param pt The cake patch data.
 * @param a The a local coordinate, ranging from (-1, 1)
 * @param b The b local coordiante, ranging from (-1, 1)
 * @param c the c local coordiante, ranging from (-1, 1)
 * @return The value of the core derivative with respect to c.
 */
inline CCTK_REAL local_to_global_cake_core_db_db(const PatchTransformations &pt,
                                                 CCTK_REAL a, CCTK_REAL b,
                                                 CCTK_REAL c) {
  using std::sqrt;

  const auto r0 = pt.cake_inner_boundary_radius;
  const auto r1 = pt.cake_outer_boundary_radius;

  const auto factor_1 = 2 + Power(a, 2) * (1 + c) + Power(b, 2) * (1 + c);
  const auto factor_2 = 2 * factor_1;
  const auto denominator = sqrt(factor_1) * Power(factor_2, 3);
  const auto numerator = Sqrt(2) * (1 + c) * ((-1 + c) * r0 - (1 + c) * r1) *
                         (6 * Power(b, 2) * (1 + c) - Power(factor_2, 2));

  return -numerator / denominator;
}

/**
 * Second derivative of the local -> global core functions of the cake
 * coordinate transformation with respect to b and c.
 *
 * @param pt The cake patch data.
 * @param a The a local coordinate, ranging from (-1, 1)
 * @param b The b local coordiante, ranging from (-1, 1)
 * @param c the c local coordiante, ranging from (-1, 1)
 * @return The value of the core derivative with respect to c.
 */
inline CCTK_REAL local_to_global_cake_core_db_dc(const PatchTransformations &pt,
                                                 CCTK_REAL a, CCTK_REAL b,
                                                 CCTK_REAL c) {
  using std::sqrt;

  const auto r0 = pt.cake_inner_boundary_radius;
  const auto r1 = pt.cake_outer_boundary_radius;

  const auto factor_1 = 2 + Power(a, 2) * (1 + c) + Power(b, 2) * (1 + c);
  const auto factor_2 = 2 * factor_1;
  const auto denominator = sqrt(factor_1) * Power(factor_2, 3);
  const auto numerator = Sqrt(2) * b *
                         (3 * (Power(a, 2) + Power(b, 2)) * (1 + c) *
                              ((-1 + c) * r0 - (1 + c) * r1) +
                          2 * Power(factor_2, 2) * (r1 + c * (-r0 + r1)));

  return -numerator / denominator;
}

/**
 * Second derivative of the local -> global core function of the cake coordinate
 * transformation with respect to c and c.
 *
 * @param pt The cake patch data.
 * @param a The a local coordinate, ranging from (-1, 1)
 * @param b The b local coordiante, ranging from (-1, 1)
 * @param c the c local coordiante, ranging from (-1, 1)
 * @return The value of the core derivative with respect to c.
 */
inline CCTK_REAL local_to_global_cake_core_dc_dc(const PatchTransformations &pt,
                                                 CCTK_REAL a, CCTK_REAL b,
                                                 CCTK_REAL c) {
  const auto r0 = pt.cake_inner_boundary_radius;
  const auto r1 = pt.cake_outer_boundary_radius;

  const auto factor = 4 + 2 * Power(a, 2) * (1 + c) + 2 * Power(b, 2) * (1 + c);
  const auto denominator = sqrt(factor) * Power(factor, 3);
  const auto numerator =
      (Power(a, 2) + Power(b, 2)) *
      (-3 * (Power(a, 2) + Power(b, 2)) * ((-1 + c) * r0 - (1 + c) * r1) +
       2 * (r0 - r1) * Power(factor, 2));

  return numerator / denominator;
}

/**
 * The local to global coordinate transformation implementation.
 * This is where the actual coordinate transformations are performed but it is
 * also suitable for passing to CarpetX.
 *
 * @param pt The patch data
 * @param patch The index of the patch to transform.
 * @param local_vars The values of the local variables (a,b,c)
 * @return A spatial vector containing the coordinate transformations.
 */
CCTK_DEVICE CCTK_HOST svec_u local2global(const PatchTransformations &pt,
                                          int patch, const svec_u &local_vars) {

  svec_u global_vars = {0.0, 0.0, 0.0};

  const auto a = local_vars(0);
  const auto b = local_vars(1);
  const auto c = local_vars(2);

  const auto base_x = local_to_global_cake_core(pt, a, b, c);
  const auto base_y = b * base_x;
  const auto base_z = a * base_x;

  switch (patch) {
  case static_cast<int>(patch_piece::cartesian):
    global_vars = local_vars;
    break;
  case static_cast<int>(patch_piece::plus_x):
    global_vars = {base_x, base_y, base_z};
    break;
  case static_cast<int>(patch_piece::minus_x):
    global_vars = {-base_x, -base_y, base_z};
    break;
  case static_cast<int>(patch_piece::plus_y):
    global_vars = {-base_y, base_x, base_z};
    break;
  case static_cast<int>(patch_piece::minus_y):
    global_vars = {base_y, -base_x, base_z};
    break;
  case static_cast<int>(patch_piece::plus_z):
    global_vars = {-base_z, base_y, base_x};
    break;
  case static_cast<int>(patch_piece::minus_z):
    global_vars = {base_z, base_y, -base_x};
    break;
  default:
    CCTK_VERROR("No local -> global transformations available for patch %s",
                piece_name(static_cast<patch_piece>(patch)).c_str());
    break;
  }

  return global_vars;
}

/**
 * The global to local coordinate transformation implementation.
 * This is where the actual coordinate transformations are performed but it is
 * also suitable to be passed to CarpetX without any wrappers.
 *
 * @param pt The patch data
 * @param global_vars The values of the local global (x, y, z)
 * @return A tuple containing the patch piece and the global coordinate triplet.
 */
std::tuple<int, svec_u> global2local(const PatchTransformations &pt,
                                     const svec_u &global_vars) {
  const auto x = global_vars(0);
  const auto y = global_vars(1);
  const auto z = global_vars(2);

  const auto r0 = pt.cake_inner_boundary_radius;
  const auto r1 = pt.cake_outer_boundary_radius;

  const auto piece = get_owner_patch(pt, global_vars);
  svec_u local_vars = {0, 0, 0};

  const auto x2 = x * x;
  const auto y2 = y * y;
  const auto z2 = z * z;
  const auto r02 = r0 * r0;
  const auto r12 = r1 * r1;

  switch (static_cast<int>(piece)) {

  case static_cast<int>(patch_piece::cartesian):
    local_vars = global_vars;
    break;

  case static_cast<int>(patch_piece::plus_x):
    local_vars = {
        z / x, y / x,
        (r02 - r12 + y2 + z2 +
         sqrt(4 * r12 * x2 + Power(y2 + z2, 2) + 4 * r02 * (x2 + y2 + z2) -
              4 * r0 * r1 * (2 * x2 + y2 + z2))) /
            Power(r0 - r1, 2)};
    break;

  case static_cast<int>(patch_piece::minus_x):
    local_vars = {
        -z / x, y / x,
        (r02 - r12 + y2 + z2 +
         sqrt(4 * r12 * x2 + Power(y2 + z2, 2) + 4 * r02 * (x2 + y2 + z2) -
              4 * r0 * r1 * (2 * x2 + y2 + z2))) /
            Power(r0 - r1, 2)};
    break;

  case static_cast<int>(patch_piece::plus_y):
    local_vars = {
        z / y, -x / y,
        (r02 - r12 + x2 + z2 +
         sqrt(4 * r12 * y2 + Power(x2 + z2, 2) + 4 * r02 * (x2 + y2 + z2) -
              4 * r0 * r1 * (x2 + 2 * y2 + z2))) /
            Power(r0 - r1, 2)};
    break;

  case static_cast<int>(patch_piece::minus_y):
    local_vars = {
        -z / y, -x / y,
        (r02 - r12 + x2 + z2 +
         sqrt(4 * r12 * y2 + Power(x2 + z2, 2) + 4 * r02 * (x2 + y2 + z2) -
              4 * r0 * r1 * (x2 + 2 * y2 + z2))) /
            Power(r0 - r1, 2)};
    break;

  case static_cast<int>(patch_piece::plus_z):
    local_vars = {-x / z, y / z,
                  (r02 - r12 + x2 + y2 +
                   sqrt((x2 + y2) * (4 * r0 * (r0 - r1) + x2 + y2) +
                        4 * Power(r0 - r1, 2) * z2)) /
                      Power(r0 - r1, 2)};
    break;

  case static_cast<int>(patch_piece::minus_z):
    local_vars = {-x / z, -y / z,
                  (r02 - r12 + x2 + y2 +
                   sqrt((x2 + y2) * (4 * r0 * (r0 - r1) + x2 + y2) +
                        4 * Power(r0 - r1, 2) * z2)) /
                      Power(r0 - r1, 2)};
    break;

  default:
    CCTK_VERROR("No global -> local transformations available for patch %s",
                piece_name(piece).c_str());
    break;
  }

  return std::make_tuple(piece_idx(piece), local_vars);
}

/*******************************************************************
 * TODO: Review Jacobians                                          *
 *******************************************************************/

/**
 * This function computes the local to global coordinate transformation
 * jacobian and it's derivative. This is better than computng the two
 * separatelly because the derivative of J depends on J itself. J[i][j] =
 * $J_{i j} = \frac{d x^i}{d x^j}$. dJ[i][j][k] = $dJ_{i j k} = \frac{d^2
 * x^i}{d x^j d x^k}$.
 *
 * @tparam p The patch piece to perform the coordinate transformations for.
 * @param pt The patch data
 * @param local_vars The values of the local variables (a,b,c)
 * @return A tuple, containing the the local to global transformation the
 * jacobian and it's derivative.
 */
template <patch_piece p>
inline std::tuple<svec_u, jac_t, djac_t>
jacobians(const PatchTransformations &pt, const svec_u &local_vars) {
  jac_t J = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  djac_t dJ = {zero<smat_d>()(), zero<smat_d>()(), zero<smat_d>()()};

  // Since the cartesian patch has a trivial jacobian derivative, we can
  // skip the whole calculation
  if constexpr (p == patch_piece::exterior) {
    CCTK_WARN(CCTK_WARN_DEBUG,
              "Attempted to compute the jacobians of the exterior of the grid. "
              "Returning unit jacobian and zero jacobian derivative.");

  } else if constexpr (p != patch_piece::cartesian) {

    const auto a = local_vars(0);
    const auto b = local_vars(1);
    const auto c = local_vars(2);

    const auto core = local_to_global_cake_core(pt, a, b, c);
    const auto d_core_da = local_to_global_cake_core_da(pt, a, b, c);
    const auto d_core_db = local_to_global_cake_core_db(pt, a, b, c);
    const auto d_core_dc = local_to_global_cake_core_dc(pt, a, b, c);
    const auto d2_core_da_da = local_to_global_cake_core_da_da(pt, a, b, c);
    const auto d2_core_da_db = local_to_global_cake_core_da_db(pt, a, b, c);
    const auto d2_core_da_dc = local_to_global_cake_core_da_dc(pt, a, b, c);
    const auto d2_core_db_db = local_to_global_cake_core_db_db(pt, a, b, c);
    const auto d2_core_db_dc = local_to_global_cake_core_db_dc(pt, a, b, c);
    const auto d2_core_dc_dc = local_to_global_cake_core_dc_dc(pt, a, b, c);

    if constexpr (p == patch_piece::minus_x) {
      CAKE_MINUS_X_JACOBIAN
      CAKE_MINUS_X_DJACOBIAN
    } else if constexpr (p == patch_piece::plus_x) {
      CAKE_PLUS_X_JACOBIAN
      CAKE_PLUS_X_DJACOBIAN
    } else if constexpr (p == patch_piece::minus_y) {
      CAKE_MINUS_Y_JACOBIAN
      CAKE_MINUS_Y_DJACOBIAN
    } else if constexpr (p == patch_piece::plus_y) {
      CAKE_PLUS_Y_JACOBIAN
      CAKE_PLUS_Y_DJACOBIAN
    } else if constexpr (p == patch_piece::minus_z) {
      CAKE_MINUS_Z_JACOBIAN
      CAKE_MINUS_Z_DJACOBIAN
    } else if constexpr (p == patch_piece::plus_z) {
      CAKE_PLUS_Z_JACOBIAN
      CAKE_PLUS_Z_DJACOBIAN
    }
  }

  return std::make_tuple(local2global(pt, static_cast<int>(p), local_vars), J,
                         dJ);
}

/**
 * This function retrieves the jacobian from the jacobians function
 *
 * @tparam p The patch piece whose jacobian is to be retrieved.
 * @param pt The patch data
 * @param local_vars The values of the local variables (a,b,c)
 * @return A tuple, containing the  local to global transformation and the
 * jacobian.
 */
template <patch_piece p>
inline std::tuple<svec_u, jac_t> jacobian(const PatchTransformations &pt,
                                          const svec_u &local_vars) {
  const auto x = jacobians<p>(pt, local_vars);
  return std::make_tuple(std::get<0>(x), std::get<1>(x));
}

/*******************************************************************
 * CarpetX wrappers and utilities                                  *
 *******************************************************************/

/**
 * Wrapper around jacobian to make is suitable for passing to CarpetX.
 *
 * @param pt The patch data
 * @param patch The index of the patch to transform.
 * @param local_vars The values of the local variables (a,b,c)
 * @return A tuple containing the local to global coordinate transformation
 * and the local to global jacobian matrix.
 */
std::tuple<svec_u, jac_t> dlocal_dglobal(const PatchTransformations &pt,
                                         int patch, const svec_u &local_vars) {
  switch (patch) {
  case piece_idx(patch_piece::cartesian):
    return jacobian<patch_piece::cartesian>(pt, local_vars);
  case piece_idx(patch_piece::minus_x):
    return jacobian<patch_piece::minus_x>(pt, local_vars);
  case piece_idx(patch_piece::plus_x):
    return jacobian<patch_piece::plus_x>(pt, local_vars);
  case piece_idx(patch_piece::minus_y):
    return jacobian<patch_piece::minus_y>(pt, local_vars);
  case piece_idx(patch_piece::plus_y):
    return jacobian<patch_piece::plus_y>(pt, local_vars);
  case piece_idx(patch_piece::minus_z):
    return jacobian<patch_piece::minus_z>(pt, local_vars);
  case piece_idx(patch_piece::plus_z):
    return jacobian<patch_piece::plus_z>(pt, local_vars);
  default:
    CCTK_VERROR("No dlocal_dglobal available for patch index %d", patch);
    break;
  }
}

/**
 * Wrapper around jacobian to make is suitable for passing to CarpetX.
 *
 * @param pt The patch data
 * @param patch The index of the patch to transform.
 * @param local_vars The values of the local variables (a,b,c)
 * @return A tuple containing the local to global coordinate transformation,
 * the local to global jacobian matrix and it's derivative.
 */
std::tuple<svec_u, jac_t, djac_t>
d2local_dglobal2(const PatchTransformations &pt, int patch,
                 const svec_u &local_vars) {
  switch (patch) {
  case piece_idx(patch_piece::cartesian):
    return jacobians<patch_piece::cartesian>(pt, local_vars);
  case piece_idx(patch_piece::minus_x):
    return jacobians<patch_piece::minus_x>(pt, local_vars);
  case piece_idx(patch_piece::plus_x):
    return jacobians<patch_piece::plus_x>(pt, local_vars);
  case piece_idx(patch_piece::minus_y):
    return jacobians<patch_piece::minus_y>(pt, local_vars);
  case piece_idx(patch_piece::plus_y):
    return jacobians<patch_piece::plus_y>(pt, local_vars);
  case piece_idx(patch_piece::minus_z):
    return jacobians<patch_piece::minus_z>(pt, local_vars);
  case piece_idx(patch_piece::plus_z):
    return jacobians<patch_piece::plus_z>(pt, local_vars);
  default:
    CCTK_VERROR("No dlocal_dglobal available for patch index %d", patch);
    break;
  }
}

/**
 * Creates a cake patch
 *
 * @tparam p The piece of the patc to make.
 * @param pt The patch transformation object with patch data.
 * @return The constructed patch piece.
 */
template <patch_piece p> Patch make_patch(const PatchTransformations &pt) {
  const auto cartesian_ncells_i = pt.cake_cartesian_ncells_i;
  const auto cartesian_ncells_j = pt.cake_cartesian_ncells_j;
  const auto cartesian_ncells_k = pt.cake_cartesian_ncells_k;

  const auto spherical_radial_cells = pt.cake_radial_cells;
  const auto spherical_angular_cells = pt.cake_angular_cells;

  // This is the most likelly configuration to occur. The only exception is
  // the cartesian patch where the contents of these fields are overwritten
  Patch patch;
  patch.ncells = {spherical_radial_cells, spherical_angular_cells,
                  spherical_angular_cells};

  patch.xmin = {-1, -1, -1};
  patch.xmax = {+1, +1, +1};

  patch.is_cartesian = false;

  // This is different for every patch but returning it only once eliminates
  // repetition
  PatchFace m_i = {true, -1};
  PatchFace p_i = {true, -1};
  PatchFace m_j = {true, -1};
  PatchFace p_j = {true, -1};
  PatchFace m_k = {true, -1};
  PatchFace p_k = {true, -1};

  if constexpr (p == patch_piece::cartesian) {
    const auto f = pt.cake_inner_boundary_radius / 2;

    patch.ncells = {cartesian_ncells_i, cartesian_ncells_j, cartesian_ncells_k};
    patch.xmin = {-f, -f, -f};
    patch.xmax = {f, f, f};
    patch.is_cartesian = true;

    m_i = PatchFace{false, piece_idx(patch_piece::minus_x)};
    p_i = PatchFace{false, piece_idx(patch_piece::plus_x)};
    m_j = PatchFace{false, piece_idx(patch_piece::minus_y)};
    p_j = PatchFace{false, piece_idx(patch_piece::plus_y)};
    m_k = PatchFace{false, piece_idx(patch_piece::minus_z)};
    p_k = PatchFace{false, piece_idx(patch_piece::plus_z)};

  } else if constexpr (p == patch_piece::plus_x) {
    m_i = PatchFace{false, piece_idx(patch_piece::cartesian)};
    p_i = PatchFace{true, piece_idx(patch_piece::exterior)};
    m_j = PatchFace{false, piece_idx(patch_piece::minus_y)};
    p_j = PatchFace{false, piece_idx(patch_piece::plus_y)};
    m_k = PatchFace{false, piece_idx(patch_piece::minus_z)};
    p_k = PatchFace{false, piece_idx(patch_piece::plus_z)};

  } else if constexpr (p == patch_piece::minus_x) {
    m_i = PatchFace{true, piece_idx(patch_piece::exterior)};
    p_i = PatchFace{false, piece_idx(patch_piece::cartesian)};
    m_j = PatchFace{false, piece_idx(patch_piece::minus_y)};
    p_j = PatchFace{false, piece_idx(patch_piece::plus_y)};
    m_k = PatchFace{false, piece_idx(patch_piece::minus_z)};
    p_k = PatchFace{false, piece_idx(patch_piece::plus_z)};

  } else if constexpr (p == patch_piece::plus_y) {
    m_i = PatchFace{false, piece_idx(patch_piece::minus_x)};
    p_i = PatchFace{false, piece_idx(patch_piece::plus_x)};
    m_j = PatchFace{false, piece_idx(patch_piece::cartesian)};
    p_j = PatchFace{true, piece_idx(patch_piece::exterior)};
    m_k = PatchFace{false, piece_idx(patch_piece::minus_z)};
    p_k = PatchFace{false, piece_idx(patch_piece::plus_z)};

  } else if constexpr (p == patch_piece::minus_y) {
    m_i = PatchFace{false, piece_idx(patch_piece::minus_x)};
    p_i = PatchFace{false, piece_idx(patch_piece::plus_x)};
    m_j = PatchFace{true, piece_idx(patch_piece::exterior)};
    p_j = PatchFace{false, piece_idx(patch_piece::cartesian)};
    m_k = PatchFace{false, piece_idx(patch_piece::minus_z)};
    p_k = PatchFace{false, piece_idx(patch_piece::plus_z)};

  } else if constexpr (p == patch_piece::minus_z) {
    m_i = PatchFace{false, piece_idx(patch_piece::minus_x)};
    p_i = PatchFace{false, piece_idx(patch_piece::plus_x)};
    m_j = PatchFace{false, piece_idx(patch_piece::minus_x)};
    p_j = PatchFace{false, piece_idx(patch_piece::plus_y)};
    m_k = PatchFace{true, piece_idx(patch_piece::exterior)};
    p_k = PatchFace{false, piece_idx(patch_piece::cartesian)};

  } else if constexpr (p == patch_piece::plus_z) {
    m_i = PatchFace{false, piece_idx(patch_piece::minus_x)};
    p_i = PatchFace{false, piece_idx(patch_piece::plus_x)};
    m_j = PatchFace{false, piece_idx(patch_piece::minus_x)};
    p_j = PatchFace{false, piece_idx(patch_piece::plus_y)};
    m_k = PatchFace{false, piece_idx(patch_piece::cartesian)};
    p_k = PatchFace{true, piece_idx(patch_piece::exterior)};
  }

  patch.faces = {{m_i, m_j, m_k}, {p_i, p_j, p_k}};
  return patch;
}

} // namespace Cake

/**
 * Creates a Cake patch system
 *
 * @return A PatchSystem object with Cake data and functions
 */
PatchSystem SetupCake() {
  PatchTransformations pt;
  pt.global2local = &Cake::global2local;
  pt.local2global = &Cake::local2global;
  pt.dlocal_dglobal = &Cake::dlocal_dglobal;
  pt.d2local_dglobal2 = &Cake::d2local_dglobal2;
  pt.global2local_device = nullptr;
  pt.local2global_device = nullptr;
  pt.dlocal_dglobal_device = nullptr;
  pt.d2local_dglobal2_device = nullptr;

  const auto patches =
      std::vector<Patch>{Cake::make_patch<Cake::patch_piece::cartesian>(pt),
                         Cake::make_patch<Cake::patch_piece::minus_x>(pt),
                         Cake::make_patch<Cake::patch_piece::plus_x>(pt),
                         Cake::make_patch<Cake::patch_piece::minus_y>(pt),
                         Cake::make_patch<Cake::patch_piece::plus_y>(pt),
                         Cake::make_patch<Cake::patch_piece::minus_z>(pt),
                         Cake::make_patch<Cake::patch_piece::plus_z>(pt)};

  return PatchSystem(patches, std::move(pt));
}

} // namespace MultiPatch