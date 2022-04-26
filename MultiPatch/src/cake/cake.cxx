#include "multipatch.hxx"
#include "tests.hxx"
#include "cake_cartesian_jac.hxx"
#include "cake_plus_x_jac.hxx"
#include "cake_plus_y_jac.hxx"
#include "cake_plus_z_jac.hxx"
#include "cake_minus_x_jac.hxx"
#include "cake_minus_y_jac.hxx"
#include "cake_minus_z_jac.hxx"
#include "cake.hxx"

#include <cassert>
#include <cmath>

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
CCTK_DEVICE CCTK_HOST inline CCTK_REAL
local_to_global_cake_core(const PatchTransformations &pt, CCTK_REAL a,
                          CCTK_REAL b, CCTK_REAL c) {
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
#ifndef __CUDACC__
    CCTK_VERROR("No local -> global transformations available for patch %s",
                piece_name(static_cast<patch_piece>(patch)).c_str());
#else
    assert(0);
#endif
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
CCTK_DEVICE CCTK_HOST std_tuple<int, svec_u>
global2local(const PatchTransformations &pt, const svec_u &global_vars) {
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
#ifndef __CUDACC__
    CCTK_VERROR("No global -> local transformations available for patch %s",
                piece_name(piece).c_str());
#else
    assert(0);
#endif
    break;
  }

  return std_make_tuple(static_cast<int>(piece), local_vars);
}

/**
 * This function computes the local to global coordinate
 * transformation, the jacobian and it's derivative. It can be passed directly
 * to CarpetX.
 *
 * Note that:
 * J(i)(j) = $J^{i}_{j} = \frac{d a^i}{d x^j}$.
 * dJ(i)(j,k) = $dJ^{i}_{j k} = \frac{d}{d x^k} \left( \frac{d a^i}{d x^j}
 * \right)$.
 *
 * TODO: Erik says: "You have six files that contain very similar code.
 * presumably, they differ in their permutations of x, y, and z, and a few minus
 * signs. you could create just one file, and call the functions in that file
 * six times with different arguments. this might reduce the amount of code
 * generated, and might speed up things. or it might not. don’t worry about this
 * now. this is just a thought for the future."
 *
 * @param pt The patch data
 * @param patch The index of the patch to transform.
 * @param local_vars The values of the local variables (a,b,c)
 * @return A tuple containing the local to global coordinate transformation,
 * the local to global jacobian matrix and it's derivative.
 */
CCTK_DEVICE CCTK_HOST std_tuple<svec_u, jac_t, djac_t>
d2local_dglobal2(const PatchTransformations &pt, int patch,
                 const svec_u &local_vars) {

  auto local_to_global_result = pt.local2global(pt, patch, local_vars);

  auto jacobian_results = cake_cartesian_jac(pt, local_vars);

  switch (patch) {

  case static_cast<int>(patch_piece::cartesian):
    // No action is necessary
    break;

  case static_cast<int>(patch_piece::plus_x):
    jacobian_results = cake_plus_x_jac(pt, local_vars);
    break;

  case static_cast<int>(patch_piece::minus_x):
    jacobian_results = cake_minus_x_jac(pt, local_vars);
    break;

  case static_cast<int>(patch_piece::plus_y):
    jacobian_results = cake_plus_y_jac(pt, local_vars);
    break;

  case static_cast<int>(patch_piece::minus_y):
    jacobian_results = cake_minus_y_jac(pt, local_vars);
    break;

  case static_cast<int>(patch_piece::plus_z):
    jacobian_results = cake_plus_z_jac(pt, local_vars);
    break;

  case static_cast<int>(patch_piece::minus_z):
    jacobian_results = cake_minus_z_jac(pt, local_vars);
    break;

  default:
#ifndef __CUDACC__
    CCTK_VERROR("No jacobians available for patch %s",
                piece_name(static_cast<patch_piece>(patch)).c_str());
#else
    assert(0);
#endif
    break;
  }

  return std_make_tuple(local_to_global_result, std::get<0>(jacobian_results),
                         std::get<1>(jacobian_results));
} // namespace Cake

/**
 * This function computes the local to global coordinate transformation and the
 * jacobian. It can be passed directly to CarpetX.
 *
 * Note that:
 * J(i)(j) = $J^{i}_{j} = \frac{d a^i}{d x^j}$.
 *
 * @param pt The patch data
 * @param patch The index of the patch to transform.
 * @param local_vars The values of the local variables (a,b,c)
 * @return A tuple containing the local to global coordinate transformation
 * and the local to global jacobian matrix.
 */
CCTK_DEVICE CCTK_HOST std_tuple<svec_u, jac_t>
dlocal_dglobal(const PatchTransformations &pt, int patch,
               const svec_u &local_vars) {

  const auto data = pt.d2local_dglobal2(pt, patch, local_vars);
  return std_make_tuple(std::get<0>(data), std::get<1>(data));
}

/**
 * TODO: Review. This is probably causing a segfault.
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

  PatchFace m_i = {true, -1};
  PatchFace p_i = {true, -1};
  PatchFace m_j = {true, -1};
  PatchFace p_j = {true, -1};
  PatchFace m_k = {true, -1};
  PatchFace p_k = {true, -1};

  patch.faces = {{m_i, m_j, m_k}, {p_i, p_j, p_k}};
  return patch;
}

} // namespace Cake

/**
 * TODO: Review. Are the correct device functions and signatures being used?
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
  pt.global2local_device = &Cake::global2local;
  pt.local2global_device = &Cake::local2global;
  pt.dlocal_dglobal_device = &Cake::dlocal_dglobal;
  pt.d2local_dglobal2_device = &Cake::d2local_dglobal2;

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