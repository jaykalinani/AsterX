#include "cake.hxx"

MultiPatch::Cake::patch_piece
MultiPatch::Cake::get_owner_patch(const PatchTransformations &pt,
                                  const svec_u &global_vars) {
  using MultiPatchTests::isapprox;
  using std::abs;

  const auto x = global_vars(0);
  const auto y = global_vars(1);
  const auto z = global_vars(2);

  const auto f = pt.cake_inner_boundary_radius / 2;
  const auto Rf = pt.cake_outer_boundary_radius;
  const auto r = std::sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2));

  if (within(x, f) && within(y, f) && within(z, f)) {
    return patch_piece::cartesian;
  } else if (at_boundary(x, f) || at_boundary(y, f) || at_boundary(z, f)) {
    return patch_piece::inner_boundary;
  } else if (r < Rf || !isapprox(r, Rf)) {
    /* We are in one of the 6 patches or theur inner boundaries. To determine
     * the direction of the patch, we must determine which coordinate (here
     * represented by an index of the global_vars vector) has the largest
     * absolute value:
     */
    const auto abs_x = abs(x);
    const auto abs_y = abs(y);
    const auto abs_z = abs(z);

    const auto max_pair_idx = abs_x > abs_y ? 0 : 1;
    const auto max_pair = abs(global_vars(max_pair_idx));

    const auto max_coord_idx = abs_z > max_pair ? 2 : max_pair_idx;

    /* After the direction is know, we must dtermine the sense, which is
     * represente by the sign of the coordinate with the largest absolute value.
     */
    const int max_coord_sign = global_vars(max_coord_idx) > 0.0 ? 1 : -1;

    if (max_coord_sign < 0 && max_coord_idx == 0) {
      return patch_piece::minus_x;
    } else if (max_coord_sign > 0 && max_coord_idx == 0) {
      return patch_piece::plus_x;
    } else if (max_coord_sign < 0 && max_coord_idx == 1) {
      return patch_piece::minus_y;
    } else if (max_coord_sign > 0 && max_coord_idx == 1) {
      return patch_piece::plus_y;
    } else if (max_coord_sign < 0 && max_coord_idx == 2) {
      return patch_piece::minus_z;
    } else if (max_coord_sign > 0 && max_coord_idx == 2) {
      return patch_piece::plus_z;
    }
  } else if (isapprox(r, Rf)) {
    return patch_piece::outer_boundary;
  } else {
    CCTK_VERROR("Coordinate triplet (%f, %f, %f) outside of simulation domain",
                x, y, z);
  }
  return patch_piece::exterior;
}