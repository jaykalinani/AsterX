#include "cake.hxx"

#include <cassert>

CCTK_DEVICE CCTK_HOST MultiPatch::Cake::patch_piece
MultiPatch::Cake::get_owner_patch(const PatchTransformations &pt,
                                  const svec_u &global_vars) {
  using MultiPatchTests::at_boundary;
  using MultiPatchTests::isapprox;
  using MultiPatchTests::within;
  using std::abs;
  using std::sqrt;

  const auto x = global_vars(0);
  const auto y = global_vars(1);
  const auto z = global_vars(2);

  const auto r0 = pt.cake_inner_boundary_radius;
  const auto r1 = pt.cake_outer_boundary_radius;
  const auto r = sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2));

  patch_piece piece = patch_piece::exterior;

  // Are we on the interior?
  if (r < r1 || !isapprox(r, r1)) {

    // Are we on the Cartesian core?
    if (within(x, r0) && within(y, r0) && within(z, r0)) {
      piece = patch_piece::cartesian;
    } else if (at_boundary(x, r0) || at_boundary(y, r0) || at_boundary(z, r0)) {
      piece = patch_piece::inner_boundary;
    }

    // Are we on the 6 patches?
    else if (!isapprox(r, r1)) {

      /* The patch direction is the determined by the coordinate with the
       * largest absolute value
       */
      const auto abs_x = abs(x);
      const auto abs_y = abs(y);
      const auto abs_z = abs(z);

      const auto max_pair_idx = abs_x > abs_y ? 0 : 1;
      const auto max_pair = abs(global_vars(max_pair_idx));

      const auto max_coord_idx = abs_z > max_pair ? 2 : max_pair_idx;

      /* After the direction is know, we must dtermine the sense, which is
       * represented by the sign of the coordinate with the largest absolute
       * value.
       */
      const int max_coord_sign = global_vars(max_coord_idx) > 0.0 ? 1 : -1;

      if (max_coord_sign < 0 && max_coord_idx == 0) {
        piece = patch_piece::minus_x;
      } else if (max_coord_sign > 0 && max_coord_idx == 0) {
        piece = patch_piece::plus_x;
      } else if (max_coord_sign < 0 && max_coord_idx == 1) {
        piece = patch_piece::minus_y;
      } else if (max_coord_sign > 0 && max_coord_idx == 1) {
        piece = patch_piece::plus_y;
      } else if (max_coord_sign < 0 && max_coord_idx == 2) {
        piece = patch_piece::minus_z;
      } else if (max_coord_sign > 0 && max_coord_idx == 2) {
        piece = patch_piece::plus_z;
      }
    }

    // Are we on the exterior boundary?
    else if (isapprox(r, r1)) {
      piece = patch_piece::outer_boundary;
    }

    // We dont know where we are
    else {
#ifndef __CUDACC__
      CCTK_VERROR("Coordinate triplet (%f, %f, %f) cannot be located within "
                  "the simulation domain",
                  x, y, z);
#else
      assert(0);
#endif
    }
  }

  return piece;
}
