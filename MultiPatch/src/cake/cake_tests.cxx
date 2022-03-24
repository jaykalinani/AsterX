#include "cake.hxx"
#include "tests.hxx"

#include <cctk_Parameters.h>

#include <string>
#include <utility>
#include <cmath>

namespace MultiPatch {
namespace CakeTests {

/**
 * TODO: doc
 */
std::string patch_owner_test(const PatchTransformations &pt,
                             const MultiPatch::Cake::svec_u &x,
                             MultiPatch::Cake::patch_piece expected) {

  using namespace MultiPatch::Cake;
  using namespace MultiPatchTests;

  std::string msg{"has "};

  const auto owner_patch = get_owner_patch(pt, x);

  if (owner_patch == expected) {
    msg += colored<string_color::green>("PASSED");
  } else {
    msg += colored<string_color::red>("FAILED");
    msg += ". Reason: Expected to get patch ";
    msg += piece_name(expected);
    msg += " and got ";
    msg += piece_name(owner_patch);
  }

  msg += ".";
  return msg;
}

/**
 * Tests if local2global(global2local(global)) == global
 * TODO: doc
 */
std::string global_identity_test(const PatchTransformations &pt,
                                 const MultiPatch::Cake::svec_u &global_vars) {

  using MultiPatchTests::colored;
  using MultiPatchTests::isapprox;
  using MultiPatchTests::string_color;

  std::string msg{"has "};

  const auto g2l = pt.global2local(pt, global_vars);

  const auto owner_patch_idx = std::get<0>(g2l);
  const auto local_vars = std::get<1>(g2l);

  const auto l2g = pt.local2global(pt, owner_patch_idx, local_vars);

  const auto test1 = isapprox(l2g(0), global_vars(0));
  const auto test2 = isapprox(l2g(1), global_vars(1));
  const auto test3 = isapprox(l2g(2), global_vars(2));

  if (test1 && test2 && test3) {
    msg += colored<string_color::green>("PASSED");
  } else {
    msg += colored<string_color::red>("FAILED");
    msg += ". Reason: ";

    if (!test1) {
      msg += std::to_string(l2g(0));
      msg += " =/= ";
      msg += std::to_string(global_vars(0));
      msg += ". ";
    }

    if (!test2) {
      msg += std::to_string(l2g(1));
      msg += " =/= ";
      msg += std::to_string(global_vars(1));
      msg += ". ";
    }

    if (!test3) {
      msg += std::to_string(l2g(2));
      msg += " =/= ";
      msg += std::to_string(global_vars(2));
      msg += ". ";
    }
  }

  return msg;
}

/**
 * Tests if global2local(local2global(local, patch)) == (local, patch)
 * TODO: doc
 */
std::string local_identity_test(const PatchTransformations &pt, int patch,
                                const MultiPatch::Cake::svec_u &local_point) {

  using MultiPatch::Cake::patch_piece;
  using MultiPatch::Cake::piece_name;
  using MultiPatchTests::colored;
  using MultiPatchTests::isapprox;
  using MultiPatchTests::string_color;

  std::string msg{"has "};

  const auto l2g = pt.local2global(pt, patch, local_point);
  const auto g2l = pt.global2local(pt, l2g);

  const auto computed_patch_idx = std::get<0>(g2l);
  const auto computed_local_point = std::get<1>(g2l);

  const bool test1 = computed_patch_idx == patch;
  const bool test2 = isapprox(computed_local_point(0), local_point(0));
  const bool test3 = isapprox(computed_local_point(1), local_point(1));
  const bool test4 = isapprox(computed_local_point(2), local_point(2));

  if (test1 && test2 && test3 && test4) {
    msg += colored<string_color::green>("PASSED");
  } else {
    msg += colored<string_color::red>("FAILED");
    msg += ". Reason: ";

    if (!test1) {
      msg += "the computed patch is ";
      msg += piece_name(static_cast<patch_piece>(computed_patch_idx));
      msg += ". ";
    }

    if (!test2) {
      msg += "the computed local coordinate a is ";
      msg += std::to_string(computed_local_point(0));
      msg += ". ";
    }

    if (!test3) {
      msg += "the computed local coordinate b is ";
      msg += std::to_string(computed_local_point(1));
      msg += ". ";
    }

    if (!test4) {
      msg += "the computed local coordinate c is ";
      msg += std::to_string(computed_local_point(2));
      msg += ". ";
    }
  }

  return msg;
}

} // namespace CakeTests
} // namespace MultiPatch

extern "C" void run_cake_tests() {
  DECLARE_CCTK_PARAMETERS

  using std::mt19937_64;
  using std::uniform_int_distribution;
  using std::uniform_real_distribution;

  using std::cos;
  using std::sin;

  using MultiPatch::Cake::patch_piece;
  using MultiPatch::Cake::piece_name;
  using MultiPatch::Cake::svec_u;

  using MultiPatch::dim;
  using MultiPatch::SetupCake;
  using MultiPatchTests::random_seed;

  using namespace MultiPatch::CakeTests;

  const auto ps = SetupCake();
  const auto &pt = ps.transformations;

  const auto r0 = pt.cake_inner_boundary_radius;
  const auto r1 = pt.cake_outer_boundary_radius;

  const auto cube_midpoint = r0 / 2.0;
  const auto patch_midpoint = r0 + (r1 - r0) / 2.0;

  CCTK_INFO("Running cake patch tests:");

  /*
   * Fixed point tests
   */
  const std::array<std::pair<svec_u, patch_piece>, 27> owner_test_data = {
      // Cartesian cube interior
      std::make_pair(svec_u{cube_midpoint, cube_midpoint, cube_midpoint},
                     patch_piece::cartesian),
      std::make_pair(svec_u{-cube_midpoint, cube_midpoint, cube_midpoint},
                     patch_piece::cartesian),
      std::make_pair(svec_u{cube_midpoint, -cube_midpoint, cube_midpoint},
                     patch_piece::cartesian),
      std::make_pair(svec_u{cube_midpoint, cube_midpoint, -cube_midpoint},
                     patch_piece::cartesian),

      // Cartesian cube corners
      std::make_pair(svec_u{r0, r0, r0}, patch_piece::inner_boundary),
      std::make_pair(svec_u{-r0, r0, r0}, patch_piece::inner_boundary),
      std::make_pair(svec_u{r0, -r0, r0}, patch_piece::inner_boundary),
      std::make_pair(svec_u{r0, r0, -r0}, patch_piece::inner_boundary),

      // Middle of the +x patch
      std::make_pair(svec_u{r0, 0.0, 0.0}, patch_piece::inner_boundary),
      std::make_pair(svec_u{patch_midpoint, 0.0, 0.0}, patch_piece::plus_x),
      std::make_pair(svec_u{r1, 0.0, 0.0}, patch_piece::outer_boundary),

      // Middle of the -x patch
      std::make_pair(svec_u{-r0, 0.0, 0.0}, patch_piece::inner_boundary),
      std::make_pair(svec_u{-patch_midpoint, 0.0, 0.0}, patch_piece::minus_x),
      std::make_pair(svec_u{-r1, 0.0, 0.0}, patch_piece::outer_boundary),

      // Middle of the +y patch
      std::make_pair(svec_u{0.0, r0, 0.0}, patch_piece::inner_boundary),
      std::make_pair(svec_u{0.0, patch_midpoint, 0.0}, patch_piece::plus_y),
      std::make_pair(svec_u{0.0, r1, 0.0}, patch_piece::outer_boundary),

      // Middle of the -y patch
      std::make_pair(svec_u{0.0, -r0, 0.0}, patch_piece::inner_boundary),
      std::make_pair(svec_u{0.0, -patch_midpoint, 0.0}, patch_piece::minus_y),
      std::make_pair(svec_u{0.0, -r1, 0.0}, patch_piece::outer_boundary),

      // Middle of the +z patch
      std::make_pair(svec_u{0.0, 0.0, r0}, patch_piece::inner_boundary),
      std::make_pair(svec_u{0.0, 0.0, patch_midpoint}, patch_piece::plus_z),
      std::make_pair(svec_u{0.0, 0.0, r1}, patch_piece::outer_boundary),

      // Middle of the -z patch
      std::make_pair(svec_u{0.0, 0.0, -r0}, patch_piece::inner_boundary),
      std::make_pair(svec_u{0.0, 0.0, -patch_midpoint}, patch_piece::minus_z),
      std::make_pair(svec_u{0.0, 0.0, -r1}, patch_piece::outer_boundary),

      std::make_pair(svec_u{r1 + 1.0, r1 + 1.0, r1 + 1.0},
                     patch_piece::exterior)};

  // Tests if the patch owner is correct.
  for (const auto &data : owner_test_data) {
    CCTK_VINFO("  Patch owner test at point (%f, %f, %f) %s", data.first(0),
               data.first(1), data.first(2),
               patch_owner_test(pt, data.first, data.second).c_str());
  }

  /*
   * Random point tests
   */
  mt19937_64 engine(random_seed);

  constexpr const double pi = 3.14159265358979323846;

  uniform_real_distribution<CCTK_REAL> r_distrib(0.0, r1);
  uniform_real_distribution<CCTK_REAL> theta_distrib(0.0, pi);
  uniform_real_distribution<CCTK_REAL> phi_distrib(0.0, 2.0 * pi);

  uniform_real_distribution<CCTK_REAL> local_distrib(-1.0, 1.0);

  uniform_int_distribution<int> patch_distrib(0, 6);

  auto global_point = svec_u{0.0, 0.0, 0.0};
  auto local_point = svec_u{0.0, 0.0, 0.0};
  double r = 0.0, theta = 0.0, phi = 0.0;
  int patch = 0;

  // Tests if local2global(global2local(global)) == global
  for (int i = 0; i < repeat_tests; i++) {
    r = r_distrib(engine);
    theta = theta_distrib(engine);
    phi = phi_distrib(engine);

    global_point = {r * sin(theta) * cos(phi), r * sin(theta) * sin(phi),
                    r * cos(theta)};

    CCTK_VINFO("  local2global(global2local(global)) transformation test at "
               "point (%f, %f, %f) %s",
               global_point(0), global_point(1), global_point(2),
               global_identity_test(pt, global_point).c_str());
  }

  // Tests if global2local(local2global(local, patch)) == (local, patch)
  for (int i = 0; i < repeat_tests; i++) {
    local_point = {local_distrib(engine), local_distrib(engine),
                   local_distrib(engine)};

    patch = patch_distrib(engine);

    CCTK_VINFO("  global2local(local2global(local, patch)) transformation test "
               "at point (%f, %f, %f) patch %s %s",
               local_point(0), local_point(1), local_point(2),
               piece_name(static_cast<patch_piece>(patch)).c_str(),
               local_identity_test(pt, patch, local_point).c_str());
  }
}