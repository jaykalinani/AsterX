#include "cake.hxx"
#include "tests.hxx"

#include <cctk_Parameters.h>

#include <utility>

namespace MultiPatch {
namespace CakeTests {

std::string patch_owner(const PatchTransformations &pt,
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

std::string global2local_test(const PatchTransformations &pt,
                              const MultiPatch::Cake::svec_u &x) {
  std::string msg{"has returnd "};

  const auto g2l = pt.global2local(pt, x);

  msg += "(";
  msg += std::to_string(std::get<1>(g2l)(0));
  msg += ",";
  msg += std::to_string(std::get<1>(g2l)(1));
  msg += ",";
  msg += std::to_string(std::get<1>(g2l)(2));
  msg += ")";

  return msg;
}

std::string local2global_test(const PatchTransformations &pt,
                              const MultiPatch::Cake::svec_u &x) {
  using MultiPatch::Cake::patch_piece;
  std::string msg{"has returnd "};

  const auto l2g = pt.local2global(pt, piece_idx(patch_piece::plus_z), x);

  msg += "(";
  msg += std::to_string(l2g(0));
  msg += ",";
  msg += std::to_string(l2g(1));
  msg += ",";
  msg += std::to_string(l2g(2));
  msg += ")";

  return msg;
}

std::string glg_test(const PatchTransformations &pt,
                     const MultiPatch::Cake::svec_u &original_x) {
  std::string msg{"has returnd "};

  const auto g2l = pt.global2local(pt, original_x);

  MultiPatch::Cake::svec_u local_coords = {
      (std::get<1>(g2l))(0), (std::get<1>(g2l))(1), (std::get<1>(g2l))(2)};

  const auto owner_patch = MultiPatch::Cake::get_owner_patch(pt, original_x);

  const auto l2g = pt.local2global(pt, piece_idx(owner_patch), local_coords);

  msg += "(";
  msg += std::to_string(l2g(0));
  msg += ",";
  msg += std::to_string(l2g(1));
  msg += ",";
  msg += std::to_string(l2g(2));
  msg += ")";

  return msg;
}

} // namespace CakeTests
} // namespace MultiPatch

extern "C" void run_cake_tests() {
  DECLARE_CCTK_PARAMETERS

  using MultiPatch::Cake::patch_piece;
  using MultiPatch::Cake::svec_u;

  using std::mt19937_64;
  using std::uniform_real_distribution;

  using namespace MultiPatch::CakeTests;

  using MultiPatch::dim;
  using MultiPatch::SetupCake;
  using MultiPatchTests::random_seed;

  const auto ps = SetupCake();
  const auto &pt = ps.transformations;

  mt19937_64 engine(random_seed);

  const auto r0 = pt.cake_inner_boundary_radius;
  const auto r1 = pt.cake_outer_boundary_radius;

  const auto cube_midpoint = r0 / 2.0;
  const auto patch_midpoint = r1 / 2.0;

  uniform_real_distribution<CCTK_REAL> x_distrib(-r1, r1);
  uniform_real_distribution<CCTK_REAL> y_distrib(-r1, r1);
  uniform_real_distribution<CCTK_REAL> z_distrib(-r1, r1);

  uniform_real_distribution<CCTK_REAL> local_distrib(-1.0, 1.0);

  auto global_point = svec_u{0.0, 0.0, 0.0};
  auto local_point = svec_u{0.0, 0.0, 0.0};

  CCTK_INFO("Running cake patch tests:");

  // The patch_owner tests do not take random numbers as input because we
  // need to compare results with known values.
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

  for (const auto &data : owner_test_data) {
    CCTK_VINFO("  Patch owner test at point (%f, %f, %f) %s", data.first(0),
               data.first(1), data.first(2),
               patch_owner(pt, data.first, data.second).c_str());
  }

  // These tests take random numbers as input
  for (int i = 0; i < repeat_tests; i++) {
    global_point = {x_distrib(engine), y_distrib(engine), z_distrib(engine)};
    local_point = {local_distrib(engine), local_distrib(engine),
                   local_distrib(engine)};

    CCTK_VINFO("  Global -> local -> global transformation test at point (%f, "
               "%f, %f) %s",
               local_point(0), local_point(1), local_point(2),
               glg_test(pt, global_point).c_str());
  }
}