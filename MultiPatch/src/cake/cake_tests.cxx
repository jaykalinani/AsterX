#include "cake.hxx"
#include "tests.hxx"

#include <cctk_Parameters.h>

namespace MultiPatch {
namespace CakeTests {

std::string patch_owner(const PatchTransformations &pt,
                        const MultiPatch::Cake::svec_u &x) {
  using namespace MultiPatch::Cake;

  std::string msg{"has returnd "};

  const auto owner_patch = get_owner_patch(pt, x);

  return piece_name(owner_patch);
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

  const auto max = pt.cake_outer_boundary_radius;
  const auto min = 0.0;

  uniform_real_distribution<CCTK_REAL> x_distrib(min, max);
  uniform_real_distribution<CCTK_REAL> y_distrib(min, max);
  uniform_real_distribution<CCTK_REAL> z_distrib(min, max);

  uniform_real_distribution<CCTK_REAL> local_distrib(-1.0, 1.0);

  auto global_point = svec_u{0.0, 0.0, 0.0};
  auto local_point = svec_u{0.0, 0.0, 0.0};

  CCTK_INFO("Running cake patch tests:");

  for (int i = 0; i < repeat_tests; i++) {
    global_point = {x_distrib(engine), y_distrib(engine), z_distrib(engine)};
    local_point = {local_distrib(engine), local_distrib(engine),
                   local_distrib(engine)};

    CCTK_VINFO("  Patch owner test at point (%f, %f, %f) %s", global_point(0),
               global_point(1), global_point(2),
               patch_owner(pt, global_point).c_str());

    CCTK_VINFO("  Global to local transformation test at point (%f, %f, %f) %s",
               global_point(0), global_point(1), global_point(2),
               global2local_test(pt, global_point).c_str());

    CCTK_VINFO("  Local to global transformation test at point (%f, %f, %f) %s",
               local_point(0), local_point(1), local_point(2),
               local2global_test(pt, local_point).c_str());

    CCTK_VINFO("  Global -> local -> global transformation test at point (%f, "
               "%f, %f) %s",
               local_point(0), local_point(1), local_point(2),
               glg_test(pt, global_point).c_str());
  }
}