#include "tests.hxx"
#include "multipatch.hxx"

#include <cctk_Parameters.h>
#include <string>

namespace MultiPatch {
namespace CartesianTests {

std::string local2global(const PatchTransformations &pt,
                         const vec<CCTK_REAL, dim, UP> &x) {

  using MultiPatchTests::colored;
  using MultiPatchTests::isapprox;
  using MultiPatchTests::string_color;

  std::string msg{"has "};

  const auto l2g = pt.local2global(pt, 0, x);
  const bool eq_x = isapprox(l2g(0), x(0));
  const bool eq_y = isapprox(l2g(1), x(1));
  const bool eq_z = isapprox(l2g(2), x(2));

  if (eq_x && eq_y && eq_z) {
    msg += colored<string_color::green>("PASSED");
  } else {
    msg += colored<string_color::green>("FAILED");

    if (!eq_x) {
      msg += " The result in the x direction is ";
      msg += std::to_string(l2g(0));
      msg += ".";
    }

    if (!eq_y) {
      msg += " The result in the y direction is ";
      msg += std::to_string(l2g(1));
      msg += ".";
    }

    if (!eq_z) {
      msg += " The result in the z direction is ";
      msg += std::to_string(l2g(02));
      msg += ".";
    }
  }

  return msg;
}

std::string global2local(const PatchTransformations &pt,
                         const vec<CCTK_REAL, dim, UP> &x) {

  using MultiPatchTests::isapprox;

  std::string msg{"has "};

  const auto g2l = pt.global2local(pt, x);
  const bool eq_x = isapprox(std::get<1>(g2l)(0), x(0));
  const bool eq_y = isapprox(std::get<1>(g2l)(1), x(1));
  const bool eq_z = isapprox(std::get<1>(g2l)(2), x(2));

  if (eq_x && eq_y && eq_z) {
    msg += "\033[32;1mPASSED\033[0m.";
  } else {
    msg += "\033[31;1mFAILED\033[0m. Reason(s):";

    if (!eq_x) {
      msg += " The result in the x direction is ";
      msg += std::to_string(std::get<1>(g2l)(0));
      msg += ".";
    }

    if (!eq_y) {
      msg += " The result in the y direction is ";
      msg += std::to_string(std::get<1>(g2l)(1));
      msg += ".";
    }

    if (!eq_z) {
      msg += " The result in the z direction is ";
      msg += std::to_string(std::get<1>(g2l)(02));
      msg += ".";
    }
  }

  return msg;
}

} // namespace CartesianTests

} // namespace MultiPatch

extern "C" void run_cartesian_tests() {
  DECLARE_CCTK_PARAMETERS

  using std::mt19937_64;
  using std::uniform_real_distribution;

  using namespace MultiPatch::CartesianTests;

  using MultiPatch::dim;
  using MultiPatch::SetupCartesian;
  using MultiPatchTests::random_seed;

  using Arith::UP;
  using Arith::vec;

  const auto ps = SetupCartesian();
  const auto &pt = ps.transformations;

  mt19937_64 engine(random_seed);
  uniform_real_distribution<CCTK_REAL> x_distrib(pt.cartesian_xmin,
                                                 pt.cartesian_xmax);
  uniform_real_distribution<CCTK_REAL> y_distrib(pt.cartesian_ymin,
                                                 pt.cartesian_ymax);
  uniform_real_distribution<CCTK_REAL> z_distrib(pt.cartesian_zmin,
                                                 pt.cartesian_zmax);
  auto point = vec<CCTK_REAL, dim, UP>{0.0, 0.0, 0.0};

  CCTK_INFO("Running cartesian patch tests:");

  for (int i = 0; i < repeat_tests; i++) {
    point = {x_distrib(engine), y_distrib(engine), z_distrib(engine)};
    CCTK_VINFO("  Local to global transformation test at point (%f, %f, %f) %s",
               point(0), point(1), point(2), local2global(pt, point).c_str());
  }

  for (int i = 0; i < repeat_tests; i++) {
    point = {x_distrib(engine), y_distrib(engine), z_distrib(engine)};
    CCTK_VINFO("  Global to local transformation test at point (%f, %f, %f) %s",
               point(0), point(1), point(2), global2local(pt, point).c_str());
  }
}