#include "mat.hxx"

#include <cctk.h>

#include <functional>

namespace Arith {
using namespace std;

// This function is compiled, but not executed. The tests are "run" at
// compile time. If this function compiles, the tests pass.
void TestMat() {
  constexpr equal_to<CCTK_REAL> eq;

  using M3D = mat<CCTK_REAL, 3, DN, DN>;
  constexpr equal_to<M3D> eqm;

  static_assert(eq(mat<CCTK_REAL, 3, DN, DN>()(0, 0), 0));
  static_assert(eq(mat<CCTK_REAL, 3, DN, DN>()(0, 1), 0));
  static_assert(eq(mat<CCTK_REAL, 3, DN, DN>()(0, 2), 0));
  static_assert(eq(mat<CCTK_REAL, 3, UP, UP>()(1, 1), 0));
  static_assert(eq(mat<CCTK_REAL, 3, UP, UP>()(1, 2), 0));
  static_assert(eq(mat<CCTK_REAL, 3, UP, UP>()(2, 2), 0));

  static_assert(eq(mat<CCTK_REAL, 3, DN, DN>{1, 2, 3, 4, 5, 6}(0, 0), 1));
  static_assert(eq(mat<CCTK_REAL, 3, DN, DN>{1, 2, 3, 4, 5, 6}(0, 1), 2));
  static_assert(eq(mat<CCTK_REAL, 3, DN, DN>{1, 2, 3, 4, 5, 6}(0, 2), 3));
  static_assert(eq(mat<CCTK_REAL, 3, DN, DN>{1, 2, 3, 4, 5, 6}(1, 0), 2));
  static_assert(eq(mat<CCTK_REAL, 3, DN, DN>{1, 2, 3, 4, 5, 6}(1, 1), 4));
  static_assert(eq(mat<CCTK_REAL, 3, DN, DN>{1, 2, 3, 4, 5, 6}(1, 2), 5));
  static_assert(eq(mat<CCTK_REAL, 3, DN, DN>{1, 2, 3, 4, 5, 6}(2, 0), 3));
  static_assert(eq(mat<CCTK_REAL, 3, DN, DN>{1, 2, 3, 4, 5, 6}(2, 1), 5));
  static_assert(eq(mat<CCTK_REAL, 3, DN, DN>{1, 2, 3, 4, 5, 6}(2, 2), 6));

  static_assert(eqm(mat<CCTK_REAL, 3, DN, DN>::iota1(), {0, 0, 0, 1, 1, 2}));
  static_assert(eqm(mat<CCTK_REAL, 3, DN, DN>::iota2(), {0, 1, 2, 1, 2, 2}));

  static_assert(eqm(mat<CCTK_REAL, 3, DN, DN>::unit(0, 0), {1, 0, 0, 0, 0, 0}));
  static_assert(eqm(mat<CCTK_REAL, 3, DN, DN>::unit(0, 1), {0, 1, 0, 0, 0, 0}));
  static_assert(eqm(mat<CCTK_REAL, 3, DN, DN>::unit(0, 2), {0, 0, 1, 0, 0, 0}));
  static_assert(eqm(mat<CCTK_REAL, 3, DN, DN>::unit(1, 0), {0, 1, 0, 0, 0, 0}));
  static_assert(eqm(mat<CCTK_REAL, 3, DN, DN>::unit(1, 1), {0, 0, 0, 1, 0, 0}));
  static_assert(eqm(mat<CCTK_REAL, 3, DN, DN>::unit(1, 2), {0, 0, 0, 0, 1, 0}));
  static_assert(eqm(mat<CCTK_REAL, 3, DN, DN>::unit(2, 0), {0, 0, 1, 0, 0, 0}));
  static_assert(eqm(mat<CCTK_REAL, 3, DN, DN>::unit(2, 1), {0, 0, 0, 0, 1, 0}));
  static_assert(eqm(mat<CCTK_REAL, 3, DN, DN>::unit(2, 2), {0, 0, 0, 0, 0, 1}));

  constexpr M3D x = {71, 90, 14, 50, 41, 65};
  static_assert(eqm(+x, {71, 90, 14, 50, 41, 65}));
  static_assert(eqm(-x, {-71, -90, -14, -50, -41, -65}));

  constexpr M3D y = {77, 32, 67, 5, 81, 10};
  static_assert(eqm(x + y, {148, 122, 81, 55, 122, 75}));
  static_assert(eqm(x - y, {-6, 58, -53, 45, -40, 55}));
  static_assert(eqm(3 * x, {213, 270, 42, 150, 123, 195}));
  static_assert(eqm(x * 3, {213, 270, 42, 150, 123, 195}));
}

} // namespace Arith
