#include "vec.hxx"

#include <cctk.h>

#include <functional>

namespace Arith {
using namespace std;

// This function is compiled, but not executed. The tests are "run" at
// compile time. If this function compiles, the tests pass.
void TestVec() {
  constexpr equal_to<CCTK_REAL> eq;

  using V3D = vec<CCTK_REAL, 3, DN>;
  constexpr equal_to<V3D> eqv;

  static_assert(eq(vec<CCTK_REAL, 3, DN>()(0), 0));
  static_assert(eq(vec<CCTK_REAL, 3, DN>()(1), 0));
  static_assert(eq(vec<CCTK_REAL, 3, DN>()(2), 0));
  static_assert(eq(vec<CCTK_REAL, 3, UP>()(0), 0));
  static_assert(eq(vec<CCTK_REAL, 3, UP>()(1), 0));
  static_assert(eq(vec<CCTK_REAL, 3, UP>()(2), 0));

  static_assert(eq(vec<CCTK_REAL, 3, DN>{1, 2, 3}(0), 1));
  static_assert(eq(vec<CCTK_REAL, 3, DN>{1, 2, 3}(1), 2));
  static_assert(eq(vec<CCTK_REAL, 3, DN>{1, 2, 3}(2), 3));

  static_assert(eqv(vec<CCTK_REAL, 3, DN>::iota(), {0, 1, 2}));

  static_assert(eqv(vec<CCTK_REAL, 3, DN>::unit(0), {1, 0, 0}));
  static_assert(eqv(vec<CCTK_REAL, 3, DN>::unit(1), {0, 1, 0}));
  static_assert(eqv(vec<CCTK_REAL, 3, DN>::unit(2), {0, 0, 1}));

  constexpr V3D x = {5, 6, 8};
  static_assert(eqv(+x, {5, 6, 8}));
  static_assert(eqv(-x, {-5, -6, -8}));

  constexpr V3D y = {10, 4, 6};
  static_assert(eqv(x + y, {15, 10, 14}));
  static_assert(eqv(x - y, {-5, 2, 2}));

  static_assert(eqv(3 * x, {15, 18, 24}));
  static_assert(eqv(x * 3, {15, 18, 24}));
}

} // namespace Arith
