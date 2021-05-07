#include "vec.hxx"

#include <cctk.h>

#include <cmath>
#include <functional>

namespace Arith {
using namespace std;

template <typename T, typename U> constexpr bool eq(const T &x, const U &y) {
  using std::isnan;
  return equal_to<CCTK_REAL>()(x, y) || (isnan(x) && isnan(y));
}

template <typename T, int D, dnup_t dnup>
constexpr bool eqv(const vec<T, D, dnup> &x, const vec<T, D, dnup> &y) {
  using std::isnan;
  return equal_to<vec<T, D, dnup> >()(x, y) || (isnan(x) && isnan(y));
}

// This function is compiled, but not executed. The tests are "run" at
// compile time. If this function compiles, the tests pass.
void TestVec() {
  // nvcc V11.1.74 doesn't accept this as "constexpr" values
#ifndef __CUDACC__
  using V3D = vec<CCTK_REAL, 3, DN>;

  constexpr CCTK_REAL N = nan<CCTK_REAL>();
  static_assert(eq(vec<CCTK_REAL, 3, DN>()(0), N));
  static_assert(eq(vec<CCTK_REAL, 3, DN>()(1), N));
  static_assert(eq(vec<CCTK_REAL, 3, DN>()(2), N));
  static_assert(eq(vec<CCTK_REAL, 3, UP>()(0), N));
  static_assert(eq(vec<CCTK_REAL, 3, UP>()(1), N));
  static_assert(eq(vec<CCTK_REAL, 3, UP>()(2), N));

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
#endif
}

} // namespace Arith
