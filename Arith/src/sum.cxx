#include "sum.hxx"

#include <cctk.h>

#include <functional>

namespace Arith {
using namespace std;

namespace {
constexpr CCTK_REAL f0(int i) { return 0; }
constexpr CCTK_REAL f1(int i) { return 1; }
constexpr CCTK_REAL fi(int i) { return i; }

constexpr CCTK_REAL g0(int i, int j) { return 0; }
constexpr CCTK_REAL g1(int i, int j) { return 1; }
constexpr CCTK_REAL gi(int i, int j) { return i; }
constexpr CCTK_REAL gj(int i, int j) { return j; }
constexpr CCTK_REAL gij(int i, int j) { return i + 10 * j; }

constexpr CCTK_REAL h0(int i, int j, int k) { return 0; }
constexpr CCTK_REAL h1(int i, int j, int k) { return 1; }
constexpr CCTK_REAL hi(int i, int j, int k) { return i; }
constexpr CCTK_REAL hj(int i, int j, int k) { return j; }
constexpr CCTK_REAL hk(int i, int j, int k) { return j; }
constexpr CCTK_REAL hijk(int i, int j, int k) { return i + 10 * j + 100 * k; }
} // namespace

// This function is compiled, but not executed. The tests are "run" at
// compile time. If this function compiles, the tests pass.
void TestSum() {
  constexpr equal_to<CCTK_REAL> eq;

  static_assert(eq(sum<3>(f0), 0.0));
  static_assert(eq(sum<3>(f1), 3.0));
  static_assert(eq(sum<3>(fi), 3.0));

  static_assert(eq(sum<4>(f0), 0.0));
  static_assert(eq(sum<4>(f1), 4.0));
  static_assert(eq(sum<4>(fi), 6.0));

  static_assert(eq(sum<3>(g0), 0.0));
  static_assert(eq(sum<3>(g1), 9.0));
  static_assert(eq(sum<3>(gi), 9.0));
  static_assert(eq(sum<3>(gj), 9.0));
  static_assert(eq(sum<3>(gij), 99.0));

  static_assert(eq(sum<3>(h0), 0.0));
  static_assert(eq(sum<3>(h1), 27.0));
  static_assert(eq(sum<3>(hi), 27.0));
  static_assert(eq(sum<3>(hj), 27.0));
  static_assert(eq(sum<3>(hk), 27.0));
  static_assert(eq(sum<3>(hijk), 2997.0));
}

} // namespace Arith
