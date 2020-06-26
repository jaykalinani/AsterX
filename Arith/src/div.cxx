#include "div.hxx"

namespace Arith {

template <typename T> constexpr bool test_div(T x, T y) {
  T d = div(x, y);
  T m = mod(x, y);
  return d * y + m == x && m >= 0 && m < abs(y);
}

// This function is compiled, but not executed. The tests are "run" at
// compile time. If this function compiles, the tests pass.
void TestDiv() {
  static_assert(test_div(0, 1));
  static_assert(test_div(0, 2));
  static_assert(test_div(0, 3));

  static_assert(test_div(1, 1));
  static_assert(test_div(1, 2));
  static_assert(test_div(1, 3));
  static_assert(test_div(2, 1));
  static_assert(test_div(2, 2));
  static_assert(test_div(2, 3));
  static_assert(test_div(3, 1));
  static_assert(test_div(3, 2));
  static_assert(test_div(3, 3));

  static_assert(test_div(-1, 1));
  static_assert(test_div(-1, 2));
  static_assert(test_div(-1, 3));
  static_assert(test_div(-2, 1));
  static_assert(test_div(-2, 2));
  static_assert(test_div(-2, 3));
  static_assert(test_div(-3, 1));
  static_assert(test_div(-3, 2));
  static_assert(test_div(-3, 3));
}

} // namespace Arith
