#include "rational.hxx"

#include <numeric>

namespace Arith {
using namespace std;

template <typename I> constexpr bool isgood(const rational<I> &x) {
  return x.den >= 0 && gcd(x.num, x.den) == 1;
}

// This function is compiled, but not executed. The tests are "run" at
// compile time. If this function compiles, the tests pass.
void TestRational() {
  static_assert(rational().num == 0);
  static_assert(rational(1).num == 1);
  static_assert(rational(1, 2).num == 1);
  static_assert(rational(2, 4).num == 1);
  static_assert(rational().den == 1);
  static_assert(rational(1).den == 1);
  static_assert(rational(1, 2).den == 2);
  static_assert(rational(2, 4).den == 2);

  constexpr auto n = rational();
  constexpr auto e = rational(1);
  constexpr auto x = rational(1, 2);
  constexpr auto y = rational(2, 3);
  constexpr auto z = rational(3, 4);
  constexpr auto a = 4;

  +x;
  -x;

  x + y;
  x - y;
  x *y;
  x / y;

  a + y;
  a - y;
  a *y;
  a / y;

  x + a;
  x - a;
  x *a;
  x / a;

  auto r = x;
  r = x;
  r += y;
  r = x;
  r -= y;
  r = x;
  r *= y;
  r = x;
  r /= y;

  r = x;
  r += a;
  r = x;
  r -= a;
  r = x;
  r *= a;
  r = x;
  r /= a;

  x == y;
  x != y;
  x < y;
  x > y;
  x <= y;
  x >= y;

  a == y;
  a != y;
  a < y;
  a > y;
  a <= y;
  a >= y;

  x == a;
  x != a;
  x < a;
  x > a;
  x <= a;
  x >= a;

  static_assert(isgood(abs(x)));

  static_assert(isgood(max(x, y)));
  static_assert(isgood(min(x, y)));

  static_assert(isgood(pow(x, a)));

  static_assert(+x == n + x);
  static_assert(x + n == x);
  static_assert(n + x == x);
  static_assert(x + (y + z) == (x + y) + z);
  static_assert(x + y == y + x);

  static_assert(-x == n - x);
  static_assert(-(-x) == x);
  static_assert(x + (-y) == x - y);
  static_assert(x - x == n);

  static_assert(x * e == x);
  static_assert(e * x == x);
  static_assert(x * (y * z) == (x * y) * z);
  static_assert(x * y == y * x);
  static_assert(n * x == n);

  static_assert(1 / (1 / x) == x);
  static_assert(x / x == e);

  static_assert(x * (y + z) == x * y + x * z);
}

} // namespace Arith
