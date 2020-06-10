#ifndef RATIONAL_HXX
#define RATIONAL_HXX

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <type_traits>

namespace Arith {
using namespace std;

template <typename I = int64_t> struct rational {
  I num, den;

  struct no_normalize {};

  constexpr void normalize() {
    const I x = gcd(num, den);
    num /= x;
    den /= x;
  }

  constexpr rational() : num(0), den(1) {}

  constexpr rational(const rational &) = default;
  constexpr rational(rational &&) = default;
  constexpr rational &operator=(const rational &) = default;
  constexpr rational &operator=(rational &&) = default;

  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  constexpr rational(const J &i) : num(i), den(1) {}
  template <typename J, typename K, enable_if_t<is_integral_v<J> > * = nullptr,
            enable_if_t<is_integral_v<K> > * = nullptr>
  constexpr rational(const J &num, const K &den, no_normalize)
      : num(num), den(den) {}
  template <typename J, typename K, enable_if_t<is_integral_v<J> > * = nullptr,
            enable_if_t<is_integral_v<K> > * = nullptr>
  constexpr rational(const J &num, const K &den) : num(num), den(den) {
    if (this->den < 0) {
      this->num = -this->num;
      this->den = -this->den;
    }
    normalize();
  }

  template <typename F, enable_if_t<is_floating_point_v<F> > * = nullptr>
  constexpr rational(const F &f) : num(0), den(1) {
    const F f1 = nextafter(f);
    const F df = f1 - f;
    const F df1 = 1 / df;
    assert(df1 == rint(df1));
    den = llrint(df1);
    num = llrint(f * df1);
    normalize();
  }

  template <typename F, enable_if_t<is_floating_point_v<F> > * = nullptr>
  constexpr operator F() const {
    return F(num) / F(den);
  }

  friend constexpr rational operator+(const rational &x) {
    return rational(+x.num, x.den, no_normalize());
  }
  friend constexpr rational operator-(const rational &x) {
    return rational(-x.num, x.den, no_normalize());
  }

  friend constexpr rational operator+(const rational &x, const rational &y) {
    return rational(x.num * y.den + x.den * y.num, x.den * y.den);
  }
  friend constexpr rational operator-(const rational &x, const rational &y) {
    return rational(x.num * y.den - x.den * y.num, x.den * y.den);
  }
  friend constexpr rational operator*(const rational &x, const rational &y) {
    return rational(x.num * y.num, x.den * y.den);
  }
  friend constexpr rational operator/(const rational &x, const rational &y) {
    return rational(x.num * y.den, x.den * y.num);
  }

  constexpr rational &operator+=(const rational &x) {
    return *this = *this + x;
  }
  constexpr rational &operator-=(const rational &x) {
    return *this = *this - x;
  }
  constexpr rational &operator*=(const rational &x) {
    return *this = *this * x;
  }
  constexpr rational &operator/=(const rational &x) {
    return *this = *this / x;
  }

  friend constexpr bool operator==(const rational &x, const rational &y) {
    return x.num * y.den == x.den * y.num;
  }
  friend constexpr bool operator!=(const rational &x, const rational &y) {
    return !(x == y);
  }
  friend constexpr bool operator<(const rational &x, const rational &y) {
    return x.num * y.den < x.den * y.num;
  }
  friend constexpr bool operator>(const rational &x, const rational &y) {
    return y < x;
  }
  friend constexpr bool operator<=(const rational &x, const rational &y) {
    return !(x > y);
  }
  friend constexpr bool operator>=(const rational &x, const rational &y) {
    return y <= x;
  }

  friend ostream &operator<<(ostream &os, const rational &x) {
    return os << "(" << x.num << "/" << x.den << ")";
  }
};

template <typename I> constexpr rational<I> abs(const rational<I> &x) {
  // std::abs is not constexpr
  return rational<I>(x.num >= 0 ? x.num : -x.num, x.den);
}

template <typename I>
constexpr rational<I> max(const rational<I> &x, const rational<I> &y) {
  return x >= y ? x : y;
}

template <typename I>
constexpr rational<I> min(const rational<I> &x, const rational<I> &y) {
  return x <= y ? x : y;
}

template <typename I, typename J, enable_if_t<is_integral_v<J> > * = nullptr>
constexpr rational<I> pow(rational<I> x, J a) {
  if (a < 0) {
    x = 1 / x;
    a = -a;
  }
  rational<I> r = 1;
  while (a > 0) {
    if (a % 2)
      r *= x;
    x *= x;
    a /= 2;
  }
  return r;
}
template <typename I>
constexpr double pow(const rational<I> &x, const rational<I> &y) {
  return std::pow(double(x), double(y));
}

} // namespace Arith

#endif // #ifndef RATIONAL_HXX
