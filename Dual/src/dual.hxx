#include <cassert>
#include <cmath>
#include <functional>
#include <ostream>

namespace Dual {
using namespace std;

template <typename T> struct dual {
  T val, eps;

  dual(const dual &) = default;
  dual(dual &&) = default;
  dual &operator=(const dual &) = default;
  dual &operator=(dual &&) = default;

  constexpr dual() : val(), eps() {}
  constexpr dual(const T &x) : val(x), eps() {}
  constexpr dual(const T &x, const T &y) : val(x), eps(y) {}

  friend constexpr dual operator+(const dual &x) { return {+x.val, +x.eps}; }
  friend constexpr dual operator-(const dual &x) { return {-x.val, -x.eps}; }

  friend constexpr dual operator+(const dual &x, const dual &y) {
    return {x.val + y.val, x.eps + y.eps};
  }
  friend constexpr dual operator-(const dual &x, const dual &y) {
    return {x.val - y.val, x.eps - y.eps};
  }
  friend constexpr dual operator+(const dual &x, const T &y) {
    return {x.val + y, x.eps};
  }
  friend constexpr dual operator-(const dual &x, const T &y) {
    return {x.val - y, x.eps};
  }
  friend constexpr dual operator*(const dual &x, const T &y) {
    return {x.val * y, x.eps * y};
  }
  friend constexpr dual operator*(const T &x, const dual &y) {
    return {x * y.val, x * y.eps};
  }
  friend constexpr dual operator/(const dual &x, const T &y) {
    return {x.val / y, x.eps / y};
  }
  friend constexpr dual operator/(const dual &x, const dual &y) {
    return {x.val / y.val, (x.eps * y.val - x.val * y.eps) / (y.val * y.val)};
  }

  friend constexpr dual operator*(const dual &x, const dual &y) {
    return {x.val * y.val, x.val * y.eps + x.eps * y.val};
  }

  dual &operator+=(const dual &x) { return *this = *this + x; }
  dual &operator-=(const dual &x) { return *this = *this - x; }
  dual &operator*=(const T &x) { return *this = *this * x; }
  dual &operator/=(const T &x) { return *this = *this / x; }
  dual &operator*=(const dual &x) { return *this = *this * x; }
  dual &operator/=(const dual &x) { return *this = *this / x; }

  friend constexpr bool operator==(const dual &x, const dual &y) {
    return x.val == y.val;
  };
  friend constexpr bool operator<(const dual &x, const dual &y) {
    return x.val < y.val;
  };

  friend constexpr bool operator!=(const dual &x, const dual &y) {
    return !(x == y);
  }
  friend constexpr bool operator>(const dual &x, const dual &y) {
    return y < x;
  }
  friend constexpr bool operator<=(const dual &x, const dual &y) {
    return !(x > y);
  }
  friend constexpr bool operator>=(const dual &x, const dual &y) {
    return !(x < y);
  }

  friend ostream &operator<<(ostream &os, const dual &x) {
    return os << x.val << "+eps*" << x.val;
  }
};

} // namespace Dual
namespace std {
template <typename T> using dual = Dual::dual<T>;

template <typename T> struct equal_to<dual<T> > {
  constexpr bool operator()(const dual<T> &x, const dual<T> &y) const {
    return equal_to<T>()(x.val, y.val) && equal_to<T>()(x.eps, y.eps);
  }
};

template <typename T> struct less<dual<T> > {
  constexpr bool operator()(const dual<T> &x, const dual<T> &y) const {
    return less<T>(x.val, y.val) ||
           (equal_to<T>(x.val, y.val) && less<T>(x.eps, y.eps));
  }
};

template <typename T> constexpr dual<T> cos(const dual<T> &x);
template <typename T> constexpr dual<T> exp(const dual<T> &x);
template <typename T> constexpr dual<T> fabs(const dual<T> &x);
template <typename T> constexpr dual<T> pow2(const dual<T> &x);
template <typename T> constexpr dual<T> sin(const dual<T> &x);
template <typename T> constexpr dual<T> sqrt(const dual<T> &x);

template <typename T> constexpr dual<T> cos(const dual<T> &x) {
  return {cos(x.val), -sin(x.val) * x.eps};
}

template <typename T> constexpr dual<T> exp(const dual<T> &x) {
  auto r = exp(x.val);
  return {r, r * x.eps};
}

template <typename T> constexpr dual<T> fabs(const dual<T> &x) {
  return sqrt(pow2(x));
}

template <typename T> constexpr dual<T> pow2(const dual<T> &x) { return x * x; }

template <typename T> constexpr dual<T> sin(const dual<T> &x) {
  return {sin(x.val), cos(x.val) * x.eps};
}

template <typename T> constexpr dual<T> sqrt(const dual<T> &x) {
  auto r = sqrt(x.val);
  return {r, x.eps / (2 * r)};
}
} // namespace std
