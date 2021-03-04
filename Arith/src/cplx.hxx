#ifndef CPLX_HXX
#define CPLX_HXX

#ifdef __CUDACC__
#define CPLX_DEVICE __device__
#define CPLX_HOST __host__
#else
#define CPLX_DEVICE
#define CPLX_HOST
#endif

#include "vect.hxx"

#include <cassert>
#include <cmath>
#include <functional>
#include <ostream>

namespace Arith {
using namespace std;

template <typename T> struct cplx;

template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST Arith::cplx<T>
abs(const Arith::cplx<T> &x);
template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST Arith::cplx<T>
cbrt(const Arith::cplx<T> &x);
template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST Arith::cplx<T>
conj(const Arith::cplx<T> &x);
template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST Arith::cplx<T>
cos(const Arith::cplx<T> &x);
template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST Arith::cplx<T>
exp(const Arith::cplx<T> &x);
template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST T
imag(const Arith::cplx<T> &x);
template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST bool
isnan(const Arith::cplx<T> &x);
template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST T
norm(const Arith::cplx<T> &x);
template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST T
real(const Arith::cplx<T> &x);
template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST Arith::cplx<T>
pow(Arith::cplx<T> x, int n);
template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST Arith::cplx<T>
pow2(const Arith::cplx<T> &x);
template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST Arith::cplx<T>
sin(const Arith::cplx<T> &x);
template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST Arith::cplx<T>
sqrt(const Arith::cplx<T> &x);

template <typename T> struct cplx {
  T real, imag;

  CCTK_ATTRIBUTE_ALWAYS_INLINE cplx(const cplx &) = default;
  CCTK_ATTRIBUTE_ALWAYS_INLINE cplx(cplx &&) = default;
  CCTK_ATTRIBUTE_ALWAYS_INLINE cplx &operator=(const cplx &) = default;
  CCTK_ATTRIBUTE_ALWAYS_INLINE cplx &operator=(cplx &&) = default;

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST cplx()
      : real(), imag() {}
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST cplx(const T &x)
      : real(x), imag() {}
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST cplx(const T &x,
                                                                    const T &y)
      : real(x), imag(y) {}

  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST cplx
  operator+(const cplx &x) {
    return {+x.real, +x.imag};
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST cplx
  operator-(const cplx &x) {
    return {-x.real, -x.imag};
  }

  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST cplx
  operator+(const cplx &x, const cplx &y) {
    return {x.real + y.real, x.imag + y.imag};
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST cplx
  operator-(const cplx &x, const cplx &y) {
    return {x.real - y.real, x.imag - y.imag};
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST cplx
  operator+(const cplx &x, const T &y) {
    return {x.real + y, x.imag};
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST cplx
  operator-(const cplx &x, const T &y) {
    return {x.real - y, x.imag};
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST cplx
  operator*(const cplx &x, const T &y) {
    return {x.real * y, x.imag * y};
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST cplx
  operator*(const T &x, const cplx &y) {
    return {x * y.real, x * y.imag};
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST cplx
  operator/(const cplx &x, const T &y) {
    return {x.real / y, x.imag / y};
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST cplx
  operator*(const cplx &x, const cplx &y) {
    return {x.real * y.real - x.imag * y.imag,
            x.real * y.imag + x.imag * y.real};
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST cplx
  operator/(const cplx &x, const cplx &y) {
    return x * conj(y) / norm(y);
  }

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST cplx &
  operator+=(const cplx &x) {
    return *this = *this + x;
  }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST cplx &
  operator-=(const cplx &x) {
    return *this = *this - x;
  }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST cplx &
  operator*=(const T &x) {
    return *this = *this * x;
  }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST cplx &
  operator/=(const T &x) {
    return *this = *this / x;
  }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST cplx &
  operator*=(const cplx &x) {
    return *this = *this * x;
  }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST cplx &
  operator/=(const cplx &x) {
    return *this = *this / x;
  }

  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST bool
  operator==(const cplx &x, const cplx &y) {
    return x.real == y.real;
  };
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST bool
  operator<(const cplx &x, const cplx &y) {
    return x.real < y.real;
  };

  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST bool
  operator!=(const cplx &x, const cplx &y) {
    return !(x == y);
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST bool
  operator>(const cplx &x, const cplx &y) {
    return y < x;
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST bool
  operator<=(const cplx &x, const cplx &y) {
    return !(x > y);
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST bool
  operator>=(const cplx &x, const cplx &y) {
    return !(x < y);
  }

  friend ostream &operator<<(ostream &os, const cplx &x) {
    return os << x.real << "+" << x.real << "*i";
  }
};

template <typename T> struct zero<cplx<T> > {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST cplx<T>
  operator()() const {
    return cplx<T>(zero<T>()(), zero<T>()());
  }
};

template <typename T> struct one<cplx<T> > {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST cplx<T>
  operator()() const {
    return cplx<T>(one<T>()(), zero<T>()());
  }
};

} // namespace Arith
namespace std {
template <typename T> struct equal_to<Arith::cplx<T> > {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST bool
  operator()(const Arith::cplx<T> &x, const Arith::cplx<T> &y) const {
    return equal_to<T>()(x.real, y.real) && equal_to<T>()(x.imag, y.imag);
  }
};

template <typename T> struct less<Arith::cplx<T> > {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST bool
  operator()(const Arith::cplx<T> &x, const Arith::cplx<T> &y) const {
    if (less<T>(x.real, y.real))
      return true;
    if (less<T>(y.real, x.real))
      return false;
    if (less<T>(x.imag, y.imag))
      return true;
    if (less<T>(y.imag, x.imag))
      return false;
    return false;
  }
};
} // namespace std

namespace Arith {

template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST Arith::cplx<T>
abs(const Arith::cplx<T> &x) {
  return sqrt(norm(x));
}

// template <typename T> constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE
// CPLX_HOST Arith::cplx<T> cbrt(const Arith::cplx<T> &x)
// {
//   using std::cbrt;
//   const T r = cbrt(x.real);
//   return {r, r / (3 * x.real) * x.imag};
// }

template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST Arith::cplx<T>
conj(const Arith::cplx<T> &x) {
  return {x.real, -x.imag};
}

// template <typename T> constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE
// CPLX_HOST Arith::cplx<T> cos(const Arith::cplx<T> &x) {
//   using std::cos, std::sin;
//   return {cos(x.real), -sin(x.real) * x.imag};
// }

// template <typename T> constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE
// CPLX_HOST Arith::cplx<T> exp(const Arith::cplx<T> &x) {
//   using std::exp;
//   const T r = exp(x.real);
//   return {r, r * x.imag};
// }

template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST T
imag(const Arith::cplx<T> &x) {
  return x.imag;
}

template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST bool
isnan(const Arith::cplx<T> &x) {
  using std::isnan;
  return isnan(x.real) || isnan(x.imag);
}

template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST T
norm(const Arith::cplx<T> &x) {
  // using std::min;
  // return min(T(0), real(x * conj(x)));
  return x.real * x.real + x.imag * x.imag;
}

// template <typename T> constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE
// CPLX_HOST Arith::cplx<T> pow(Arith::cplx<T> x, int n) {
//   Arith::cplx<T> r = 1;
//   while (n > 0) {
//     if (n % 2 != 0)
//       r *= x;
//     x *= x;
//     n /= 2;
//   }
//   return r;
// }
template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST Arith::cplx<T>
pow(Arith::cplx<T> x, int n) {
  if (n == 0)
    return 1;
  if (n == 1)
    return x;
  const auto y = pow(x * x, n / 2);
  if (n % 2 == 0)
    return y;
  return x * y;
}

template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST Arith::cplx<T>
pow2(const Arith::cplx<T> &x) {
  return x * x;
}

template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE CPLX_HOST T
real(const Arith::cplx<T> &x) {
  return x.real;
}

// template <typename T> constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE
// CPLX_HOST Arith::cplx<T> sin(const Arith::cplx<T> &x) {
//   using std::cos, std::sin;
//   return {sin(x.real), cos(x.real) * x.imag};
// }

// template <typename T> constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CPLX_DEVICE
// CPLX_HOST Arith::cplx<T> sqrt(const Arith::cplx<T> &x)
// {
//   using std::sqrt;
//   const T r = sqrt(x.real);
//   return {r, x.imag / (2 * r)};
// }

} // namespace Arith

#endif // #ifndef CPLX_HXX
