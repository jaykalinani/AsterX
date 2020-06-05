#ifndef FIXMATH_HXX
#define FIXMATH_HXX

// Include this file BEFORE including <cctk.h>, best before including any other
// include files

#include <cmath>

namespace std {
template <typename T> constexpr T copysign1(T x, T y) {
  return std::copysign(x, y);
}
template <typename T> constexpr int fpclassify1(T x) {
  return std::fpclassify(x);
}
template <typename T> constexpr bool isfinite1(T x) { return std::isfinite(x); }
template <typename T> constexpr bool isinf1(T x) { return std::isinf(x); }
template <typename T> constexpr bool isnan1(T x) { return std::isnan(x); }
template <typename T> constexpr bool isnormal1(T x) { return std::isnormal(x); }
template <typename T> constexpr bool signbit1(T x) { return std::signbit(x); }
} // namespace std

#include <cctk.h>

#undef copysign
#undef fpclassify
#undef isfinite
#undef isinf
#undef isnan
#undef isnormal
#undef signbit

#endif // #ifndef FIXMATH_HXX
