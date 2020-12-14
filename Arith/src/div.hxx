#ifndef DIV_HXX
#define DIV_HXX

#include <cassert>
#include <cmath>

#include <cctk.h>

#ifdef CCTK_DEBUG
#define ARITH_INLINE
#else
#define ARITH_INLINE CCTK_ATTRIBUTE_ALWAYS_INLINE
#endif

namespace Arith {
using namespace std;

template <typename T> constexpr T div_floor(T x, T y) ARITH_INLINE;
template <typename T> constexpr T div_floor(T x, T y) {
  // C++ division truncates; we want to round towards -infinity instead
  x = y < 0 ? -x : x;
  y = y < 0 ? -y : y;
  return (x >= 0 ? x : x - y + 1) / y;
}

template <typename T> constexpr T mod_floor(T x, T y) ARITH_INLINE;
template <typename T> constexpr T mod_floor(T x, T y) {
  return x - div_floor(x, y) * y;
}

} // namespace Arith

#undef ARITH_INLINE

#endif // #ifndef DIV_HXX
