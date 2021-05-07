#ifndef DIV_HXX
#define DIV_HXX

#include "defs.hxx"

namespace Arith {

template <typename T> constexpr ARITH_INLINE T div_floor(T x, T y) {
  // C++ division truncates; we want to round towards -infinity instead
  x = y < 0 ? -x : x;
  y = y < 0 ? -y : y;
  return (x >= 0 ? x : x - y + 1) / y;
}

template <typename T> constexpr ARITH_INLINE T mod_floor(T x, T y) {
  return x - div_floor(x, y) * y;
}

} // namespace Arith

#endif // #ifndef DIV_HXX
