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

template <typename T> constexpr T div(T x, T y) ARITH_INLINE {
  // C++ division truncates; we want to round towards -infinity instead
  if (y < 0) {
    x = -x;
    y = -y;
  }
  return x >= 0 ? x / y : (x - y + 1) / y;
}

template <typename T> constexpr T mod(const T &x, const T &y) ARITH_INLINE {
  return x - div(x, y) * y;
}

} // namespace Arith

#undef ARITH_INLINE

#endif // #ifndef DIV_HXX
