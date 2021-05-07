#ifndef DEFS_HH
#define DEFS_HH

#include <fixmath.hxx>
#include <simd.hxx>

#include <array>
#include <cmath>
#include <tuple>

namespace Hydro {
using namespace Arith;
using namespace std;

constexpr int dim = 3;

template <typename T> inline T fmax3(T x0, T x1, T x2) {
  T x01 = fmax(x0, x1);
  T x012 = fmax(x2, x01);
  return x012;
}

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE T fmax5(T x0, T x1, T x2, T x3, T x4) {
  T x01 = fmax(x0, x1);
  T x23 = fmax(x2, x3);
  T x014 = fmax(x4, x01);
  T x01234 = fmax(x23, x014);
  return x01234;
}

typedef simd<CCTK_REAL> CCTK_REALVEC;
typedef simdl<CCTK_REAL> CCTK_BOOLVEC;
constexpr int vsize = tuple_size_v<CCTK_REALVEC>;

} // namespace Hydro

#endif // #ifndef DEFS_HH
