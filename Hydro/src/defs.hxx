#ifndef DEFS_HH
#define DEFS_HH

#include <fixmath.hxx>
#include <vectors.h>

#include <array>
#include <cmath>

namespace Hydro {
using namespace std;

constexpr int dim = 3;

template <typename T> inline CCTK_ATTRIBUTE_ALWAYS_INLINE T pow2(T x) {
  return x * x;
}

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

typedef vectype<CCTK_REAL> CCTK_REALVEC;
constexpr int vsize = vecprops<CCTK_REAL>::size();

inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_REALVEC
vloadu(const CCTK_REAL &restrict x) {
  return CCTK_REALVEC::loadu(x);
}

inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_REALVEC viota() {
  array<CCTK_REAL, vsize> iota_;
  for (int i = 0; i < vsize; ++i)
    iota_[i] = i;
  return CCTK_REALVEC::loadu(iota_[0]);
}

inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_REALVEC vsin(CCTK_REALVEC x) {
  return ksin(x);
}

} // namespace Hydro

#endif // #ifndef DEFS_HH
