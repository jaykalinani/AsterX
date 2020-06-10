#ifndef SUM_HXX
#define SUM_HXX

#include <fixmath.hxx> // include this before <cctk.h>
#include <cctk.h>

#include <functional>

#ifdef CCTK_DEBUG
#define ARITH_INLINE
#else
#define ARITH_INLINE CCTK_ATTRIBUTE_ALWAYS_INLINE
#endif

namespace Arith {
using namespace std;

template <int D, typename F>
constexpr ARITH_INLINE remove_cv_t<remove_reference_t<result_of_t<F(int)> > >
sum(F f) {
  using R = remove_cv_t<remove_reference_t<result_of_t<F(int)> > >;
  R s{0};
  for (int x = 0; x < D; ++x)
    s += f(x);
  return s;
}

template <int D, typename F>
constexpr
    ARITH_INLINE remove_cv_t<remove_reference_t<result_of_t<F(int, int)> > >
    sum(F f) {
  using R = remove_cv_t<remove_reference_t<result_of_t<F(int, int)> > >;
  R s{0};
  for (int x = 0; x < D; ++x)
    for (int y = 0; y < D; ++y)
      s += f(x, y);
  return s;
}

template <int D, typename F>
constexpr ARITH_INLINE
    remove_cv_t<remove_reference_t<result_of_t<F(int, int, int)> > >
    sum(F f) {
  using R = remove_cv_t<remove_reference_t<result_of_t<F(int, int, int)> > >;
  R s{0};
  for (int x = 0; x < D; ++x)
    for (int y = 0; y < D; ++y)
      for (int z = 0; z < D; ++z)
        s += f(x, y, z);
  return s;
}

} // namespace Arith

#undef ARITH_INLINE

#endif
