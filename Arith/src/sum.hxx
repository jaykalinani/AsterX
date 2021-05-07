#ifndef SUM_HXX
#define SUM_HXX

#include "defs.hxx"
#include "vect.hxx"

#include <functional>

namespace Arith {
using namespace std;

template <int D, typename F>
constexpr ARITH_INLINE remove_cv_t<remove_reference_t<result_of_t<F(int)> > >
sum(F f) {
  typedef remove_cv_t<remove_reference_t<result_of_t<F(int)> > > R;
  R s = zero<R>();
  for (int x = 0; x < D; ++x)
    s += f(x);
  return s;
}

template <int D, typename F>
constexpr
    ARITH_INLINE remove_cv_t<remove_reference_t<result_of_t<F(int, int)> > >
    sum(F f) {
  typedef remove_cv_t<remove_reference_t<result_of_t<F(int, int)> > > R;
  R s = zero<R>();
  for (int x = 0; x < D; ++x)
    for (int y = 0; y < D; ++y)
      s += f(x, y);
  return s;
}

template <int D, typename F>
constexpr ARITH_INLINE
    remove_cv_t<remove_reference_t<result_of_t<F(int, int, int)> > >
    sum(F f) {
  typedef remove_cv_t<remove_reference_t<result_of_t<F(int, int, int)> > > R;
  R s = zero<R>();
  for (int x = 0; x < D; ++x)
    for (int y = 0; y < D; ++y)
      for (int z = 0; z < D; ++z)
        s += f(x, y, z);
  return s;
}

template <int D, typename F>
constexpr ARITH_INLINE
    remove_cv_t<remove_reference_t<result_of_t<F(int, int, int, int)> > >
    sum(F f) {
  typedef remove_cv_t<remove_reference_t<result_of_t<F(int, int, int, int)> > >
      R;
  R s = zero<R>();
  for (int x = 0; x < D; ++x)
    for (int y = 0; y < D; ++y)
      for (int z = 0; z < D; ++z)
        for (int w = 0; w < D; ++w)
          s += f(x, y, z, w);
  return s;
}

template <int D, typename F>
constexpr
    ARITH_INLINE remove_cv_t<remove_reference_t<result_of_t<F(int, int)> > >
    sum_symm(F f) {
  typedef remove_cv_t<remove_reference_t<result_of_t<F(int, int)> > > R;
  R s = zero<R>();
  for (int x = 0; x < D; ++x)
    for (int y = x; y < D; ++y)
      s += (x == y ? 1 : 2) * f(x, y);
  return s;
}

template <int D, typename F>
constexpr ARITH_INLINE
    remove_cv_t<remove_reference_t<result_of_t<F(int, int, int)> > >
    sum_symm(F f) {
  typedef remove_cv_t<remove_reference_t<result_of_t<F(int, int, int)> > > R;
  R s = zero<R>();
  for (int x = 0; x < D; ++x)
    for (int y = x; y < D; ++y)
      for (int z = y; z < D; ++z)
        s += (x == y && x == z             ? 1
              : x == y || x == z || y == z ? 2
                                           : 3) *
             f(x, y, z);
  return s;
}

} // namespace Arith

#endif
