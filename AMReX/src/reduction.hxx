#ifndef REDUCTION_HXX
#define REDUCTION_HXX

#include <cctk.h>

#include <cmath>

namespace AMReX {
using namespace std;

template <typename T> struct reduction {
  T min, max, sum, sum2;
  T vol, maxabs, sumabs, sum2abs;
  reduction();
  reduction(const T &V, const T &x);
  reduction(const reduction &x, const reduction &y);
  reduction operator+(const reduction &y) const { return reduction(*this, y); }
  reduction &operator+=(const reduction &y) { return *this = *this + y; }

  T avg() const { return sum / vol; }
  T sdv() const { return sqrt(fmax(T(0.0), vol * sum2 - sum * sum)); }
  T norm0() const { return vol; }
  T norm1() const { return sumabs / vol; }
  T norm2() const { return sqrt(sum2abs / vol); }
  T norm_inf() const { return maxabs; }
};

template <typename T>
reduction<T>::reduction()
    : min(1.0 / 0.0), max(-1.0 / 0.0), sum(0.0), sum2(0.0), vol(0.0),
      maxabs(0.0), sumabs(0.0), sum2abs(0.0) {}

template <typename T>
reduction<T>::reduction(const T &V, const T &x)
    : min(x), max(x), sum(V * x), sum2((V * x) * (V * x)), vol(V),
      maxabs(fabs(x)), sumabs(fabs(V * x)), sum2abs(fabs(V * x) * fabs(V * x)) {
}

template <typename T>
reduction<T>::reduction(const reduction &x, const reduction &y)
    : min(fmin(x.min, y.min)), max(fmax(x.max, y.max)), sum(x.sum + y.sum),
      sum2(x.sum2 + y.sum2), vol(x.vol + y.vol),
      maxabs(fmax(x.maxabs, y.maxabs)), sumabs(x.sumabs + y.sumabs),
      sum2abs(x.sum2abs + y.sum2abs) {}

typedef reduction<CCTK_REAL> reduction_CCTK_REAL;
#pragma omp declare reduction(reduction:reduction_CCTK_REAL : omp_out += omp_in)

reduction<CCTK_REAL> reduce(int gi, int vi, int tl);

} // namespace AMReX

#endif // #ifndef REDUCTION_HXX
