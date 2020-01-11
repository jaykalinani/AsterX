#ifndef REDUCTION_HXX
#define REDUCTION_HXX

#include <cctk.h>
#undef copysign
#undef fpclassify
#undef isfinite
#undef isinf
#undef isnan
#undef isnormal
#undef signbit

#include <mpi.h>

#include <cmath>
#include <ostream>

namespace CarpetX {
using namespace std;

template <typename T> struct mpi_datatype;
template <> struct mpi_datatype<float> {
  static constexpr MPI_Datatype value = MPI_FLOAT;
};
template <> struct mpi_datatype<double> {
  static constexpr MPI_Datatype value = MPI_DOUBLE;
};
template <> struct mpi_datatype<long double> {
  static constexpr MPI_Datatype value = MPI_LONG_DOUBLE;
};

template <typename T> T pow21(T x) { return x * x; }

template <typename T> T fmax1(T x, T y) {
  if (CCTK_isnan(CCTK_REAL(x)))
    return x;
  if (CCTK_isnan(CCTK_REAL(y)))
    return y;
  return fmax(x, y);
}
template <typename T> T fmin1(T x, T y) {
  if (CCTK_isnan(CCTK_REAL(x)))
    return x;
  if (CCTK_isnan(CCTK_REAL(y)))
    return y;
  return fmin(x, y);
}

template <typename T> struct reduction {
  // TODO: contains_inf, contains_nan?
  // TODO: minloc, maxloc
  T min, max, sum, sum2;
  T vol, maxabs, sumabs, sum2abs;
  reduction();
  reduction(const T &V, const T &x);
  reduction(const reduction &x, const reduction &y);
  reduction operator+(const reduction &y) const { return reduction(*this, y); }
  reduction &operator+=(const reduction &y) { return *this = *this + y; }

  T avg() const { return sum / vol; }
  T sdv() const { return sqrt(fmax1(T(0), vol * sum2 - pow21(sum))); }
  T norm0() const { return vol; }
  T norm1() const { return sumabs / vol; }
  T norm2() const { return sqrt(sum2abs / vol); }
  T norm_inf() const { return maxabs; }

  friend ostream &operator<<(ostream &os, const reduction &red) {
    return os << "reduction{min:" << red.min << ",max:" << red.max
              << ",sum:" << red.sum << ",sum2:" << red.sum2
              << ",vol:" << red.vol << ",maxabs:" << red.maxabs
              << ",sumabs:" << red.sumabs << ",sum2abs:" << red.sum2abs << "}";
  }
};

template <typename T>
reduction<T>::reduction()
    : min(1.0 / 0.0), max(-1.0 / 0.0), sum(0.0), sum2(0.0), vol(0.0),
      maxabs(0.0), sumabs(0.0), sum2abs(0.0) {}

template <typename T>
reduction<T>::reduction(const T &V, const T &x)
    : min(x), max(x), sum(V * x), sum2(V * pow21(x)), vol(V), maxabs(fabs(x)),
      sumabs(V * fabs(x)), sum2abs(V * pow21(fabs(x))) {}

template <typename T>
reduction<T>::reduction(const reduction &x, const reduction &y)
    : min(fmin1(x.min, y.min)), max(fmax1(x.max, y.max)), sum(x.sum + y.sum),
      sum2(x.sum2 + y.sum2), vol(x.vol + y.vol),
      maxabs(fmax1(x.maxabs, y.maxabs)), sumabs(x.sumabs + y.sumabs),
      sum2abs(x.sum2abs + y.sum2abs) {}

typedef reduction<CCTK_REAL> reduction_CCTK_REAL;
#pragma omp declare reduction(reduction:reduction_CCTK_REAL : omp_out += omp_in)

MPI_Datatype reduction_mpi_datatype_CCTK_REAL();
MPI_Op reduction_mpi_op();

reduction<CCTK_REAL> reduce(int gi, int vi, int tl);

} // namespace CarpetX

#endif // #ifndef REDUCTION_HXX
