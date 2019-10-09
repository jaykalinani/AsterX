#ifndef LOOP_HXX
#define LOOP_HXX

#include <cctk.h>
#undef copysign
#undef fpclassify
#undef isfinite
#undef isinf
#undef isnan
#undef isnormal
#undef signbit

#include <array>
#include <iostream>
#include <string>

namespace Loop {
using namespace std;

constexpr int dim = 3;

////////////////////////////////////////////////////////////////////////////////

template <typename T, int D> struct vect {
  array<T, D> elts;

  vect() {
    for (int d = 0; d < D; ++d)
      elts[d] = 0;
  }

  vect(const array<T, D> &arr) : elts(arr) {}

  static vect unit(int dir) {
    vect r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = d == dir;
    return r;
  }

  const T &operator[](int d) const { return elts[d]; }
  T &operator[](int d) { return elts[d]; }

  friend vect operator+(const vect &x) {
    vect r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = +x.elts[d];
    return r;
  }
  friend vect operator-(const vect &x) {
    vect r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = -x.elts[d];
    return r;
  }

  friend vect operator+(const vect &x, const vect &y) {
    vect r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = x.elts[d] + y.elts[d];
    return r;
  }
  friend vect operator-(const vect &x, const vect &y) {
    vect r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = x.elts[d] - y.elts[d];
    return r;
  }

  friend vect operator*(T a, const vect &x) {
    vect r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = a * x.elts[d];
    return r;
  }
  friend vect operator*(const vect &x, T a) {
    vect r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = x.elts[d] * a;
    return r;
  }

  friend ostream &operator<<(ostream &os, vect &x) {
    os << "[";
    for (int d = 0; d < D; ++d) {
      if (d > 0)
        os << ",";
      os << x.elts[d];
    }
    os << "]";
    return os;
  }
};

////////////////////////////////////////////////////////////////////////////////

template <typename T, int CI, int CJ, int CK> struct GF3D {
  static_assert(CI == 0 || CI == 1, "");
  static_assert(CJ == 0 || CJ == 1, "");
  static_assert(CK == 0 || CK == 1, "");
  typedef T value_type;
  T *restrict ptr;
  static constexpr int di = 1;
  int dj, dk;
  int ni, nj, nk;
  static constexpr array<int, dim> indextype() { return {CI, CJ, CK}; }
  inline GF3D(const cGH *restrict cctkGH, T *restrict ptr)
      : ptr(ptr), dj(di * (cctkGH->cctk_ash[0] + 1 - CI)),
        dk(dj * (cctkGH->cctk_ash[1] + 1 - CJ)),
        ni(cctkGH->cctk_lsh[0] + 1 - CI), nj(cctkGH->cctk_lsh[1] + 1 - CJ),
        nk(cctkGH->cctk_lsh[2] + 1 - CK) {}
  inline int offset(int i, int j, int k) const {
    // These index checks prevent vectorization. We thus only enable
    // them in debug mode.
#ifdef CCTK_DEBUG
    assert(i >= 0 && i < ni);
    assert(j >= 0 && j < nj);
    assert(k >= 0 && k < nk);
#endif
    return i * di + j * dj + k * dk;
  }
  inline T &restrict operator()(int i, int j, int k) const {
    return ptr[offset(i, j, k)];
  }
  inline T &restrict operator()(const vect<int, dim> &I) const {
    return ptr[offset(I[0], I[1], I[2])];
  }
};

////////////////////////////////////////////////////////////////////////////////

struct PointDesc {
  int i, j, k;
  CCTK_REAL x, y, z;
  int idx;
  static constexpr int di = 1;
  int dj, dk;
  vect<int, dim> I;
  vect<int, dim> DI(int d) const { return vect<int, dim>::unit(d); }
};

struct GridDescBase {
  array<int, dim> gsh;
  array<int, dim> lbnd, ubnd;
  array<int, dim> lsh;
  array<int, dim> ash;
  array<int, 2 * dim> bbox;
  array<int, dim> nghostzones;
  array<int, dim> tmin, tmax;

  array<CCTK_REAL, dim> x0;
  array<CCTK_REAL, dim> dx;

  template <typename T, size_t N>
  static void output(ostream &os, const string &str, const array<T, N> &arr) {
    os << str << ":[";
    for (size_t n = 0; n < N; ++n) {
      if (n > 0)
        os << ",";
      os << arr[n];
    }
    os << "]";
  }
  friend ostream &operator<<(ostream &os, const GridDescBase &grid) {
    os << "GridDescBase{";
    output(os, "gsh", grid.gsh);
    output(os, ",lbnd", grid.lbnd);
    output(os, ",ubnd", grid.ubnd);
    output(os, ",lsh", grid.lsh);
    output(os, ",bbox", grid.bbox);
    output(os, ",nghostzones", grid.nghostzones);
    output(os, ",tmin", grid.tmin);
    output(os, ",tmax", grid.tmax);
    os << "}";
    return os;
  }

protected:
  GridDescBase();

public:
  GridDescBase(const cGH *cctkGH);

  // Loop over a given box
  template <int CI, int CJ, int CK, typename F>
  void loop_box(const F &f, const array<int, dim> &restrict imin,
                const array<int, dim> &restrict imax) const {
    static_assert(CI == 0 || CI == 1, "");
    static_assert(CJ == 0 || CJ == 1, "");
    static_assert(CK == 0 || CK == 1, "");

    for (int d = 0; d < dim; ++d)
      if (imin[d] >= imax[d])
        return;

    constexpr int di = 1;
    const int dj = di * (ash[0] + 1 - CI);
    const int dk = dj * (ash[1] + 1 - CJ);

    for (int k = imin[2]; k < imax[2]; ++k) {
      for (int j = imin[1]; j < imax[1]; ++j) {
#pragma omp simd
        for (int i = imin[0]; i < imax[0]; ++i) {
          CCTK_REAL x = x0[0] + (lbnd[0] + i + CCTK_REAL(CI - 1) / 2) * dx[0];
          CCTK_REAL y = x0[1] + (lbnd[1] + j + CCTK_REAL(CJ - 1) / 2) * dx[1];
          CCTK_REAL z = x0[2] + (lbnd[2] + k + CCTK_REAL(CK - 1) / 2) * dx[2];
          int idx = i * di + j * dj + k * dk;
          vect<int, dim> I{{i, j, k}};
          const PointDesc p{i, j, k, x, y, z, idx, dj, dk, I};
          f(p);
        }
      }
    }
  }

  // Loop over all points
  template <int CI, int CJ, int CK, typename F>
  void loop_all(const F &f) const {
    constexpr array<int, dim> offset{!CI, !CJ, !CK};
    array<int, dim> imin, imax;
    for (int d = 0; d < dim; ++d) {
      imin[d] = max(tmin[d], 0);
      imax[d] = min(tmax[d] + (tmax[d] >= lsh[d] ? offset[d] : 0),
                    lsh[d] + offset[d]);
    }

    loop_box<CI, CJ, CK>(f, imin, imax);
  }

  // Loop over all interior points
  template <int CI, int CJ, int CK, typename F>
  void loop_int(const F &f) const {
    constexpr array<int, dim> offset{!CI, !CJ, !CK};
    array<int, dim> imin, imax;
    for (int d = 0; d < dim; ++d) {
      imin[d] = max(tmin[d], nghostzones[d]);
      imax[d] = min(tmax[d] + (tmax[d] >= lsh[d] ? offset[d] : 0),
                    lsh[d] + offset[d] - nghostzones[d]);
    }

    loop_box<CI, CJ, CK>(f, imin, imax);
  }

  // Loop over all outer boundary points. This excludes ghost faces, but
  // includes ghost edges/corners on non-ghost faces.
  template <int CI, int CJ, int CK, typename F>
  void loop_bnd(const F &f) const {
    constexpr array<int, dim> offset{!CI, !CJ, !CK};

    for (int dir = 0; dir < dim; ++dir) {
      for (int face = 0; face < 2; ++face) {
        if (bbox[2 * dir + face]) {

          array<int, dim> imin, imax;
          for (int d = 0; d < dim; ++d) {
            // by default, include interior and outer boundaries and ghosts
            imin[d] = 0;
            imax[d] = lsh[d] + offset[d];

            // avoid covering edges and corners multiple times
            if (d < dir) {
              if (bbox[2 * d])
                imin[d] = nghostzones[d]; // only interior
              if (bbox[2 * d + 1])
                imax[d] = lsh[d] + offset[d] - nghostzones[d]; // only interior
            }
          }
          // only one face on outer boundary
          if (face == 0)
            imax[dir] = nghostzones[dir];
          else
            imin[dir] = lsh[dir] + offset[dir] - nghostzones[dir];

          for (int d = 0; d < dim; ++d) {
            imin[d] = max(tmin[d], imin[d]);
            imax[d] =
                min(tmax[d] + (tmax[d] >= lsh[d] ? offset[d] : 0), imax[d]);
          }

          loop_box<CI, CJ, CK>(f, imin, imax);
        }
      }
    }
  }

  template <typename F>
  void loop_all(const array<int, dim> &indextype, const F &f) const {
    // typedef void (GridDescBase::*funptr)(const F &f) const;
    // constexpr array<funptr, 8> funptrs{
    //     &GridDescBase::loop_all<0, 0, 0, F>,
    //     &GridDescBase::loop_all<1, 0, 0, F>,
    //     &GridDescBase::loop_all<0, 1, 0, F>,
    //     &GridDescBase::loop_all<1, 1, 0, F>,
    //     &GridDescBase::loop_all<0, 0, 1, F>,
    //     &GridDescBase::loop_all<1, 0, 1, F>,
    //     &GridDescBase::loop_all<0, 1, 1, F>,
    //     &GridDescBase::loop_all<1, 1, 1, F>,
    // };
    // return (this->*funptrs[indextype[0] + 2 * indextype[1] + 4 *
    // indextype[2]])(
    //     f);

    switch (indextype[0] + 2 * indextype[1] + 4 * indextype[2]) {
    case 0b000:
      return loop_all<0, 0, 0>(f);
    case 0b001:
      return loop_all<1, 0, 0>(f);
    case 0b010:
      return loop_all<0, 1, 0>(f);
    case 0b011:
      return loop_all<1, 1, 0>(f);
    case 0b100:
      return loop_all<0, 0, 1>(f);
    case 0b101:
      return loop_all<1, 0, 1>(f);
    case 0b110:
      return loop_all<0, 1, 1>(f);
    case 0b111:
      return loop_all<1, 1, 1>(f);
    default:
      assert(0);
    }
  }

  template <typename F>
  void loop_int(const array<int, dim> &indextype, const F &f) const {
    switch (indextype[0] + 2 * indextype[1] + 4 * indextype[2]) {
    case 0b000:
      return loop_int<0, 0, 0>(f);
    case 0b001:
      return loop_int<1, 0, 0>(f);
    case 0b010:
      return loop_int<0, 1, 0>(f);
    case 0b011:
      return loop_int<1, 1, 0>(f);
    case 0b100:
      return loop_int<0, 0, 1>(f);
    case 0b101:
      return loop_int<1, 0, 1>(f);
    case 0b110:
      return loop_int<0, 1, 1>(f);
    case 0b111:
      return loop_int<1, 1, 1>(f);
    default:
      assert(0);
    }
  }

  template <typename F>
  void loop_bnd(const array<int, dim> &indextype, const F &f) const {
    switch (indextype[0] + 2 * indextype[1] + 4 * indextype[2]) {
    case 0b000:
      return loop_bnd<0, 0, 0>(f);
    case 0b001:
      return loop_bnd<1, 0, 0>(f);
    case 0b010:
      return loop_bnd<0, 1, 0>(f);
    case 0b011:
      return loop_bnd<1, 1, 0>(f);
    case 0b100:
      return loop_bnd<0, 0, 1>(f);
    case 0b101:
      return loop_bnd<1, 0, 1>(f);
    case 0b110:
      return loop_bnd<0, 1, 1>(f);
    case 0b111:
      return loop_bnd<1, 1, 1>(f);
    default:
      assert(0);
    }
  }
};

template <int CI, int CJ, int CK, typename F>
void loop_all(const cGH *cctkGH, const F &f) {
  GridDescBase(cctkGH).loop_all<CI, CJ, CK>(f);
}

template <int CI, int CJ, int CK, typename F>
void loop_int(const cGH *cctkGH, const F &f) {
  GridDescBase(cctkGH).loop_int<CI, CJ, CK>(f);
}

template <int CI, int CJ, int CK, typename F>
void loop_bnd(const cGH *cctkGH, const F &f) {
  GridDescBase(cctkGH).loop_bnd<CI, CJ, CK>(f);
}

template <typename F>
void loop_all(const cGH *cctkGH, const array<int, dim> &indextype, const F &f) {
  GridDescBase(cctkGH).loop_all(indextype, f);
}

template <typename F>
void loop_int(const cGH *cctkGH, const array<int, dim> &indextype, const F &f) {
  GridDescBase(cctkGH).loop_int(indextype, f);
}

template <typename F>
void loop_bnd(const cGH *cctkGH, const array<int, dim> &indextype, const F &f) {
  GridDescBase(cctkGH).loop_bnd(indextype, f);
}

} // namespace Loop

#endif // #ifndef LOOP_HXX
