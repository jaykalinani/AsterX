#ifndef LOOP_HXX
#define LOOP_HXX

#include <vect.hxx>

#include <cctk.h>

#include <algorithm>
#include <array>
#include <functional>
#include <initializer_list>
#include <limits>
#include <ostream>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

namespace Loop {
using namespace std;

using Arith::vect;

constexpr int dim = 3;

enum class where_t { everywhere, interior, boundary, ghosts_inclusive, ghosts };

struct PointDesc {
  int imin, imax;
  int i, j, k;
  CCTK_REAL x, y, z;
  CCTK_REAL dx, dy, dz;
  int idx;
  static constexpr int di = 1;
  int dj, dk, np;
  vect<int, dim> I;
  vect<int, dim> NI; // outward boundary normal, or zero
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect<int, dim> DI(int d) const {
    return vect<int, dim>::unit(d);
  }
  vect<CCTK_REAL, dim> X;
  vect<CCTK_REAL, dim> DX;

  friend ostream &operator<<(ostream &os, const PointDesc &p) {
    os << "PointDesc{"
       << "imin,imax:{" << p.imin << "," << p.imax << "}, "
       << "ijk:"
       << "{" << p.i << "," << p.j << "," << p.k << "}, "
       << "xyz:"
       << "{" << p.x << "," << p.y << "," << p.z << "}, "
       << "dxyz:"
       << "{" << p.dx << "," << p.dy << "," << p.dz << "}, "
       << "idx:" << p.idx << ", "
       << "dijk:"
       << "{" << p.di << "," << p.dj << "," << p.dk << "},"
       << "np:" << p.np << ","
       << "nijk:"
       << "{" << p.NI[0] << "," << p.NI[1] << "," << p.NI[2] << "}"
       << "}";
    return os;
  }
};

struct GridDescBase {
  array<int, dim> gsh;
  array<int, dim> lbnd, ubnd;
  array<int, dim> lsh;
  array<int, dim> ash;
  array<int, 2 * dim> bbox;
  array<int, dim> nghostzones;
  array<int, dim> tmin, tmax;

  // for current level
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

  template <int CI, int CJ, int CK>
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE PointDesc
  point_desc(const vect<int, dim> &restrict NI, const int imin, const int imax,
             const int i, const int j, const int k) const {
    constexpr int di = 1;
    const int dj = di * (ash[0] + !CI);
    const int dk = dj * (ash[1] + !CJ);
    const int np = dk * (ash[2] + !CK);

    const CCTK_REAL x = x0[0] + (lbnd[0] + i - CCTK_REAL(!CI) / 2) * dx[0];
    const CCTK_REAL y = x0[1] + (lbnd[1] + j - CCTK_REAL(!CJ) / 2) * dx[1];
    const CCTK_REAL z = x0[2] + (lbnd[2] + k - CCTK_REAL(!CK) / 2) * dx[2];
    const int idx = i * di + j * dj + k * dk;
    return PointDesc{imin, imax,      i,     j,         k,   x,  y,
                     z,    dx[0],     dx[1], dx[2],     idx, dj, dk,
                     np,   {i, j, k}, NI,    {x, y, z}, dx};
  }

  template <int CI, int CJ, int CK>
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE PointDesc
  point_desc(const PointDesc &p) const {
    return point_desc<CI, CJ, CK>(p.NI, p.imin, p.imax, p.i, p.j, p.k);
  }

  // Loop over a given box
  template <int CI, int CJ, int CK, typename F>
  CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_box(const F &f, const array<int, dim> &restrict imin,
           const array<int, dim> &restrict imax,
           const array<int, dim> &restrict inormal) const {
    static_assert(CI == 0 || CI == 1, "");
    static_assert(CJ == 0 || CJ == 1, "");
    static_assert(CK == 0 || CK == 1, "");

    for (int d = 0; d < dim; ++d)
      if (imin[d] >= imax[d])
        return;

    const auto kernel{[&](const int i, const int j,
                          const int k) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      f(point_desc<CI, CJ, CK>(inormal, imin[0], imax[0], i, j, k));
    }};

    array<bool, dim> bforward;
    for (int d = 0; d < dim; ++d)
      bforward[d] = inormal[d] >= 0;
    bool all_forward = true;
    for (int d = 0; d < dim; ++d)
      all_forward &= bforward[d];

    if (all_forward) {

      for (int k = imin[2]; k < imax[2]; ++k) {
        for (int j = imin[1]; j < imax[1]; ++j) {
#pragma omp simd
          for (int i = imin[0]; i < imax[0]; ++i) {
            kernel(i, j, k);
          }
        }
      }

    } else {
      // At least one direction is reversed; loop from the inside out

      for (int k0 = 0; k0 < imax[2] - imin[2]; ++k0) {
        for (int j0 = 0; j0 < imax[1] - imin[1]; ++j0) {
#pragma omp simd
          for (int i0 = 0; i0 < imax[0] - imin[0]; ++i0) {
            const int i = bforward[0] ? imin[0] + i0 : imax[0] - 1 - i0;
            const int j = bforward[1] ? imin[1] + j0 : imax[1] - 1 - j0;
            const int k = bforward[2] ? imin[2] + k0 : imax[2] - 1 - k0;
            kernel(i, j, k);
          }
        }
      }
    }
  }

  // Box including all points
  template <int CI, int CJ, int CK>
  void box_all(const array<int, dim> &group_nghostzones,
               vect<int, dim> &restrict imin,
               vect<int, dim> &restrict imax) const {
    const array<int, dim> offset{!CI, !CJ, !CK};
    for (int d = 0; d < dim; ++d) {
      int ghost_offset = nghostzones[d] - group_nghostzones[d];
      imin[d] = std::max(tmin[d], ghost_offset);
      imax[d] = std::min(tmax[d] + (tmax[d] >= lsh[d] ? offset[d] : 0),
                         lsh[d] + offset[d] - ghost_offset);
    }
  }

  // Box including all interior points
  template <int CI, int CJ, int CK>
  void box_int(const array<int, dim> &group_nghostzones,
               vect<int, dim> &restrict imin,
               vect<int, dim> &restrict imax) const {
    const array<int, dim> offset{!CI, !CJ, !CK};
    for (int d = 0; d < dim; ++d) {
      imin[d] = std::max(tmin[d], nghostzones[d]);
      imax[d] = std::min(tmax[d] + (tmax[d] >= lsh[d] ? offset[d] : 0),
                         lsh[d] + offset[d] - nghostzones[d]);
    }
  }

  // Loop over all points
  template <int CI, int CJ, int CK, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_all(const array<int, dim> &group_nghostzones, const F &f) const {
    const array<int, dim> offset{!CI, !CJ, !CK};
    array<int, dim> imin, imax;
    for (int d = 0; d < dim; ++d) {
      int ghost_offset = nghostzones[d] - group_nghostzones[d];
      imin[d] = std::max(tmin[d], ghost_offset);
      imax[d] = std::min(tmax[d] + (tmax[d] >= lsh[d] ? offset[d] : 0),
                         lsh[d] + offset[d] - ghost_offset);
    }
    const array<int, dim> inormal{0, 0, 0};

    loop_box<CI, CJ, CK>(f, imin, imax, inormal);
  }

  // Loop over all interior points
  template <int CI, int CJ, int CK, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_int(const array<int, dim> &group_nghostzones, const F &f) const {
    const array<int, dim> offset{!CI, !CJ, !CK};
    array<int, dim> imin, imax;
    for (int d = 0; d < dim; ++d) {
      imin[d] = std::max(tmin[d], nghostzones[d]);
      imax[d] = std::min(tmax[d] + (tmax[d] >= lsh[d] ? offset[d] : 0),
                         lsh[d] + offset[d] - nghostzones[d]);
    }
    const array<int, dim> inormal{0, 0, 0};

    loop_box<CI, CJ, CK>(f, imin, imax, inormal);
  }

  // Loop over all outer boundary points. This excludes ghost faces, but
  // includes ghost edges/corners on non-ghost faces. Loop over faces first,
  // then edges, then corners.
  template <int CI, int CJ, int CK, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_bnd(const array<int, dim> &group_nghostzones, const F &f) const {
    const array<int, dim> offset{!CI, !CJ, !CK};

    for (int rank = dim - 1; rank >= 0; --rank) {

      for (int nk = -1; nk <= +1; ++nk) {
        for (int nj = -1; nj <= +1; ++nj) {
          for (int ni = -1; ni <= +1; ++ni) {
            if ((ni == 0) + (nj == 0) + (nk == 0) == rank) {

              if ((ni != 0 && bbox[0 + (ni == -1 ? 0 : 1)]) ||
                  (nj != 0 && bbox[2 + (nj == -1 ? 0 : 1)]) ||
                  (nk != 0 && bbox[4 + (nk == -1 ? 0 : 1)])) {

                const array<int, dim> inormal{ni, nj, nk};

                array<int, dim> imin, imax;
                for (int d = 0; d < dim; ++d) {
                  const int ghost_offset =
                      nghostzones[d] - group_nghostzones[d];
                  const int begin_bnd = ghost_offset;
                  const int begin_int = nghostzones[d];
                  const int end_int = lsh[d] + offset[d] - nghostzones[d];
                  const int end_bnd = lsh[d] + offset[d] - ghost_offset;
                  switch (inormal[d]) {
                  case -1: // lower boundary
                    imin[d] = begin_bnd;
                    imax[d] = begin_int;
                    break;
                  case 0: // interior
                    imin[d] = begin_int;
                    imax[d] = end_int;
                    break;
                  case +1: // upper boundary
                    imin[d] = end_int;
                    imax[d] = end_bnd;
                    break;
                  default:
                    assert(0);
                  }

                  imin[d] = std::max(tmin[d], imin[d]);
                  imax[d] = std::min(
                      tmax[d] + (tmax[d] >= lsh[d] ? offset[d] : 0), imax[d]);
                }

                loop_box<CI, CJ, CK>(f, imin, imax, inormal);
              }
            } // if rank
          }
        }
      }

    } // for rank
  }

  // Loop over all outer ghost points. This includes ghost edges/corners on
  // non-ghost faces. Loop over faces first, then edges, then corners.
  template <int CI, int CJ, int CK, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_ghosts_inclusive(const array<int, dim> &group_nghostzones,
                        const F &f) const {
    const array<int, dim> offset{!CI, !CJ, !CK};

    for (int rank = dim - 1; rank >= 0; --rank) {

      for (int nk = -1; nk <= +1; ++nk) {
        for (int nj = -1; nj <= +1; ++nj) {
          for (int ni = -1; ni <= +1; ++ni) {
            if ((ni == 0) + (nj == 0) + (nk == 0) == rank) {

              if ((ni != 0 && !bbox[0 + (ni == -1 ? 0 : 1)]) ||
                  (nj != 0 && !bbox[2 + (nj == -1 ? 0 : 1)]) ||
                  (nk != 0 && !bbox[4 + (nk == -1 ? 0 : 1)])) {

                const array<int, dim> inormal{ni, nj, nk};

                array<int, dim> imin, imax;
                for (int d = 0; d < dim; ++d) {
                  const int ghost_offset =
                      nghostzones[d] - group_nghostzones[d];
                  const int begin_bnd = ghost_offset;
                  const int begin_int = nghostzones[d];
                  const int end_int = lsh[d] + offset[d] - nghostzones[d];
                  const int end_bnd = lsh[d] + offset[d] - ghost_offset;
                  switch (inormal[d]) {
                  case -1: // lower boundary
                    imin[d] = begin_bnd;
                    imax[d] = begin_int;
                    break;
                  case 0: // interior
                    imin[d] = begin_int;
                    imax[d] = end_int;
                    break;
                  case +1: // upper boundary
                    imin[d] = end_int;
                    imax[d] = end_bnd;
                    break;
                  default:
                    assert(0);
                  }

                  imin[d] = std::max(tmin[d], imin[d]);
                  imax[d] = std::min(
                      tmax[d] + (tmax[d] >= lsh[d] ? offset[d] : 0), imax[d]);
                }

                loop_box<CI, CJ, CK>(f, imin, imax, inormal);
              }
            } // if rank
          }
        }
      }
    } // for rank
  }

  // Loop over all outer ghost points. This excludes ghost edges/corners on
  // non-ghost faces. Loop over faces first, then edges, then corners.
  template <int CI, int CJ, int CK, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_ghosts(const array<int, dim> &group_nghostzones, const F &f) const {
    const array<int, dim> offset{!CI, !CJ, !CK};

    for (int rank = dim - 1; rank >= 0; --rank) {

      for (int nk = -1; nk <= +1; ++nk) {
        for (int nj = -1; nj <= +1; ++nj) {
          for (int ni = -1; ni <= +1; ++ni) {
            if ((ni == 0) + (nj == 0) + (nk == 0) == rank) {

              if ((ni == 0 || !bbox[0 + (ni == -1 ? 0 : 1)]) &&
                  (nj == 0 || !bbox[2 + (nj == -1 ? 0 : 1)]) &&
                  (nk == 0 || !bbox[4 + (nk == -1 ? 0 : 1)])) {

                const array<int, dim> inormal{ni, nj, nk};

                array<int, dim> imin, imax;
                for (int d = 0; d < dim; ++d) {
                  const int ghost_offset =
                      nghostzones[d] - group_nghostzones[d];
                  const int begin_bnd = ghost_offset;
                  const int begin_int = nghostzones[d];
                  const int end_int = lsh[d] + offset[d] - nghostzones[d];
                  const int end_bnd = lsh[d] + offset[d] - ghost_offset;
                  switch (inormal[d]) {
                  case -1: // lower boundary
                    imin[d] = begin_bnd;
                    imax[d] = begin_int;
                    break;
                  case 0: // interior
                    imin[d] = begin_int;
                    imax[d] = end_int;
                    break;
                  case +1: // upper boundary
                    imin[d] = end_int;
                    imax[d] = end_bnd;
                    break;
                  default:
                    assert(0);
                  }

                  imin[d] = std::max(tmin[d], imin[d]);
                  imax[d] = std::min(
                      tmax[d] + (tmax[d] >= lsh[d] ? offset[d] : 0), imax[d]);
                }

                loop_box<CI, CJ, CK>(f, imin, imax, inormal);
              }
            } // if rank
          }
        }
      }
    } // for rank
  }

  template <int CI, int CJ, int CK, typename F>
  void loop(where_t where, const array<int, dim> &group_nghostzones,
            const F &f) const {
    switch (where) {
    case where_t::everywhere:
      return loop_all<CI, CJ, CK>(group_nghostzones, f);
    case where_t::interior:
      return loop_int<CI, CJ, CK>(group_nghostzones, f);
    case where_t::boundary:
      return loop_bnd<CI, CJ, CK>(group_nghostzones, f);
    case where_t::ghosts_inclusive:
      return loop_ghosts_inclusive<CI, CJ, CK>(group_nghostzones, f);
    case where_t::ghosts:
      return loop_ghosts<CI, CJ, CK>(group_nghostzones, f);
    default:
      assert(0);
    }
  }

  template <int CI, int CJ, int CK, typename F>
  void loop(where_t where, const F &f) const {
    loop<CI, CJ, CK>(where, nghostzones, f);
  }

  template <typename F>
  void loop_idx(where_t where, const array<int, dim> &indextype,
                const array<int, dim> &group_nghostzones, const F &f) const {
    switch (indextype[0] + 2 * indextype[1] + 4 * indextype[2]) {
    case 0b000:
      return loop<0, 0, 0>(where, group_nghostzones, f);
    case 0b001:
      return loop<1, 0, 0>(where, group_nghostzones, f);
    case 0b010:
      return loop<0, 1, 0>(where, group_nghostzones, f);
    case 0b011:
      return loop<1, 1, 0>(where, group_nghostzones, f);
    case 0b100:
      return loop<0, 0, 1>(where, group_nghostzones, f);
    case 0b101:
      return loop<1, 0, 1>(where, group_nghostzones, f);
    case 0b110:
      return loop<0, 1, 1>(where, group_nghostzones, f);
    case 0b111:
      return loop<1, 1, 1>(where, group_nghostzones, f);
    default:
      assert(0);
    }
  }

  template <typename F>
  void loop_idx(where_t where, const array<int, dim> &indextype,
                const F &f) const {
    loop_idx(where, indextype, nghostzones, f);
  }
};

template <typename F>
void loop_idx(const cGH *cctkGH, where_t where,
              const array<int, dim> &indextype,
              const array<int, dim> &nghostzones, const F &f) {
  GridDescBase(cctkGH).loop_idx(where, indextype, nghostzones, f);
}

template <typename F>
void loop_idx(const cGH *cctkGH, where_t where,
              const array<int, dim> &indextype, const F &f) {
  GridDescBase(cctkGH).loop_idx(where, indextype, f);
}

template <int CI, int CJ, int CK, typename F>
void loop(const cGH *cctkGH, where_t where, const F &f) {
  GridDescBase(cctkGH).loop<CI, CJ, CK>(where, f);
}

// Keep these for convenience
template <int CI, int CJ, int CK, typename F>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE void loop_all(const cGH *cctkGH,
                                                  const F &f) {
  loop<CI, CJ, CK>(cctkGH, where_t::everywhere, f);
}

template <int CI, int CJ, int CK, typename F>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE void loop_int(const cGH *cctkGH,
                                                  const F &f) {
  loop<CI, CJ, CK>(cctkGH, where_t::interior, f);
}

template <int CI, int CJ, int CK, typename F>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE void loop_bnd(const cGH *cctkGH,
                                                  const F &f) {
  loop<CI, CJ, CK>(cctkGH, where_t::boundary, f);
}

template <int CI, int CJ, int CK, typename F>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
loop_ghosts_inclusive(const cGH *cctkGH, const F &f) {
  loop<CI, CJ, CK>(cctkGH, where_t::ghosts_inclusive, f);
}

template <int CI, int CJ, int CK, typename F>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE void loop_ghosts(const cGH *cctkGH,
                                                     const F &f) {
  loop<CI, CJ, CK>(cctkGH, where_t::ghosts, f);
}

////////////////////////////////////////////////////////////////////////////////

extern bool CarpetX_poison_undefined_values;

template <typename T, int CI, int CJ, int CK> struct GF3D {
  static_assert(CI == 0 || CI == 1, "");
  static_assert(CJ == 0 || CJ == 1, "");
  static_assert(CK == 0 || CK == 1, "");
  typedef T value_type;
  static constexpr int di = 1;
  const int dj, dk, np;
  const int ni, nj, nk;
  T *restrict const ptr;
  static constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE array<int, dim> indextype() {
    return {CI, CJ, CK};
  }
  GF3D() = delete;
  GF3D(const GF3D &gf) = default;
  GF3D(GF3D &&) = default;
  GF3D &operator=(const GF3D &gf) = default;
  GF3D &operator=(GF3D &&) = default;
  GF3D(const cGH *restrict cctkGH, T *restrict ptr)
      : dj(di * (cctkGH->cctk_ash[0] + !CI)),
        dk(dj * (cctkGH->cctk_ash[1] + !CJ)),
        np(dk * (cctkGH->cctk_ash[2] + !CK)), ni(cctkGH->cctk_lsh[0] + !CI),
        nj(cctkGH->cctk_lsh[1] + !CJ), nk(cctkGH->cctk_lsh[2] + !CK), ptr(ptr) {
  }
  GF3D(const cGH *restrict cctkGH, vector<vector<T> > &buffers)
      : dj(di * (cctkGH->cctk_ash[0] + !CI)),
        dk(dj * (cctkGH->cctk_ash[1] + !CJ)),
        np(dk * (cctkGH->cctk_ash[2] + !CK)), ni(cctkGH->cctk_lsh[0] + !CI),
        nj(cctkGH->cctk_lsh[1] + !CJ), nk(cctkGH->cctk_lsh[2] + !CK),
        ptr(alloc_buffer(buffers, np)) {}

private:
  static T *restrict alloc_buffer(vector<vector<T> > &buffers, const int np) {
    if (CarpetX_poison_undefined_values)
      buffers.emplace_back(np, numeric_limits<T>::quiet_NaN());
    else
      buffers.emplace_back(np);
    return buffers.back().data();
  }

public:
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE int offset(int i, int j, int k) const {
    // These index checks prevent vectorization. We thus only enable
    // them in debug mode.
#ifdef CCTK_DEBUG
    assert(i >= 0 && i < ni);
    assert(j >= 0 && j < nj);
    assert(k >= 0 && k < nk);
#endif
    return i * di + j * dj + k * dk;
  }
  inline int offset(const vect<int, dim> &I) const {
    return offset(I[0], I[1], I[2]);
  }
  inline T &restrict operator()(int i, int j, int k) const {
    return ptr[offset(i, j, k)];
  }
  inline T &restrict operator()(const vect<int, dim> &I) const {
    return ptr[offset(I)];
  }
};

template <typename T> struct GF3D1 {
  typedef T value_type;
  T *restrict ptr;
#ifdef CCTK_DEBUG
  array<int, dim> imin, imax;
  array<int, dim> ash;
#endif
  static constexpr int di = 1;
  int dj, dk, np;
  int off;
  GF3D1() = delete;
  GF3D1(const GF3D1 &) = delete;
  GF3D1(GF3D1 &&) = default;
  GF3D1 &operator=(const GF3D1 &) = delete;
  GF3D1 &operator=(GF3D1 &&) = default;
  GF3D1(T *restrict ptr, const array<int, dim> &imin,
        const array<int, dim> &imax, const array<int, dim> &ash)
      : ptr(ptr),
#ifdef CCTK_DEBUG
        imin(imin), imax(imax), ash(ash),
#endif
        dj(di * ash[0]), dk(dj * ash[1]), np(dk * ash[2]),
        off(imin[0] * di + imin[1] * dj + imin[2] * dk) {
  }
  GF3D1(const cGH *restrict cctkGH, const array<int, dim> &indextype,
        const array<int, dim> &nghostzones, T *restrict ptr) {
    for (int d = 0; d < dim; ++d)
      assert(indextype[d] == 0 || indextype[d] == 1);
    for (int d = 0; d < dim; ++d) {
      assert(nghostzones[d] >= 0);
      assert(nghostzones[d] <= cctkGH->cctk_nghostzones[d]);
    }
    array<int, dim> imin, imax;
    for (int d = 0; d < dim; ++d) {
      imin[d] = cctkGH->cctk_nghostzones[d] - nghostzones[d];
      imax[d] = cctkGH->cctk_lsh[d] + (1 - indextype[d]) -
                (cctkGH->cctk_nghostzones[d] - nghostzones[d]);
    }
    array<int, dim> ash;
    for (int d = 0; d < dim; ++d)
      ash[d] = cctkGH->cctk_ash[d] + (1 - indextype[d]) -
               2 * (cctkGH->cctk_nghostzones[d] - nghostzones[d]);
    *this = GF3D1(ptr, imin, imax, ash);
  }
  inline int offset(int i, int j, int k) const {
    // These index checks prevent vectorization. We thus only enable
    // them in debug mode.
#ifdef CCTK_DEBUG
    assert(i >= imin[0] && i < imax[0]);
    assert(j >= imin[1] && j < imax[1]);
    assert(k >= imin[2] && k < imax[2]);
#endif
    return i * di + j * dj + k * dk - off;
  }
  inline T &restrict operator()(int i, int j, int k) const {
    return ptr[offset(i, j, k)];
  }
  inline T &restrict operator()(const vect<int, dim> &I) const {
    return ptr[offset(I[0], I[1], I[2])];
  }
  inline T &restrict operator()(const PointDesc &p) const {
    return (*this)(p.I);
  }
};

} // namespace Loop

#endif // #ifndef LOOP_HXX
