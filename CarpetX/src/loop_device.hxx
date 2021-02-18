#ifndef LOOP_DEVICE_HXX
#define LOOP_DEVICE_HXX

#include "loop.hxx"

#ifdef __CUDACC__
#define CCTK_DEVICE __device__
#define CCTK_HOST __host__
#else
#define CCTK_DEVICE
#define CCTK_HOST
#endif

#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#ifdef AMREX_USE_GPU
#include <AMReX_Gpu.H>
#endif

namespace Loop {

// The struct GridDescBaseDevice is very similar to GridDescBase, except that it
// requires AMReX include files, and thus "REQUIRES CarpetX" in the using
// thorns' configuration.ccl file.
struct GridDescBaseDevice : GridDescBase {

protected:
  GridDescBaseDevice() : GridDescBase() {}

public:
  GridDescBaseDevice(const cGH *cctkGH) : GridDescBase(cctkGH) {}

  // Loop over a given box
  template <int CI, int CJ, int CK, int VS = 1, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_box_device(const F &f, const array<int, dim> &restrict imin,
                  const array<int, dim> &restrict imax,
                  const array<int, dim> &restrict inormal) const {
    static_assert(CI == 0 || CI == 1, "");
    static_assert(CJ == 0 || CJ == 1, "");
    static_assert(CK == 0 || CK == 1, "");
    // TODO static_assert(VS is power of 2);

    for (int d = 0; d < dim; ++d)
      assert(!(imin[d] >= imax[d]));

    // For some reason, the argument inormal cannot be captured, but a copy of
    // inormal can
    const auto inormal1 = inormal;
    const auto kernel =
        [=, *this] CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST(
            const int i, const int j, const int k) {
          f(point_desc<CI, CJ, CK>(inormal1, imin[0], imax[0], i, j, k));
        };

    // array<bool, dim> bforward;
    // for (int d = 0; d < dim; ++d)
    //   bforward[d] = inormal[d] >= 0;
    // bool all_forward = true;
    // for (int d = 0; d < dim; ++d)
    //   all_forward &= bforward[d];

    // TODO: introduce loop_bnd_stepped, have stepping loop there
    constexpr array<bool, dim> bforward{false, false, false};
    constexpr bool all_forward = false;

#ifndef AMREX_USE_GPU
    // Run on CPU

    if (all_forward) {

      // const int imin0_aligned = imin[0] & ~(VS - 1);

      for (int k = imin[2]; k < imax[2]; ++k) {
        for (int j = imin[1]; j < imax[1]; ++j) {
          // #pragma omp simd
          // for (int i = imin[0]; i < imax[0]; ++i) {
          // for (int i = imin0_aligned; i < imax[0]; i += VS) {
          for (int i = imin[0]; i < imax[0]; i += VS) {
            kernel(i, j, k);
          }
        }
      }

    } else {
      // At least one direction is reversed; loop from the inside out
      assert(VS == 1);

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

#else
    // Run on GPU
    static_assert(VS == 1, "");

    // Convert to AMReX box
    const amrex::Box box(
        amrex::IntVect(imin[0], imin[1], imin[2]),
        amrex::IntVect(imax[0] - 1, imax[1] - 1, imax[2] - 1),
        amrex::IntVect(CI ? amrex::IndexType::CELL : amrex::IndexType::NODE,
                       CJ ? amrex::IndexType::CELL : amrex::IndexType::NODE,
                       CK ? amrex::IndexType::CELL : amrex::IndexType::NODE));
    amrex::launch(box, [=] CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE(
                           const amrex::Box &box) {
#if 0
      const amrex::IntVect bmin = box.smallEnd();
      const amrex::IntVect bend = box.bigEnd();
      const array<int, dim> imin{bmin[0], bmin[1], bmin[2]};
      const array<int, dim> imax{bend[0] + 1, bend[1] + 1, bend[2] + 1};

      for (int k = imin[2]; k < imax[2]; ++k) {
        for (int j = imin[1]; j < imax[1]; ++j) {
          for (int i = imin[0]; i < imax[0]; ++i) {
            kernel(i, j, k);
          }
        }
      }
#endif
      assert(box.bigEnd()[0] == box.smallEnd()[0] &&
             box.bigEnd()[1] == box.smallEnd()[1] &&
             box.bigEnd()[2] == box.smallEnd()[2]);
      const int i = box.smallEnd()[0];
      const int j = box.smallEnd()[1];
      const int k = box.smallEnd()[2];
      kernel(i, j, k);
    });

#endif
  }

  // Loop over all points
  template <int CI, int CJ, int CK, int VS = 1, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_all_device(const array<int, dim> &group_nghostzones, const F &f) const {
    const array<int, dim> offset{!CI, !CJ, !CK};
    array<int, dim> imin, imax;
    for (int d = 0; d < dim; ++d) {
      int ghost_offset = nghostzones[d] - group_nghostzones[d];
      imin[d] = std::max(tmin[d], ghost_offset);
      imax[d] = std::min(tmax[d], lsh[d] - offset[d] - ghost_offset);
    }
    const array<int, dim> inormal{0, 0, 0};

    loop_box_device<CI, CJ, CK, VS>(f, imin, imax, inormal);
  }

  // Loop over all interior points
  template <int CI, int CJ, int CK, int VS = 1, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_int_device(const array<int, dim> &group_nghostzones, const F &f) const {
    const array<int, dim> offset{!CI, !CJ, !CK};
    array<int, dim> imin, imax;
    for (int d = 0; d < dim; ++d) {
      imin[d] = std::max(tmin[d], nghostzones[d]);
      imax[d] = std::min(tmax[d], lsh[d] - offset[d] - nghostzones[d]);
    }
    const array<int, dim> inormal{0, 0, 0};

    loop_box_device<CI, CJ, CK, VS>(f, imin, imax, inormal);
  }

  // Loop over all outer boundary points. This excludes ghost faces, but
  // includes ghost edges/corners on non-ghost faces. Loop over faces first,
  // then edges, then corners.
  template <int CI, int CJ, int CK, int VS = 1, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_bnd_device(const array<int, dim> &group_nghostzones, const F &f) const {
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
                  const int end_int = lsh[d] - offset[d] - nghostzones[d];
                  const int end_bnd = lsh[d] - offset[d] - ghost_offset;
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
                  imax[d] = std::min(tmax[d], imax[d]);
                }

                loop_box_device<CI, CJ, CK, VS>(f, imin, imax, inormal);
              }
            } // if rank
          }
        }
      }

    } // for rank
  }
};

} // namespace Loop

#endif // #ifndef LOOP_DEVICE_HXX
