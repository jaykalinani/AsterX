#ifndef LOOP_DEVICE_HXX
#define LOOP_DEVICE_HXX

#include "loop.hxx"

#include <cctk.h>

#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#ifdef AMREX_USE_GPU
#include <AMReX_Gpu.H>
#endif

#include <cassert>

#ifndef HAVE_CAPABILITY_AMReX
#error                                                                         \
    "Using #include <loop_device.hxx> requires the capability 'AMReX' in the calling thorn. Add this to the thorn's 'configuration.ccl' file."
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
#ifndef AMREX_USE_GPU
    return this->template loop_box<CI, CJ, CK, VS>(f, imin, imax, inormal);
#else

    static_assert(CI == 0 || CI == 1, "");
    static_assert(CJ == 0 || CJ == 1, "");
    static_assert(CK == 0 || CK == 1, "");
    static_assert(VS > 0, "");
    // TODO static_assert(VS is power of 2);

    for (int d = 0; d < dim; ++d)
      assert(!(imin[d] >= imax[d]));

    // array<bool, dim> bforward;
    // for (int d = 0; d < dim; ++d)
    //   bforward[d] = inormal[d] >= 0;
    // bool all_forward = true;
    // for (int d = 0; d < dim; ++d)
    //   all_forward &= bforward[d];

    // Run on GPU
    static_assert(VS == 1, "Only vertex centered code is supported on GPUs");

    // For some reason, the arguments imin, imax, and inormal cannot
    // be captured correctly in CUDA, but copies of them can
    const auto imin1 = imin;
    const auto imax1 = imax;
    const auto inormal1 = inormal;

    // Convert to AMReX box
    const amrex::Box box(
        amrex::IntVect(imin[0], imin[1], imin[2]),
        amrex::IntVect(imax[0] - 1, imax[1] - 1, imax[2] - 1),
        amrex::IntVect(CI ? amrex::IndexType::CELL : amrex::IndexType::NODE,
                       CJ ? amrex::IndexType::CELL : amrex::IndexType::NODE,
                       CK ? amrex::IndexType::CELL : amrex::IndexType::NODE));
    amrex::launch(box, [=, *this] CCTK_DEVICE(
                           const amrex::Box &box) CCTK_ATTRIBUTE_ALWAYS_INLINE {
#ifdef CCTK_DEBUG
      assert(box.bigEnd()[0] == box.smallEnd()[0] &&
             box.bigEnd()[1] == box.smallEnd()[1] &&
             box.bigEnd()[2] == box.smallEnd()[2]);
#endif
      const int i = box.smallEnd()[0];
      const int j = box.smallEnd()[1];
      const int k = box.smallEnd()[2];
      f(point_desc<CI, CJ, CK>(inormal1, imin1[0], imax1[0], i, j, k));
    });
#endif

#ifdef AMREX_USE_GPU
    static bool have_gpu_sync_after_every_kernel = false;
    static bool gpu_sync_after_every_kernel;
    if (!have_gpu_sync_after_every_kernel) {
      int type;
      const void *const gpu_sync_after_every_kernel_ptr =
          CCTK_ParameterGet("gpu_sync_after_every_kernel", "CarpetX", &type);
      assert(gpu_sync_after_every_kernel_ptr);
      assert(type == PARAMETER_BOOLEAN);
      gpu_sync_after_every_kernel =
          *static_cast<const CCTK_INT *>(gpu_sync_after_every_kernel_ptr);
    }
    if (gpu_sync_after_every_kernel) {
      amrex::Gpu::synchronize();
      AMREX_GPU_ERROR_CHECK();
    }
#endif
  }

  // Loop over all points
  template <int CI, int CJ, int CK, int VS = 1, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_all_device(const array<int, dim> &group_nghostzones, const F &f) const {
    const array<int, dim> offset{CI, CJ, CK};
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
    const array<int, dim> offset{CI, CJ, CK};
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
    const array<int, dim> offset{CI, CJ, CK};

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

  // Loop over all outer ghost points. This excludes ghost edges/corners on
  // non-ghost faces. Loop over faces first, then edges, then corners.
  template <int CI, int CJ, int CK, int VS = 1, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_ghosts_device(const array<int, dim> &group_nghostzones,
                     const F &f) const {
    const array<int, dim> offset{CI, CJ, CK};

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
