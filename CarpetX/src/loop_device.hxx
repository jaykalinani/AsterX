#ifndef LOOP_DEVICE_HXX
#define LOOP_DEVICE_HXX

#include "loop.hxx"

#warning "TODO: try removing these"
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
  template <int CI, int CJ, int CK, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_box_device(const F &f, const array<int, dim> &restrict imin,
                  const array<int, dim> &restrict imax,
                  const array<int, dim> &restrict inormal) const {
    static_assert(CI == 0 || CI == 1, "");
    static_assert(CJ == 0 || CJ == 1, "");
    static_assert(CK == 0 || CK == 1, "");

    for (int d = 0; d < dim; ++d)
      assert(!(imin[d] >= imax[d]));

    constexpr int di = 1;
    const int dj = di * (ash[0] + !CI);
    const int dk = dj * (ash[1] + !CJ);

    // For some reason, the argument inormal cannot be captured, but a copy of
    // inormal can
    const auto inormal1 = inormal;
    const auto kernel =
        [=, *this] CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST(
            const int i, const int j, const int k) {
          const CCTK_REAL x =
              x0[0] + (lbnd[0] + i - CCTK_REAL(!CI) / 2) * dx[0];
          const CCTK_REAL y =
              x0[1] + (lbnd[1] + j - CCTK_REAL(!CJ) / 2) * dx[1];
          const CCTK_REAL z =
              x0[2] + (lbnd[2] + k - CCTK_REAL(!CK) / 2) * dx[2];
          const int idx = i * di + j * dj + k * dk;
          const Loop::PointDesc p{i,         j,        k,         x,   y,  z,
                                  dx[0],     dx[1],    dx[2],     idx, dj, dk,
                                  {i, j, k}, inormal1, {x, y, z}, dx};
          f(p);
        };

    array<bool, dim> bforward;
    for (int d = 0; d < dim; ++d)
      bforward[d] = inormal[d] >= 0;
    bool all_forward = true;
    for (int d = 0; d < dim; ++d)
      all_forward &= bforward[d];
    assert(all_forward);

#ifndef AMREX_USE_GPU
    // Run on CPU

    for (int k = imin[2]; k < imax[2]; ++k) {
      for (int j = imin[1]; j < imax[1]; ++j) {
#pragma omp simd
        for (int i = imin[0]; i < imax[0]; ++i) {
          kernel(i, j, k);
        }
      }
    }

#else
    // Run on GPU

    // Convert to AMReX box
    const amrex::Box box(
        amrex::IntVect(imin[0], imin[1], imin[2]),
        amrex::IntVect(imax[0] - 1, imax[1] - 1, imax[2] - 1),
        amrex::IntVect(CI ? amrex::IndexType::CELL : amrex::IndexType::NODE,
                       CJ ? amrex::IndexType::CELL : amrex::IndexType::NODE,
                       CK ? amrex::IndexType::CELL : amrex::IndexType::NODE));
    amrex::launch(box, [=] CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE(
                           const amrex::Box &box) {
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
    });
#endif
  }

  // Loop over all interior points
  template <int CI, int CJ, int CK, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_int_device(const array<int, dim> &group_nghostzones, const F &f) const {
    const array<int, dim> offset{!CI, !CJ, !CK};
    array<int, dim> imin, imax;
    for (int d = 0; d < dim; ++d) {
      imin[d] = std::max(tmin[d], nghostzones[d]);
      imax[d] = std::min(tmax[d] + (tmax[d] >= lsh[d] ? offset[d] : 0),
                         lsh[d] + offset[d] - nghostzones[d]);
    }
    const array<int, dim> inormal{0, 0, 0};

    loop_box_device<CI, CJ, CK>(f, imin, imax, inormal);
  }
};

} // namespace Loop

#endif // #ifndef LOOP_DEVICE_HXX
