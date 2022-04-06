#include <fixmath.hxx>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cmath>

namespace GRHydroToyGPU {
using namespace std;

// Computes the determinant of spatial metric
template <typename T>
CCTK_DEVICE CCTK_HOST const T calc_detg(const T &gxx, const T &gxy, const T &gxz,
		                  const T &gyy, const T &gyz, const T &gzz) {
	return -gxz * gxz * gyy + 2.0 * gxy * gxz * gyz -
                gxx * gyz * gyz - gxy * gxy * gzz +
                gxx * gyy * gzz;
  }

// Computes the upper spatial metric components
template <typename T>
CCTK_DEVICE CCTK_HOST const array<T, 6> calc_upperg(const T &gxx, const T &gxy, const T &gxz,
                                  const T &gyy, const T &gyz, const T &gzz, const T &detg) {
	return {(-gyz*gyz + gyy*gzz)/detg,  //uxx
                (gxz*gyz - gxy*gzz)/detg,   //uxy
                (-gxz*gxz + gxx*gzz)/detg,  //uyy
                (-gxz*gyy + gxy*gyz)/detg,  //uxz
                (gxy*gxz - gxx*gyz)/detg,   //uyz
                (-gxy*gxy + gxx*gyy)/detg   //uzz
               };
  }

// Second order finite difference
template <typename T>
CCTK_DEVICE CCTK_HOST const T calc_fd2(const T &qp1, const T &qm1, const T &dx) {
        return (0.5/dx)*(qp1-qm1);
  }

// Fourth order finite difference
template <typename T>
CCTK_DEVICE CCTK_HOST const T calc_fd4(const T &qp2, const T &qp1, const T &qm1, 
		                   const T &qpm2, const T &dx) {
        return (1.0/(12.0*dx))*(-qp2 + 8.0*qp1 - 8.0*qm1 + qm2);
  }

}// namespace GRHydroToyGPU
