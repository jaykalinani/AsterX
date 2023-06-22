#ifndef RECONX_MP5_HXX
#define RECONX_MP5_HXX

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>

#include "reconx_utils.hxx"
#include "minmod.hxx"

namespace ReconX {

using namespace std;
using namespace Arith;
using namespace Loop;


// Compute the median of three numbers
template<typename T> inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
T median(T &x, T &y, T &z) {
    return x + minmod(y-x, z-x);
}

/* Fifth-order monotonicity preserving (MP5) scheme of Suresh and Huynh 1997.
Paper: "Accurate monotonicity-preserving schemes with Runge–Kutta time stepping, J. Comput. Phys. 136 (1997) 83–99."
*/

inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST CCTK_REAL
mp5(const GF3D2<const CCTK_REAL> &gf_var,
      const array<const vect<int, dim>, 5> &cells,
      const reconstruct_params_t &reconstruct_params) {

  // Unpack all MP5 parameters
  const CCTK_REAL mp5_alpha = reconstruct_params.mp5_alpha;

  // Unpack all cells in the stencil
  const auto &Imm = cells[0];
  const auto &Im = cells[1];
  const auto &I = cells[2];
  const auto &Ip = cells[3];
  const auto &Ipp = cells[4];

  // Grid function at neighboring cells
  const CCTK_REAL &gf_Imm = gf_var(Imm);
  const CCTK_REAL &gf_Im = gf_var(Im);
  const CCTK_REAL &gf_I = gf_var(I);
  const CCTK_REAL &gf_Ip = gf_var(Ip);
  const CCTK_REAL &gf_Ipp = gf_var(Ipp);

  const CCTK_REAL ul = (2*gf_Imm - 13*gf_Im + 47*gf_I +
                    27*gf_Ip - 3*gf_Imm) / 60.0; 

  const CCTK_REAL deltam = gf_I-gf_Im;
  const CCTK_REAL deltap = gf_Ip-gf_I;

  const CCTK_REAL ump = gf_I + minmod(deltap, mp5_alpha*deltam );
  
  if((ul - gf_I)*(ul - ump) < 0) {
     return ul;
  }
  else {
     const CCTK_REAL dm = gf_Imm + gf_I - 2*gf_Im;
     const CCTK_REAL d = gf_Im + gf_Ip - 2*gf_I;
     const CCTK_REAL dp = gf_I + gf_Ipp - 2*gf_Ip;

     const CCTK_REAL dmp = minmod(minmod(4*d-dp, 4*dp-d),
                        minmod(d, dp));
     const CCTK_REAL dmm = minmod(minmod(4*dm - d, 4*d - dm),
                        minmod(dm, d));

     const CCTK_REAL ulc = gf_I + 0.5*deltam + 4.0/3.0*dmm;
     const CCTK_REAL umd = 0.5*(gf_I + gf_Ip) - 0.5*dmp;

     const CCTK_REAL uul = gf_I + mp5_alpha*deltam;

     const CCTK_REAL umin = std::max(
                        std::min(gf_I, std::min(gf_Ip, umd)),
                        std::min(gf_I, std::min(uul, ulc)));
     const CCTK_REAL umax = std::min(
                        std::max(gf_I, std::max(gf_Ip, umd)),
                        std::max(gf_I, std::max(uul, ulc)));

     return median(umin, ul, umax);

  }
  return 0;

}

} // namespace ReconX

#endif // RECONX_MP5_HXX
