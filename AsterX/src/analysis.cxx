#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include "aster_utils.hxx"

namespace AsterX {
using namespace std;
using namespace Loop;
using namespace Arith;
using namespace AsterUtils;

extern "C" void AsterX_Analysis(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Analysis;
  DECLARE_CCTK_PARAMETERS;

  const vec<GF3D2<const CCTK_REAL>, 3> gf_beta{betax, betay, betaz};
  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};

  /* Loop over the entire grid (0 to n-1 cells in each direction) */
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        /* Computing metric components at cell centers */
        const CCTK_REAL alp_avg = calc_avg_v2c(alp, p);
        const vec<CCTK_REAL, 3> beta_avg(
            [&](int i) ARITH_INLINE { return calc_avg_v2c(gf_beta(i), p); });
        const smat<CCTK_REAL, 3> g_avg([&](int i, int j) ARITH_INLINE {
          return calc_avg_v2c(gf_g(i, j), p);
        });

        /* Compute w_lor */
        const vec<CCTK_REAL, 3> v_up{velx(p.I), vely(p.I), velz(p.I)};
        const vec<CCTK_REAL, 3> v_low = calc_contraction(g_avg, v_up);
        const vec<CCTK_REAL, 3> velshift = v_up - beta_avg / alp_avg;
        const CCTK_REAL w_lorentz = calc_wlorentz(v_low, v_up);

        /* cell-centered B^i and B_i */
        const vec<CCTK_REAL, 3> B_up{Bvecx(p.I), Bvecy(p.I), Bvecz(p.I)};
        const vec<CCTK_REAL, 3> B_low = calc_contraction(g_avg, B_up);
        /* b^0 */
        const CCTK_REAL bs0 =
            w_lorentz * calc_contraction(B_up, v_low) / alp_avg;
        /* b^2 */
        b2small(p.I) = (calc_contraction(B_up, B_low) + pow2(alp_avg * bs0)) /
                       pow2(w_lorentz);

        total_energy(p.I) = rho(p.I) * (1.0 + eps(p.I)) + b2small(p.I) / 2.0;
      }); // end of loop over grid
}

} // namespace AsterX
