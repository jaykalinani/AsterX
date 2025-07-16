#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cmath>

#include "aster_utils.hxx"

namespace AsterAnalysis {
using namespace std;
using namespace Loop;
using namespace EOSX;
using namespace AsterUtils;

extern "C" void AsterAnalysis_MHD(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterAnalysis_MHD;
  DECLARE_CCTK_PARAMETERS;

  /* grid functions */
  const vec<GF3D2<const CCTK_REAL>, 3> gf_beta{betax, betay, betaz};
  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};
  const vec<GF3D2<const CCTK_REAL>, 3> gf_Avec{Avec_x, Avec_y, Avec_z};

  /* Loop over the entire grid (0 to n-1 cells in each direction) */
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

       /* Computing metric components at cell centers */
       const CCTK_REAL alp_avg = calc_avg_v2c(alp, p);
       const vec<CCTK_REAL, 3> beta_avg(
            [&](int i) ARITH_INLINE { return calc_avg_v2c(gf_beta(i), p); });
       const smat<CCTK_REAL, 3> g_avg([&](int i, int j) ARITH_INLINE {
          return calc_avg_v2c(gf_g(i, j), p);
       });

       /* Determinant of spatial metric */
       const CCTK_REAL detg = calc_det(g_avg);
       const CCTK_REAL sqrt_detg = sqrt(detg);
       
       /* Upper metric */
       const smat<CCTK_REAL, 3> ug_avg = calc_inv(g_avg, detg);

       /* compute w_lorentz */
        vec<CCTK_REAL, 3> v_up;
        vec<CCTK_REAL, 3> v_low;
        if (use_v_vec) {

          v_up(0) = velx(p.I);
          v_up(1) = vely(p.I);
          v_up(2) = velz(p.I);
          v_low = calc_contraction(g_avg, v_up);
          w_lorentz(p.I) = calc_wlorentz(v_low, v_up);

        } else {

          const vec<CCTK_REAL, 3> z_up{zvec_x(p.I), zvec_y(p.I), zvec_z(p.I)};
          const vec<CCTK_REAL, 3> z_low = calc_contraction(g_avg, z_up);
          w_lorentz(p.I) = calc_wlorentz_zvec(z_low, z_up);
          v_up = z_up / w_lorentz(p.I);
          v_low = z_low / w_lorentz(p.I);
        }
  
        /* v^i - beta^i / alpha */
        const vec<CCTK_REAL, 3> velshift = v_up - beta_avg / alp_avg;

        /* cell-centered B^i and B_i */
        const vec<CCTK_REAL, 3> B_up{Bvecx(p.I), Bvecy(p.I), Bvecz(p.I)};
        const vec<CCTK_REAL, 3> B_low = calc_contraction(g_avg, B_up);
        
        B_norm(p.I) = sqrt(calc_contraction(B_low, B_up));

        /* cell-centered A^i and A_i */
        const vec<CCTK_REAL, 3> A_low(
            [&](int i) ARITH_INLINE { return calc_avg_e2c(gf_Avec(i), p, i); });
        const vec<CCTK_REAL, 3> A_up = calc_contraction(ug_avg, A_up);
        
        A_norm(p.I) = sqrt(calc_contraction(A_low, A_up));
 
        /* b^mu: b^0 and b^i */
        const CCTK_REAL bs0 =
            w_lorentz(p.I) * calc_contraction(B_up, v_low) / alp_avg;
        const vec<CCTK_REAL, 3> bs =
            (B_up / w_lorentz(p.I) + alp_avg * bs0 * velshift);

        /* b^2 */
        b2small(p.I) =
            (calc_contraction(B_up, B_low) + pow2(alp_avg * bs0)) /
            pow2(w_lorentz(p.I));

        mhd_press_ratio(p.I) = b2small(p.I)/(2.0*press(p.I));
        mhd_energy_ratio(p.I) = b2small(p.I)/(2.0*rho(p.I));

        

  }); // end of loop over grid

}

} // namespace AsterAnalysis
