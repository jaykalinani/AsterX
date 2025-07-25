#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>
#include <limits>

#include "aster_utils.hxx"

namespace AsterAnalysis {
using namespace std;
using namespace Loop;
using namespace EOSX;
using namespace AsterUtils;

extern "C" void AsterAnalysis_MHD(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterAnalysis_MHD;
  DECLARE_CCTK_PARAMETERS;

  /* handy aliases for grid functions */
  const vec<GF3D2<const CCTK_REAL>, 3> gf_beta{betax, betay, betaz};
  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};
  const vec<GF3D2<const CCTK_REAL>, 3> gf_Avec{Avec_x, Avec_y, Avec_z};
  const vec<GF3D2<const CCTK_REAL>, 3> gf_Bvec{Bvecx, Bvecy, Bvecz};
  const vec<GF3D2<const CCTK_REAL>, 3> gf_vel{velx, vely, velz};

  /* inverse grid spacings (assumed uniform) */
  const CCTK_REAL idx = 1.0 / CCTK_DELTA_SPACE(0);
  const CCTK_REAL idy = 1.0 / CCTK_DELTA_SPACE(1);
  const CCTK_REAL idz = 1.0 / CCTK_DELTA_SPACE(2);

  /* main device loop over interior + ghost zones */
  grid.loop_all_device<1, 1, 1>(grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

        /* ---- metric averages at cell center ---- */
        const CCTK_REAL alp_avg = calc_avg_v2c(alp, p);
        const vec<CCTK_REAL, 3> beta_avg([&](int i) ARITH_INLINE {
          return calc_avg_v2c(gf_beta(i), p);
        });
        const smat<CCTK_REAL, 3> g_avg([&](int i, int j) ARITH_INLINE {
          return calc_avg_v2c(gf_g(i, j), p);
        });
        const CCTK_REAL detg = calc_det(g_avg);
        const CCTK_REAL sqrt_detg = sqrt(detg);

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

        /* ---- magnetic field ---- */
        const vec<CCTK_REAL, 3> B_up{
            Bvecx(p.I),
            Bvecy(p.I),
            Bvecz(p.I)};
        const vec<CCTK_REAL, 3> B_low = calc_contraction(g_avg, B_up);
        const CCTK_REAL B2big = calc_contraction(B_up, B_low);
        B_norm(p.I) = sqrt(B2big);

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

        /* inverse Beta and magnetization */
        mhd_press_ratio(p.I) = b2small(p.I)/(2.0*press(p.I));
        mhd_energy_ratio(p.I) = b2small(p.I)/(2.0*rho(p.I));
        
        /* ---- shift & four‑velocity ---- */
        const vec<CCTK_REAL, 3> beta_avg_low = calc_contraction(g_avg, beta_avg);
        const CCTK_REAL beta2 = calc_contraction(beta, beta_avg_low);

        const vec<CCTK_REAL, 3> u = velshift * w_lorentz(p.I);                     /* u^i */
        const CCTK_REAL utlow = (-alp_avg + calc_contraction(v_up, beta_avg_low)) * w_lorentz(p.I); /* u_0 */

        /* b_i & b_0 */
        const vec<CCTK_REAL, 3> bs_low = bs0 * beta_avg_low + calc_contraction(g_avg, bs);
        const CCTK_REAL bstlow = bs0 * (-pow2(alp_avg) + beta2) + calc_contraction(bs, beta_avg_low);

        /* ---- Poynting flux S^i ---- */
        const vec<CCTK_REAL, 3> S = -alp_avg * (b2small(p.I) * u * utlow - bs * bstlow);
        poynting_flux_x(p.I) = S(0);
        poynting_flux_y(p.I) = S(1);
        poynting_flux_z(p.I) = S(2);
        const vec<CCTK_REAL, 3> S_low = calc_contraction(g_avg, S);
        poynting_flux_norm(p.I) = calc_contraction(S, S_low);

        /* ---- alternative Poynting (E x B) ---- */
        const vec<CCTK_REAL, 3> vtilde_low = alp_avg * v_low - beta_avg_low;
        const vec<CCTK_REAL, 3> E{
            -(vtilde_low(1) * B_low(2) - vtilde_low(2) * B_low(1)),
            -(vtilde_low(2) * B_low(0) - vtilde_low(0) * B_low(2)),
            -(vtilde_low(0) * B_low(1) - vtilde_low(1) * B_low(0))};
        const vec<CCTK_REAL, 3> S_alt{
            E(1) * B_up(2) - E(2) * B_up(1),
            E(2) * B_up(0) - E(0) * B_up(2),
            E(0) * B_up(1) - E(1) * B_up(0)};
        poynting_vector_x(p.I) = S_alt(0);
        poynting_vector_y(p.I) = S_alt(1);
        poynting_vector_z(p.I) = S_alt(2);

        /* normal to sphere */
        const CCTK_REAL rx = x(p.I);
        const CCTK_REAL ry = y(p.I);
        const CCTK_REAL rz = z(p.I);
        const CCTK_REAL rmag = r(p.I) > 0.0 ? r(p.I) : sqrt(rx * rx + ry * ry + rz * rz);
        const vec<CCTK_REAL, 3> n{rx / rmag, ry / rmag, rz / rmag};
        poynting_scalar(p.I) = calc_contraction(S_alt, n);

        /* ---- Lambda_MRI ---- */
        const CCTK_REAL omega_x = alp_avg * velx(p.I) - beta_avg(0);
        const CCTK_REAL omega_y = alp_avg * vely(p.I) - beta_avg(1);
        CCTK_REAL local_omega = 0.0;
        if (pow2(rx) + pow2(ry) > 0.0) {
          local_omega = (rx * omega_y - ry * omega_x) / (pow2(rx) + pow2(ry));
        }

        const CCTK_REAL Wtot = rho(p.I) * (1.0 + eps(p.I)) + press(p.I) + b2;
        auto lambda_MRI = [&](CCTK_REAL bi, CCTK_REAL bi_low, CCTK_REAL &lambda) {
          const CCTK_REAL va = sqrt(bi * bi_low / Wtot);
          lambda = (va > numeric_limits<CCTK_REAL>::epsilon() && fabs(local_omega) > 1.0e-15)
                  ? 2.0 * M_PI * fabs(va / local_omega)
                  : 0.0;
        };
        lambda_MRI(bs(0), bs_low(0), lambda_MRI_x(p.I));
        lambda_MRI(bs(1), bs_low(1), lambda_MRI_y(p.I));
        lambda_MRI(bs(2), bs_low(2), lambda_MRI_z(p.I));

        /* ---- divB & divA ---- */
        const auto off = [&](int di,int dj,int dk){ return p.offset(di,dj,dk); };

        /* B is zone‑centred so plain centred FD */
        const CCTK_REAL divB_val =
              0.5*idx*(dBx(off( 1,0,0)) - dBx(off(-1,0,0))) +
              0.5*idy*(dBy(off(0, 1,0)) - dBy(off(0,-1,0))) +
              0.5*idz*(dBz(off(0,0, 1)) - dBz(off(0,0,-1)));
        divB(p.I) = divB_val;

        /* A_x lives on x‑edges, etc.  Same integer offsets pick ±½ edges */
        const CCTK_REAL divA_val =
              0.5*idx*(Avec_x(off( 1,0,0)) - Avec_x(off(-1,0,0))) +
              0.5*idy*(Avec_y(off(0, 1,0)) - Avec_y(off(0,-1,0))) +
              0.5*idz*(Avec_z(off(0,0, 1)) - Avec_z(off(0,0,-1)));
        divA(p.I) = divA_val;

      }); /* end grid loop */
}

} /* namespace AsterAnalysis */

