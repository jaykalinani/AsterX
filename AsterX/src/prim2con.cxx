#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <loop_device.hxx>

#include "prim2con.hxx"

namespace AsterX {
using namespace Loop;

extern "C" void AsterX_Prim2Con_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Prim2Con_Initial;
  DECLARE_CCTK_PARAMETERS;

  // Loop over the entire grid (0 to n-1 cells in each direction)
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Interpolate metric terms from vertices to center
        const smat<CCTK_REAL, 3> g{calc_avg_v2c(gxx, p), calc_avg_v2c(gxy, p),
                                   calc_avg_v2c(gxz, p), calc_avg_v2c(gyy, p),
                                   calc_avg_v2c(gyz, p), calc_avg_v2c(gzz, p)};

        prim pv;
        pv.rho = rho(p.I);
        pv.vel(0) = velx(p.I);
        pv.vel(1) = vely(p.I);
        pv.vel(2) = velz(p.I);
        pv.eps = eps(p.I);
        pv.press = press(p.I);
        pv.Bvec(0) = Bvecx(p.I);
        pv.Bvec(1) = Bvecy(p.I);
        pv.Bvec(2) = Bvecz(p.I);

        cons cv;
        prim2con(g, pv, cv);

        dens(p.I) = cv.dens;
        momx(p.I) = cv.mom(0);
        momy(p.I) = cv.mom(1);
        momz(p.I) = cv.mom(2);
        tau(p.I) = cv.tau;
        dBx(p.I) = cv.dBvec(0);
        dBy(p.I) = cv.dBvec(1);
        dBz(p.I) = cv.dBvec(2);

        DEnt(p.I) = entropy(p.I)*cv.dens;

        saved_rho(p.I) = pv.rho;
        saved_velx(p.I) = pv.vel(0);
        saved_vely(p.I) = pv.vel(1);
        saved_velz(p.I) = pv.vel(2);
        saved_eps(p.I) = pv.eps;

	const vec<CCTK_REAL, 3> v_up{pv.vel(0),pv.vel(1),pv.vel(2)};
        const vec<CCTK_REAL, 3> v_low = calc_contraction(g, v_up);
        CCTK_REAL wlor = calc_wlorentz(v_low, v_up);

	zvec_x(p.I) = wlor * pv.vel(0);
	zvec_y(p.I) = wlor * pv.vel(1);
	zvec_z(p.I) = wlor * pv.vel(2);

	svec_x(p.I) = (pv.rho+pv.rho*pv.eps+pv.press)*wlor*wlor*pv.vel(0);
	svec_y(p.I) = (pv.rho+pv.rho*pv.eps+pv.press)*wlor*wlor*pv.vel(1);
	svec_z(p.I) = (pv.rho+pv.rho*pv.eps+pv.press)*wlor*wlor*pv.vel(2);

      });

  /* Initilaize Psi to 0.0 */
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Psi(p.I) = 0.0; });
}

} // namespace AsterX
