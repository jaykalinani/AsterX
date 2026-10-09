#include <AMReX.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuMemory.H>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

#include "util_Table.h"
#include "seeds_utils.hxx"
#include "setup_eos.hxx"

namespace AsterSeeds {
using namespace std;
using namespace Loop;
using namespace AsterUtils;
using namespace EOSX;

extern "C" void AsterSeeds_InterpolateNSVelocity(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterSeeds_InterpolateNSVelocity;
  DECLARE_CCTK_PARAMETERS;

  // Get NS velocities
  const int nPoints = 2;
  const int nInputArrays = 3;
  CCTK_REAL nsx[nPoints] = {CoM_NS1[0], CoM_NS2[0]};
  CCTK_REAL nsy[nPoints] = {CoM_NS1[1], CoM_NS2[1]};
  CCTK_REAL nsz[nPoints] = {CoM_NS1[2], CoM_NS2[2]};
  const void *interp_coords[nInputArrays] = {
      (const void *)nsx, (const void *)nsy, (const void *)nsz};
  const CCTK_INT inputArrayIndices[nInputArrays] = {
      CCTK_VarIndex("HydroBaseX::velx"), CCTK_VarIndex("HydroBaseX::vely"),
      CCTK_VarIndex("HydroBaseX::velz")};
  CCTK_REAL nsvx[nPoints], nsvy[nPoints], nsvz[nPoints];
  CCTK_POINTER outputArrays[nInputArrays] = {(void *)nsvx, (void *)nsvy,
                                             (void *)nsvz};

  // DriverInterpolate arguments that aren't currently used
  const int coordSystemHandle = 0;
  const CCTK_INT interpCoordsTypeCode = 0;
  const CCTK_INT outputArrayTypes[nInputArrays] = {0, 0, 0};

  const int interpHandle = CCTK_InterpHandle("CarpetX");
  if (interpHandle < 0) {
    CCTK_WARN(CCTK_WARN_ALERT, "Can't get interpolation handle");
    return;
  }

  // Create parameter table for interpolation
  const int paramTableHandle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
  if (paramTableHandle < 0) {
    CCTK_VERROR("Can't create parameter table: %d", paramTableHandle);
  }

  // Set interpolation order in the parameter table
  int ierr = Util_TableSetInt(paramTableHandle, 1, "order");
  if (ierr < 0) {
    CCTK_VERROR("Can't set order in parameter table: %d", ierr);
  }

  // Perform the interpolation
  ierr = DriverInterpolate(cctkGH, 3, interpHandle, paramTableHandle,
                           coordSystemHandle, nPoints, interpCoordsTypeCode,
                           interp_coords, nInputArrays, inputArrayIndices,
                           nInputArrays, outputArrayTypes, outputArrays);

  CCTK_VINFO("Interpolated (%g, %g, %g) as NS1 velocity", nsvx[0], nsvy[0],
             nsvz[0]);
  CCTK_VINFO("Interpolated (%g, %g, %g) as NS2 velocity", nsvx[1], nsvy[1],
             nsvz[1]);

  vel_NS1[0] = nsvx[0];
  vel_NS1[1] = nsvy[0];
  vel_NS1[2] = nsvz[0];
  vel_NS2[0] = nsvx[1];
  vel_NS2[1] = nsvy[1];
  vel_NS2[2] = nsvz[1];

  return;
}

extern "C" void AsterSeeds_SetInitialBetaFloor(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterSeeds_SetInitialBetaFloor;
  DECLARE_CCTK_PARAMETERS;

  auto eos_3p_tab3d = global_eos_3p_tab3d;
  if (not CCTK_EQUALS(evolution_eos, "Tabulated3d")) {
    CCTK_VERROR("Invalid evolution EOS type '%s'. Please, set "
                "EOSX::evolution_eos = \"Tabulated3d\" in your parameter file.",
                evolution_eos);
  }

  amrex::Gpu::DeviceScalar<unsigned int> failures(0);
  auto *failed = failures.dataPtr();
  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const CCTK_REAL pressL = press(p.I);
        const CCTK_REAL rhoL = rho(p.I);
        const CCTK_REAL tempL = temperature(p.I);
        const CCTK_REAL YeL = Ye(p.I);
        const CCTK_REAL epsL = eps(p.I);
        const CCTK_REAL entL = entropy(p.I);

        // Compute b^2
        /* Get covariant metric */
        const smat<CCTK_REAL, 3> glo([&](int i, int j) ARITH_INLINE {
          return calc_avg_v2c(gf_g(i, j), p);
        });

        vec<CCTK_REAL, 3> B_up{Bvecx(p.I), Bvecy(p.I), Bvecz(p.I)};
        vec<CCTK_REAL, 3> B_low = calc_contraction(glo, B_up);

        vec<CCTK_REAL, 3> v_up{velx(p.I), vely(p.I), velz(p.I)};
        vec<CCTK_REAL, 3> v_low = calc_contraction(glo, v_up);

        const CCTK_REAL wlor = calc_wlorentz(v_up, v_low);
        const CCTK_REAL alp_b0 = wlor * calc_contraction(B_up, v_low);

        const CCTK_REAL B2 = calc_contraction(B_up, B_low);
        const CCTK_REAL bsq = (B2 + alp_b0 * alp_b0) / (wlor * wlor);

        // Increase P if necessary
        CCTK_REAL press_lim = initial_beta_min * bsq * 0.5;

        if ((pressL >= press_lim) || (rhoL > initial_beta_rhocut)) {
          press(p.I) = pressL;
          rho(p.I) = rhoL;
          eps(p.I) = epsL;
          entropy(p.I) = entL;
        } else {
          if (!set_beta_floor(eos_3p_tab3d, press_lim, tempL, YeL,
                              rho(p.I), eps(p.I), press(p.I), entropy(p.I))) {
            amrex::HostDevice::Atomic::Add(failed, 1U);
            return;
          }
        }

        // TODO: The coorbiting velocity feature is not well tested. Use with
        // caution.
        if (set_coorbiting_vel) {
          const CCTK_REAL vxL = velx(p.I);
          const CCTK_REAL vyL = vely(p.I);
          const CCTK_REAL vzL = velz(p.I);
          if ((abs(vxL) < vtol) && (abs(vyL) < vtol) && (abs(vzL) < vtol)) {
            // For star 1 at minus side
            const CCTK_REAL x_local_s1 = p.x - CoM_NS1[0];
            const CCTK_REAL y_local_s1 = p.y - CoM_NS1[1];
            const CCTK_REAL cylrad2_s1 =
                x_local_s1 * x_local_s1 + y_local_s1 * y_local_s1;
            // For star 2 at minus side
            const CCTK_REAL x_local_s2 = p.x - CoM_NS2[0];
            const CCTK_REAL y_local_s2 = p.y - CoM_NS2[1];
            const CCTK_REAL cylrad2_s2 =
                x_local_s2 * x_local_s2 + y_local_s2 * y_local_s2;

            if (cylrad2_s1 < cylrad2_s2) {
              if (cylrad2_s1 < pow(3.0 * radius_NS1, 2.0)) {
                velx(p.I) = vel_NS1[0];
                vely(p.I) = vel_NS1[1];
                velz(p.I) = vel_NS1[2];
              } else {
                CCTK_REAL rfac =
                    pow(3.0 * radius_NS1, 4.0) / pow(cylrad2_s1, 2.0);
                velx(p.I) = rfac * vel_NS1[0];
                vely(p.I) = rfac * vel_NS1[1];
                velz(p.I) = rfac * vel_NS1[2];
              }
            } else {
              if (cylrad2_s2 < pow(3.0 * radius_NS2, 2.0)) {
                velx(p.I) = vel_NS2[0];
                vely(p.I) = vel_NS2[1];
                velz(p.I) = vel_NS2[2];
              } else {
                CCTK_REAL rfac =
                    pow(3.0 * radius_NS2, 4.0) / pow(cylrad2_s2, 2.0);
                velx(p.I) = rfac * vel_NS2[0];
                vely(p.I) = rfac * vel_NS2[1];
                velz(p.I) = rfac * vel_NS2[2];
              }
            }
          } else {
            velx(p.I) = vxL;
            vely(p.I) = vyL;
            velz(p.I) = vzL;
          }
        } // if set coorbiting velocity
      });
  amrex::Gpu::streamSynchronize();
  if (failures.dataValue())
    CCTK_ERROR("Initial beta floor cannot be satisfied within the EOS domain");
}

extern "C" void AsterSeeds_TestBetaFloor(CCTK_ARGUMENTS) {
  std::array<CCTK_REAL, 2> lr{log(1.0e-6), log(1.0e-2)};
  std::array<CCTK_REAL, 2> lt{log(1.0e-3), log(1.0e-1)};
  std::array<CCTK_REAL, 2> ye{0.1, 0.5};
  std::array<CCTK_REAL, 8 * NTABLES> data{};
  for (int k = 0; k < 2; ++k)
    for (int j = 0; j < 2; ++j)
      for (int i = 0; i < 2; ++i) {
        const int n = NTABLES * (i + 2 * (j + 2 * k));
        data[n + eos_3p_tabulated3d::PRESS] = lr[i] + lt[j];
        data[n + eos_3p_tabulated3d::EPS] = lt[j];
        data[n + eos_3p_tabulated3d::S] = lt[j] - lr[i];
      }
  linear_interp_uniform_ND_t<CCTK_REAL, 3, NTABLES> interp(
      data.data(), {2, 2, 2}, lr.data(), lt.data(), ye.data());
  CCTK_REAL shift = 0.05;
  eos_3p_tabulated3d eos;
  eos.interptable = &interp;
  eos.energy_shift = &shift;
  eos.rgrho = {exp(lr[0]), exp(lr[1])};
  eos.rgtemp = {exp(lt[0]), exp(lt[1])};
  eos.rgye = {ye[0], ye[1]};
  eos.rgeps = eos.compute_eps_range_full_table();
  CCTK_REAL rho = 0.0, eps = 0.0, press = 0.0, entropy = 0.0;
  if (!set_beta_floor(&eos, 6.0e-6, 0.02, 0.3, rho, eps, press, entropy) ||
      fabs(rho - 3.0e-4) > 1.0e-14 || fabs(eps + 0.03) > 1.0e-12 ||
      fabs(press - 6.0e-6) > 1.0e-15 ||
      fabs(entropy - log(0.02 / 3.0e-4)) > 1.0e-10)
    CCTK_ERROR("AsterSeeds test: inconsistent beta-floor state");
  const CCTK_REAL rho0 = rho, eps0 = eps, press0 = press, entropy0 = entropy;
  if (set_beta_floor(&eos, 1.0, 0.02, 0.3, rho, eps, press, entropy) ||
      rho != rho0 || eps != eps0 || press != press0 || entropy != entropy0)
    CCTK_ERROR("AsterSeeds test: unattainable beta floor changed the state");
  // A steep pressure interpolant amplifies the density-inversion roundoff.
  std::array<CCTK_REAL, 3> lr_small{log(1.0e-20), log(1.0e-12), log(1.0e-4)};
  std::array<CCTK_REAL, 12 * NTABLES> small_data{};
  for (int k = 0; k < 2; ++k)
    for (int j = 0; j < 2; ++j)
      for (int i = 0; i < 3; ++i) {
        const int n = NTABLES * (i + 3 * (j + 2 * k));
        small_data[n + eos_3p_tabulated3d::PRESS] = 3.0 * lr_small[i] + lt[j];
        small_data[n + eos_3p_tabulated3d::EPS] = lt[j];
      }
  linear_interp_uniform_ND_t<CCTK_REAL, 3, NTABLES> small(
      small_data.data(), {3, 2, 2}, lr_small.data(), lt.data(), ye.data());
  eos.interptable = &small;
  eos.rgrho = {exp(lr_small[0]), exp(lr_small[2])};
  const CCTK_REAL rhoL = 1.624501063509632e-14;
  const CCTK_REAL target = eos.press_from_rho_temp_ye(rhoL, 0.02, 0.3);
  if (!set_beta_floor(&eos, target, 0.02, 0.3, rho, eps, press, entropy) ||
      fabs(rho / rhoL - 1.0) > 1.0e-12 ||
      fabs(press / target - 1.0) > 1.0e-12)
    CCTK_ERROR("AsterSeeds test: attainable beta floor rejected");
  CCTK_INFO("AsterSeeds beta-floor tests passed");
}

} // namespace AsterSeeds
