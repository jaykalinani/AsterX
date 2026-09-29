#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cmath>

#include <AMReX.H>

#include <loop_device.hxx>

#include "ID_TabEOS_HydroQuantities.hxx"
#include "atmo.hxx"

#define SQ(X) ((X) * (X))

namespace ID_TabEOS_HydroQuantities {

using namespace amrex;
using namespace Con2PrimFactory;
using namespace EOSX;
using namespace Loop;

enum class TS_ID_t { Temperature, Entropy };

extern "C" void ID_TabEOS_HydroQuantities_initial_Y_e(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ID_TabEOS_HydroQuantities_initial_Y_e;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VInfo(CCTK_THORNSTRING, "Y_e initialization is ENABLED!");

  auto eos_1p_poly = global_eos_1p_poly;
  auto eos_3p_tab3d = global_eos_3p_tab3d;

  // Open the Y_e file, which should countain Y_e(rho) for the EOS table slice
  FILE *const Y_e_file = fopen(Y_e_filename, "r");

  // Check if everything is OK with the file
  if (Y_e_file == NULL) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "File \"%s\" does not exist. ABORTING", Y_e_filename);
  } else {
    // Set nrho
    const CCTK_INT nrho = eos_3p_tab3d->interptable->num_points[0];

    // Set interpolation stencil size
    const int interp_stencil_size = 5;

    Ye_reader id_ye_reader;
    id_ye_reader.init(nrho, Y_e_file, eos_3p_tab3d->interptable->x[0],
                      interp_stencil_size); // eg. logrho array

    // Close the file
    fclose(Y_e_file);

    // Set Y_e
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const CCTK_REAL radial_distance =
              sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
          const auto atmo = make_atmosphere(
              eos_1p_poly, eos_3p_tab3d, radial_distance, rho_abs_min,
              p_atmo, t_atmo, Ye_atmo, r_atmo, n_rho_atmo, n_press_atmo,
              n_temp_atmo, atmo_tol, true, false);

          if (rho(p.I) > atmo.rho_cell_reset_cut) {
            // Interpolate Y_e(rho_i) at gridpoint i
            const CCTK_REAL Y_eL =
                id_ye_reader.interpolate_1d_quantity_as_function_of_rho(
                    interp_stencil_size, nrho, rho(p.I));
            // Finally, set the Y_e gridfunction
            Ye(p.I) = limit_to_range(Y_eL, eos_3p_tab3d->rgye);
          } else {
            Ye(p.I) = atmo.ye_atmo;
          }
        });

    // loop_all_device is asynchronous. Keep the reader's managed arrays alive
    // until the device has finished using the by-value reader copy.
    Gpu::synchronize();
    id_ye_reader.release();
  }
}

// Set the initial temperature or entropy profile.
extern "C" void ID_TabEOS_HydroQuantities_initial_temp_ent(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ID_TabEOS_HydroQuantities_initial_temp_ent;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VInfo(CCTK_THORNSTRING,
             "Temperature and entropy initialization is ENABLED!");

  auto eos_1p_poly = global_eos_1p_poly;
  auto eos_3p_tab3d = global_eos_3p_tab3d;

  TS_ID_t ts_ID;

  if (CCTK_EQUALS(id_temp_ent_type, "constant temperature")) {
    ts_ID = TS_ID_t::Temperature;
  } else if (CCTK_EQUALS(id_temp_ent_type, "constant entropy")) {
    ts_ID = TS_ID_t::Entropy;
  } else {
    CCTK_ERROR("Unknown value for parameter \"id_temp_ent_type\"");
  }

  // Loop over the grid, initializing the temperature
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const CCTK_REAL radial_distance =
            std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        const auto atmo = make_atmosphere(
            eos_1p_poly, eos_3p_tab3d, radial_distance, rho_abs_min,
            p_atmo, t_atmo, Ye_atmo, r_atmo, n_rho_atmo, n_press_atmo,
            n_temp_atmo, atmo_tol, true, false);

        CCTK_REAL rhoL = rho(p.I);
        CCTK_REAL yeL = Ye(p.I);

        switch (ts_ID) {
        case TS_ID_t::Temperature: {
          const auto state = state_from_rho_temp_ye(
              eos_3p_tab3d, rhoL, atmo.temp_atmo, yeL);
          temperature(p.I) = state.temperature;
          entropy(p.I) = state.kappa;
          break;
        }
        case TS_ID_t::Entropy: {
          if (rhoL > atmo.rho_cell_reset_cut) {
            CCTK_REAL ent_val = id_entropy;
            CCTK_REAL temp_val = eos_3p_tab3d->temp_from_rho_entropy_ye(
                rhoL, ent_val, yeL);
            const auto state =
                state_from_rho_temp_ye(eos_3p_tab3d, rhoL, temp_val, yeL);
            entropy(p.I) = state.kappa;
            temperature(p.I) = state.temperature;
          } else {
            temperature(p.I) = atmo.temp_atmo;
            entropy(p.I) = atmo.entropy_atmo;
          }
          break;
        }
        default:
          assert(0);
        }

        if (temperature(p.I) < 0.0) {
          printf("Negative input for temperature at I=(%d,%d,%d) "
                 "(x=%.5e y=%.5e z=%.5e): temp=%.5e\n",
                 p.I[0], p.I[1], p.I[2], p.x, p.y, p.z, temperature(p.I));
        }
      });
}

// Now recompute all HydroQuantities, to ensure consistent initial data
extern "C" void
ID_TabEOS_HydroQuantities_recompute_HydroBase_variables(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ID_TabEOS_HydroQuantities_recompute_HydroBase_variables;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VInfo(CCTK_THORNSTRING, "Recomputing all HydroBase quantities ...");

  auto eos_1p_poly = global_eos_1p_poly;
  auto eos_3p_tab3d = global_eos_3p_tab3d;

  // Loop over the grid, recomputing the HydroBase quantities
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const CCTK_REAL radial_distance =
            std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        const auto atmo = make_atmosphere(
            eos_1p_poly, eos_3p_tab3d, radial_distance, rho_abs_min,
            p_atmo, t_atmo, Ye_atmo, r_atmo, n_rho_atmo, n_press_atmo,
            n_temp_atmo, atmo_tol, true, false);

        const bool reset_to_atmosphere =
            rho(p.I) <= atmo.rho_cell_reset_cut;
        const auto state =
            reset_to_atmosphere
                ? state_from_rho_temp_ye(eos_3p_tab3d, atmo.rho_atmo,
                                         atmo.temp_atmo, atmo.ye_atmo)
                : state_from_rho_temp_ye(eos_3p_tab3d, rho(p.I),
                                         temperature(p.I), Ye(p.I));

        rho(p.I) = state.rho;
        eps(p.I) = state.eps;
        Ye(p.I) = state.Ye;
        press(p.I) = state.press;
        temperature(p.I) = state.temperature;
        entropy(p.I) = state.kappa;

        if (reset_to_atmosphere) {
          // Initial-data atmosphere is static by construction.
          velx(p.I) = 0.0;
          vely(p.I) = 0.0;
          velz(p.I) = 0.0;
        }
      });
}

} // namespace ID_TabEOS_HydroQuantities
