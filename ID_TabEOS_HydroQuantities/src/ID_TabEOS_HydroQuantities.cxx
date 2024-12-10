#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cmath>

#include <AMReX.H>

#include <loop_device.hxx>

#include "ID_TabEOS_HydroQuantities.hxx"

#define SQ(X) ((X) * (X))

namespace ID_TabEOS_HydroQuantities {

using namespace amrex;
using namespace EOSX;
using namespace Loop;

extern "C" void ID_TabEOS_HydroQuantities__initial_Y_e(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ID_TabEOS_HydroQuantities__initial_Y_e;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VInfo(CCTK_THORNSTRING, "Y_e initialization is ENABLED!");

  auto eos_3p_tab3d = global_eos_3p_tab3d;

  // Open the Y_e file, which should countain Y_e(rho) for the EOS table slice
  FILE *Y_e_file = fopen(Y_e_filename, "r");

  // Check if everything is OK with the file
  if ((Y_e_file = fopen(Y_e_filename, "r")) == NULL) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "File \"%s\" does not exist. ABORTING", Y_e_filename);
  } else {
    // Set nrho
    const CCTK_INT nrho = eos_3p_tab3d->interptable->num_points[0];

    // Set interpolation stencil size
    const int interp_stencil_size = 5;

    Ye_reader *id_ye_reader =
        (Ye_reader *)The_Managed_Arena()->alloc(sizeof *id_ye_reader);
    id_ye_reader->init(nrho, Y_e_file, eos_3p_tab3d->interptable->x[0],
                       interp_stencil_size); // eg. logrho array

    // Close the file
    fclose(Y_e_file);

    // Set Y_e
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (rho(p.I) > id_rho_atm_max) {
            // Interpolate Y_e(rho_i) at gridpoint i
            CCTK_REAL Y_eL;
            id_ye_reader->interpolate_1d_quantity_as_function_of_rho(
                interp_stencil_size, nrho, rho(p.I), &Y_eL);
            // Finally, set the Y_e gridfunction
            Ye(p.I) = MIN(MAX(Y_eL, eos_3p_tab3d->interptable->xmin<2>()),
                          eos_3p_tab3d->interptable->xmax<2>());
          } else {
            Ye(p.I) = id_Y_e_atm;
          }
        });

    // if (id_ye_reader->interp_err) {
    //	CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
    //							"1D Interpolator Encountered Error In Y_e
    //Init");
    //}

    The_Managed_Arena()->free(id_ye_reader);
  }
}

// Set initial temperature to be constant everywhere (TODO: add other options)
extern "C" void ID_TabEOS_HydroQuantities__initial_temperature(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ID_TabEOS_HydroQuantities__initial_temperature;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VInfo(CCTK_THORNSTRING, "Temperature initialization is ENABLED!");

  auto eos_3p_tab3d = global_eos_3p_tab3d;

  // Loop over the grid, initializing the temperature
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const CCTK_REAL r_pow_T = atmo_falloff_T ? r_power_T : 0.;
        const CCTK_REAL radius = std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        const CCTK_REAL r_atmo = std::max(r_atmo_min, radius);
        const CCTK_REAL id_T_atm =
            std::max(id_T_atm_max * std::pow(r_atmo / r_atmo_min, r_pow_T),
                     eos_3p_tab3d->interptable->xmin<1>());

        temperature(p.I) = id_T_atm;
      });
}

// Now recompute all HydroQuantities, to ensure consistent initial data
extern "C" void
ID_TabEOS_HydroQuantities__recompute_HydroBase_variables(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ID_TabEOS_HydroQuantities__recompute_HydroBase_variables;
  DECLARE_CCTK_PARAMETERS;

  // Check whether or not we want to initialize the entropy in this thorn
  bool initialize_entropy;
  if (!CCTK_EQUALS(id_entropy_type, "none")) {
    CCTK_VInfo(CCTK_THORNSTRING, "Entropy initialization is ENABLED! Only "
                                 "Table-based ID type is supported ATM.");
    initialize_entropy = true;
  } else {
    CCTK_VInfo(CCTK_THORNSTRING, "Entropy initialization is DISABLED!");
    initialize_entropy = false;
  }

  CCTK_VInfo(CCTK_THORNSTRING, "Recomputing all HydroBase quantities ...");

  auto eos_3p_tab3d = global_eos_3p_tab3d;

  // Loop over the grid, recomputing the HydroBase quantities
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Find Atmospheric Density
        CCTK_REAL rhoL = rho(p.I);
        const CCTK_REAL radius = std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        const CCTK_REAL r_atmo = std::max(r_atmo_min, radius);
        const CCTK_REAL r_pow = atmo_falloff ? r_power : 0.;
        const CCTK_REAL rho_atmL =
            std::max(id_rho_atm_max * std::pow(r_atmo / r_atmo_min, r_pow),
                     eos_3p_tab3d->interptable->xmin<0>());

        if (rhoL > rho_atmL) {
          CCTK_REAL yeL = Ye(p.I);
          CCTK_REAL tempL = temperature(p.I);
          press(p.I) =
              eos_3p_tab3d->press_from_valid_rho_temp_ye(rhoL, tempL, yeL);
          eps(p.I) = eos_3p_tab3d->eps_from_valid_rho_temp_ye(rhoL, tempL, yeL);
          if (initialize_entropy)
            entropy(p.I) =
                eos_3p_tab3d->entropy_from_valid_rho_temp_ye(rhoL, tempL, yeL);
        } else {
          // Reset to atmosphere
          CCTK_REAL tempL = temperature(p.I);

          rho(p.I) = rho_atmL;
          Ye(p.I) = id_Y_e_atm;
          press(p.I) = eos_3p_tab3d->press_from_valid_rho_temp_ye(
              rho_atmL, tempL, id_Y_e_atm);
          eps(p.I) = eos_3p_tab3d->eps_from_valid_rho_temp_ye(rho_atmL, tempL,
                                                              id_Y_e_atm);
          velx(p.I) = 0.0;
          vely(p.I) = 0.0;
          velz(p.I) = 0.0;
          if (initialize_entropy)
            entropy(p.I) = eos_3p_tab3d->entropy_from_valid_rho_temp_ye(
                rhoL, tempL, id_Y_e_atm);
        }
      });
}

} // namespace ID_TabEOS_HydroQuantities
