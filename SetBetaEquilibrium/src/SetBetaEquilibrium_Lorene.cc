// Read in beta-equilibrated initial data in Lorene format
 
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstring>    // For memcpy() // FIXME: see the bottom of this file
#include <algorithm>  // For remove()
#include <vector>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Functions.h>

#include "SetBetaEquilibrium.hh"

using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::to_string;
using std::vector;
using SetBetaEquilibrium_TableToGeom::c2;
using SetBetaEquilibrium_TableToGeom::MeV_J;
using SetBetaEquilibrium_TableToGeom::rho_SItoGeom;





extern "C" void SetBetaEquilibrium_Lorene(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS_CHECKED(SetBetaEquilibrium_Lorene);
    DECLARE_CCTK_PARAMETERS;

    CCTK_INFO("Setting up beta-equilibrated initial data (Lorene EOS format)");


    // ******************** DECLARATION OF VARIABLES ***************************
    CCTK_REAL nbtorho_convfac, log10_nbtorho_convfac_mod;

    // Neutron mass in Lorene: 1.66e-27 kg
    const CCTK_REAL Mn_Mev_over_c2 = 1.66e-27*c2/MeV_J;

    CCTK_INT i, j, k, ijk, ierr;

    const CCTK_INT local_grid_size = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];

    CCTK_REAL logrho;     // Base-10 logarithm of the mass density
    CCTK_REAL T_slice;    // Constant temperature of the EOS slice

    CCTK_REAL nb_tmp;         // Baryon number density (needed to fill logrho_vec)
    CCTK_REAL ye_tmp;         // Electron fraction (needed to fill ye_vec)
    CCTK_REAL totenergy_tmp;  // Total energy (needed to fill totenergy_vec)
    CCTK_REAL press_tmp;      // Pressure (needed to fill press_vec)
    CCTK_REAL dummy;          // Only needed to skip entries when reading the EOS slice

    CCTK_REAL hot_atmo_Y_e;

    vector<CCTK_REAL> logrho_vec, ye_vec, totenergy_vec, eps_vec, press_vec;

    string ye_path, thermo_path;

    string T_slice_string, trash_string;

    ifstream ye_file, thermo_file;

    // *************************************************************************





    // ******************** READ THE ELECTRON FRACTION *************************

    /* Read the values of the beta-equilibrated electron fraction of the EOS
       slice                                                                    */
    ye_path  = ID_path;
    ye_path += ye_filename;
    ye_file.open(ye_path.c_str());

    if (ye_file.is_open()) {
        /* Create a vector containing the values of Y_e. This vector will then
           be passed to the routine 'interp1D'.                                 */
        while (!(ye_file.eof())) {
            ye_file >> ye_tmp;
            ye_vec.push_back(ye_tmp);
        }

        ye_file.close();
    }

    else CCTK_ERROR("Error while opening the file containing the beta-equilibrated electron fraction of the EOS slice. Does the file exist?");

    // *************************************************************************





    /* ***** READ THE TEMPERATURE (CONSTANT), BARYONIC NUMBER DENSITY,
             TOTAL INTERNAL ENERGY AND PRESSURE *********************************/
    thermo_path  = ID_path;
    thermo_path += Lorene_slice_filename;
    thermo_file.open(thermo_path.c_str());

    if (thermo_file.is_open()) {
        /* In thermo_file, the temperature of the EOS slice is written in the
           first line, but has a '#' character in front which must be removed.

        // FIXME: is the following comment correct?
           N.B.: the temperature of the EOS slice is not actually needed for any
                 calculation and the only purpose of reading it is to report its
                 value to the user.                                             */
        getline(thermo_file, T_slice_string);
        T_slice_string.erase(remove(T_slice_string.begin(),
                                    T_slice_string.end(), '#'),
                                T_slice_string.end());

        // Needed when setting up the temperature on the grid
        T_slice = stod(T_slice_string);

        CCTK_VINFO("Temperature will be set to %g MeV.", T_slice);


        /* In thermo_file, thermodynamic variables are written from line number
           10 on; therefore, lines number 2 (temperature has just been read at
           line number 1) to 9 must be skipped.                                 */
        getline(thermo_file, trash_string);  // Line number 2 skipped
        getline(thermo_file, trash_string);  // Line number 3 skipped
        getline(thermo_file, trash_string);  // Line number 4 skipped
        getline(thermo_file, trash_string);  // Line number 5 skipped
        getline(thermo_file, trash_string);  // Line number 6 skipped
        getline(thermo_file, trash_string);  // Line number 7 skipped
        getline(thermo_file, trash_string);  // Line number 8 skipped
        getline(thermo_file, trash_string);  // Line number 9 skipped


        /* The mass density (rho) is related to the baryonic number density (nb)
           and the neutron mass (Mn) via the relation rho = Mn·nb. Now:

           1. Mn_Mev_over_c2*MeV_J/c2 ~= 1.e-27 is the neutron mass in kg;
           2. multiply by 1.e+45 to convert the baryonic number density from
              fm^(-3) to m^(-3) (1 fm := 1.e-15 m^(-3) is the nb unit in
              Compose);
           3. multiply by rho_SItoGeom ~= 1.e-21 to convert the mass density
              from kg/(m^3) to Cactus units.

           Steps 1+2+3 give out a factor ~1.e-27*1.e+45*1.e-21 = 1.e-03 . So,
           multiply the last factor by 1.e+03, take the base-10 log (in this way
           we are taking the log10 of something of order ~1, which is more
           accurate) and subtract 3 (log10(1.e+03) = 3) from the result.

           However, nb is typically in the range [1.e-13, 1.], so ~1.e-07 on
           average. Therefore, multiplying nb by 1.e+07 is convenient in order
           to take the log10 of something of order ~1 (on average); then,
           subtract 7 from the result. Based on the discussion in the previous
           paragraph, 3+7=10 must be subracted from the overall conversion
           factor.

           Given the above, read "1.e+48" below as "1.e+45*1.e+03" and "-10" as
           "-3-7".                                                              */
        nbtorho_convfac           =       Mn_Mev_over_c2*MeV_J*1.e+45*rho_SItoGeom/c2;
        log10_nbtorho_convfac_mod = log10(Mn_Mev_over_c2*MeV_J*1.e+48*rho_SItoGeom/c2) - 10.;

        // FIXME: is it really the specific internal energy?
        /* Create vectors containing the values of the logarithmic mass density
           (built from the baryonic number density -- see comment in
           Spritz_SetBetaEquilibrium_Compose), specific internal energy and
           pressure, respectively. These vectors will then be passed to the
           routine 'Lagrangian_interp'.

           Also, build the specific internal energy (eps) as

              totenergy = (1 + eps)*rho
           => eps = (totenergy/rho) - 1                                         */
        while(!(thermo_file.eof())) {
            thermo_file >> dummy >> nb_tmp >> totenergy_tmp >> press_tmp;
            logrho_vec.push_back(log10(nb_tmp*1.e+07) + log10_nbtorho_convfac_mod);
            totenergy_vec.push_back(totenergy_tmp);
            eps_vec.push_back(totenergy_tmp/(nb_tmp*nbtorho_convfac) - 1.);
            press_vec.push_back(press_tmp);
        }

        thermo_file.close();
    }

    else CCTK_ERROR("Error while opening the file containing the beta-equilibrated baryonic number densities of the EOS slice. Does the file exist?");

    // *************************************************************************





    // *************** PUT THERMODYNAMIC VARIABLES ON THE GRID *****************
    CCTK_INFO("Filling the grid up with temperature, electron fraction and specific internl energy");

    /* Find the value of the atmospheric electron fraction based on the value of
       the atmospheric mass density                                             */
    hot_atmo_Y_e = Lagrange_interp(logrho_vec, ye_vec, log10(rho_abs_min),
                                   Lagrange_interp_order);

    /* If GRHydro is active, set its 'GRHydro_hot_atmo_Y_e' parameter to the
       value of hot_atmo_Y_e                                                    */
    if (CCTK_IsThornActive("GRHydro")) {
        ierr = CCTK_ParameterSet("GRHydro_hot_atmo_Y_e", "GRHydro",
                                 to_string(hot_atmo_Y_e).c_str());
        if (ierr)
            CCTK_ERROR("Error while setting parameter 'GRHydro_hot_atmo_Y_e'");
        else
            CCTK_VINFO("Parameter 'GRHydro_hot_atmo_Y_e' has been set to %s",
                       to_string(hot_atmo_Y_e).c_str());
    }

    else CCTK_WARN(1, "Thorn 'GRHydro' is not active, did you really mean to disable it?");



    /* Set up the temperature in the atmosphere based on the atmospheric mass
       density                                                                  */
    for (k = 0; k < cctk_lsh[2]; ++k)
        for (j = 0; j < cctk_lsh[1]; ++j)
            for (i = 0; i < cctk_lsh[0]; ++i) {
                ijk    = CCTK_GFINDEX3D(cctkGH, i, j, k);
                logrho = log10(fmax(rho[ijk], rho_abs_min));

                /* Temperature is constrained in the table interval
                   [GRHydro_hot_atmo_temp, GRHydro_max_temp].                   */
                temperature[ijk] = T_slice;
                temperature[ijk] = fmax(GRHydro_hot_atmo_temp,
                                        fmin(temperature[ijk], GRHydro_max_temp));


                // Electron fraction setup
                Y_e[ijk] = Lagrange_interp(logrho_vec, ye_vec, logrho,
                                           Lagrange_interp_order);
                /* Electron fraction is constrained in the table interval
                   [GRHydro_Y_e_min, GRHydro_Y_e_max].                          */
                Y_e[ijk] = fmax(GRHydro_Y_e_min, fmin(Y_e[ijk], GRHydro_Y_e_max));


                // Entropy setup // FIXME: is entropy really not available with Lorene?
                /*entropy[ijk] = Lagrange_interp(logrho_vec, entropy_vec, logrho
                                               Lagrange_interp_order);*/

                /* If 'get_EOS_eps_from_beta_eq_ID = "yes" in param.ccl, then
                   get the specific internal energy from beta-equilibrated
                   initial data. Otherwise, the specific internal energy will be
                   set during Prim2ConInitial (TODO: check the last statement). */
                /*if(get_EOS_eps_from_beta_eq_ID)
                    eps[ijk] = Lagrange_interp(logrho_vec, eps_vec, logrho,
                                               Lagrange_interp_order);*/

                // TODO: check the following statement
                /* N.B.: during Prim2ConInitial, temperature and entropy will be
                         overwritten using rho, eps, Y_e; pressure will be then
                         computed using the EOS.                                */


                eps[ijk] = Lagrange_interp(logrho_vec, eps_vec, logrho,
                                           Lagrange_interp_order);

                // TODO: check the following statement
                /* N.B.: during Prim2ConInitial, temperature and entropy will be
                         overwritten using rho, eps, Y_e; pressure will be then
                         computed using the EOS.                                */
            }



    /* FIXME: is the following really needed? Can't I just use
              Carpet::init_fill_timelevels = "yes"?                             */
    /* Fill up timelevels for temperature, electron fraction and specific
       internal energy                                                          */
    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::temperature") > 1) {
        memcpy(temperature_p, temperature, local_grid_size*sizeof(CCTK_REAL));
        if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::temperature") > 2)
            memcpy(temperature_p_p, temperature, local_grid_size*sizeof(CCTK_REAL));
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Y_e") > 1) {
        memcpy(Y_e_p, Y_e, local_grid_size*sizeof(CCTK_REAL));
        if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Y_e") > 2)
            memcpy(Y_e_p_p, Y_e, local_grid_size*sizeof(CCTK_REAL));
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::eps") > 1) {
        memcpy(eps_p, eps, local_grid_size*sizeof(CCTK_REAL));
        if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::eps") > 2)
            memcpy(eps_p_p, eps, local_grid_size*sizeof(CCTK_REAL));
    }

    CCTK_VINFO("Beta equilibrium set on on reflevel %d",
               GetRefinementLevel(cctkGH));

    return;
}
