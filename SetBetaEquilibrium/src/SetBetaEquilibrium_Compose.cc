// Read in beta-equilibrated initial data in Compose format

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstring>  // For memcpy() // FIXME: see the bottom of this file
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





extern "C" void SetBetaEquilibrium_Compose(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS_CHECKED(SetBetaEquilibrium_Compose);
    DECLARE_CCTK_PARAMETERS;

    CCTK_INFO("Setting up beta-equilibrated initial data (Compose EOS format)");


    // ******************** DECLARATION OF VARIABLES ***************************
    CCTK_REAL Mn_Mev_over_c2, log10_nbtorho_convfac_mod;
    CCTK_INT i, j, k, ijk, ierr;

    const CCTK_INT local_grid_size = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];

    CCTK_REAL logrho;   // Base-10 logarithm of the mass density
    CCTK_REAL T_slice;  // Constant temperature of the EOS slice if the latter is taken at constant temperature

    CCTK_REAL nb_tmp;         // Baryon number density (needed to fill logrho_vec)
    CCTK_REAL T_const_s_tmp;  // Temperature of the EOS slice if the latter is taken at constant entropy (needed to fill T_const_s_vec)
    CCTK_REAL ye_tmp;         // Electron fraction (needed to fill ye_vec)
    CCTK_REAL eps_tmp;        // Specific internal energy (needed to fill eps_vec)
    CCTK_REAL s_tmp;          // Entropy per baryon (needed to fill s_vec)
    CCTK_REAL dummy;          // Only needed to skip entries when reading the EOS slice

    CCTK_REAL hot_atmo_Y_e;

    vector<CCTK_REAL> logrho_vec, T_const_s_vec, ye_vec, eps_vec, s_vec;

    bool entropy_slice;

    string T_path, T_const_s_path, ye_path, thermo_path, nb_path;
    string T_slice_string, trash_string;

    ifstream T_file, T_const_s_file, ye_file, thermo_file, nb_file;

    // *************************************************************************





    /* *************** READ THE BARYONIC NUMBER DENSITY
                       AND BUILD THE (LOGARITHMIC) MASS DENSITY *****************/
    nb_path  = ID_path;
    nb_path += ".nb";
    nb_file.open(nb_path.c_str());

    /* Read the values of the beta-equilibrated baryonic number density of the
       EOS slice                                                                */
    if (nb_file.is_open()) {
        /* In this file, the first value of the beta-equilibrated baryonic
           number density of the EOS slice is written at line number 3.
           Therefore, the first two lines in this file must be skipped.         */
        getline(nb_file, trash_string);  // Line number 1 skipped
        getline(nb_file, trash_string);  // Line number 2 skipped

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
        log10_nbtorho_convfac_mod = log10(Mn_Mev_over_c2*MeV_J*1.e+48*rho_SItoGeom/c2) - 10.;

        while (!(nb_file.eof())) {
            nb_file >> nb_tmp;
            logrho_vec.push_back(log10(nb_tmp*1.e+07) + log10_nbtorho_convfac_mod);
        }

        nb_file.close();
    }

    else CCTK_ERROR("Error while opening the file containing the beta-equilibrated baryonic number densities of the EOS slice. Does the file exist?");

    // *************************************************************************





    // ************************* READ THE TEMPERATURE **************************

    /* Read the values of the temperature of the EOS slice if the latter is
       taken at constant entropy; in this case, a file with extension ".t_s"
       exists which contains such temperatures. Whether the slice is taken at
       constant entropy or not is a piece of information which will be used when
       setting up the temperature on the grid, so that's the reason for defining
       the boolean variable 'entropy_slice'.                                    */
    T_const_s_path  = ID_path;
    T_const_s_path += ".t_s";
    T_const_s_file.open(T_const_s_path.c_str());
    entropy_slice = T_const_s_file.is_open();

    // Read the ".t_s" file if the EOS slice is taken at constant entropy
    if (entropy_slice) {
        CCTK_INFO ("Temperature will be set at constant entropy.");

        while (!(T_const_s_file.eof())) {
            T_const_s_file >> T_const_s_tmp;
            T_const_s_vec.push_back(T_const_s_tmp);
        }

        T_const_s_file.close();

        // Consistency check
        if (T_const_s_vec.size() != logrho_vec.size())
            CCTK_VERROR("In the EOS slice, the size of the table of the temperatures at constant entropy (%d) is different from the size of the baryonic number density table (%d).",
                        T_const_s_vec.size(), logrho_vec.size());
    }



    /* If, instead, the EOS slice is taken at constant temperature, then simply
       read the temperature of the slice from the ".t" file.                    */
    else {
        T_path  = ID_path;
        T_path += ".t";
        T_file.open(T_path.c_str());

        if (T_file.is_open()) {
            /* In T_file, the temperature of the EOS slice is written at line
               number 3. Therefore, the first two lines in this file must be
               skipped.                                                         */
            getline(T_file, T_slice_string);  // Line 1 skipped
            getline(T_file, T_slice_string);  // Line 2 skipped
            getline(T_file, T_slice_string);  // Line 3 contains the temperature

            T_slice = stoi(T_slice_string);
            T_file.close();
            CCTK_VINFO("Temperature will be set to %f MeV.", T_slice);
        }

        else CCTK_ERROR("Error while opening the file containing the temperature of the EOS slice. Does the file exist?");
    }

    // *************************************************************************





    // ******************** READ THE ELECTRON FRACTION *************************

    /* Read the values of the beta-equilibrated electron fraction of the EOS
       slice                                                                    */
    ye_path  = ID_path;
    ye_path += ".yqb";
    ye_file.open(ye_path.c_str());

    if (ye_file.is_open()) {
        /* Create a vector containing the values of Y_e. This vector will then
           be passed to the routine 'Lagrangian_interp'.                        */
        while (!(ye_file.eof())) {
            ye_file >> ye_tmp;
            ye_vec.push_back(ye_tmp);
        }

        ye_file.close();
    }

    else CCTK_ERROR("Error while opening the file containing the beta-equilibrated electron fraction of the EOS slice. Does the file exist?");

    // *************************************************************************





    // ********** READ THE ENTROPY PER BARYON AND ENERGY PER BARYON ************

    /* Read the values of the beta-equilibrated entropy per baryon and specific
       internal energy of the EOS slice. Also, from the same file, read the
       neutron mass in MeV/(c^2): this will be used to build the mass density
       starting from the baryonic number density. For this reason, this file
       must read *BEFORE* reading the ".nb" file.                               */
    thermo_path  = ID_path;
    thermo_path += ".thermo";
    thermo_file.open(thermo_path.c_str());

    if (thermo_file.is_open()) {
        /* In thermo_file, the neutron mass in MeV is the first number in the
           first line.                                                          */
        thermo_file >> Mn_Mev_over_c2;
        getline(thermo_file, trash_string);  // Ignore the remaining part of the first line

        /* Create vectors containing the values of entropy per baryon and
           specific internal energy, respectively. These vectors will then be
           passed to the routine 'Lagrangian_interp'.                           */
        while (!(thermo_file.eof())) {
          thermo_file >> dummy >> dummy >> dummy >> dummy >> s_tmp >> dummy
                      >> dummy >> dummy >> dummy >> eps_tmp;
          s_vec.push_back(s_tmp);
          eps_vec.push_back(eps_tmp);
        }

        thermo_file.close();
    }

    else CCTK_ERROR("Error while opening the file containing the beta-equilibrated thermodynamic quantities of the EOS slice. Does the file exist?");

    // *************************************************************************





    // *************** PUT THERMODYNAMIC VARIABLES ON THE GRID *****************
    CCTK_INFO("Filling the grid up with temperature, electron fraction, entropy and specific internal energy (Compose format)");

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

                // Temperature setup
                entropy_slice ?
                    (temperature[ijk] = Lagrange_interp(logrho_vec,
                                                        T_const_s_vec,
                                                        logrho,
                                                        Lagrange_interp_order)) :
                    (temperature[ijk] = T_slice);

                /* Temperature is constrained in the table interval
                   [GRHydro_hot_atmo_temp, GRHydro_max_temp].                   */
                temperature[ijk] = fmax(GRHydro_hot_atmo_temp,
                                        fmin(temperature[ijk], GRHydro_max_temp));


                // Electron fraction setup
                Y_e[ijk] = Lagrange_interp(logrho_vec, ye_vec, logrho,
                                           Lagrange_interp_order);
                /* Electron fraction is constrained in the table interval
                   [GRHydro_Y_e_min, GRHydro_Y_e_max].                          */
                Y_e[ijk] = fmax(GRHydro_Y_e_min, fmin(Y_e[ijk], GRHydro_Y_e_max));


                // Entropy setup
                entropy[ijk] = Lagrange_interp(logrho_vec, s_vec, logrho,
                                               Lagrange_interp_order);

                /* If 'get_EOS_eps_from_beta_eq_ID = "yes" in param.ccl, then
                   get the specific internal energy from beta-equilibrated
                   initial data. Otherwise, the specific internal energy will be
                   set during Prim2ConInitial (TODO: check the last statement). */
                if(get_EOS_eps_from_beta_eq_ID)
                    eps[ijk] = Lagrange_interp(logrho_vec, eps_vec,
                                               logrho, Lagrange_interp_order);

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
