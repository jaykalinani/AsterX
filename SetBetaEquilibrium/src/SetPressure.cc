// Set up initial pressure based on beta-equilibrated initial data

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
//#include <cctk_Functions.h>



extern "C" void SetPressure(CCTK_ARGUMENTS) {
    /* TODO: remember to add READ/WRITES statements in schedule.ccl. Only those
             variables which are identified by READ/WRITES statements in
             schedule.ccl will be actually declared.                            */
    DECLARE_CCTK_ARGUMENTS_CHECKED(SetPressure);
    DECLARE_CCTK_PARAMETERS;

    CCTK_INFO("Setting up the initial pressure using the EOS");


    // Declaration of variables
          CCTK_INT  i, j, k, i3D;
    const CCTK_INT  local_grid_size = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];
          CCTK_INT  anyerr;
          CCTK_INT  nn = 1, xkeytemp, keyerr[1];
          CCTK_REAL xrho[1], xpress[1], xeps[1], xtemp[1], xye[1];

    // FIXME FIXME FIXME: USE RIT_EOS IN PLACE OF EOS_OMNI!!! FIXME FIXME FIXME
    CCTK_INT GRHydro_eos_handle[1];
    GRHydro_eos_handle[0] = 4;

    // TODO: comment
    for (k = 0; k < cctk_lsh[2]; ++k)
        for (j = 0; j < cctk_lsh[1]; ++j)
            for (i = 0; i < cctk_lsh[0]; ++i) {
                i3D = CCTK_GFINDEX3D(cctkGH, i, j, k);
                //xkeytemp  = 0;  // FIXME: set to 1 below
                anyerr    = 0;
                keyerr[0] = 0;
                xrho[0]   = rho[i3D];
                xeps[0]   = eps[i3D];
                xtemp[0]  = temperature[i3D];
                xye[0]    = Y_e[i3D];
                xkeytemp  = 1; //temperature has already been set, there is no need of root-finding it again.
                // keytemp=1 means that press is found interpolating in rho,Ye,T table.
                EOS_Omni_press(*GRHydro_eos_handle,/* *Spritz_eos_handle,*/ xkeytemp, GRHydro_eos_rf_prec, nn,
                               xrho, xeps, xtemp, xye, xpress, keyerr, &anyerr);
                press[i3D] = xpress[0];
        }

    // if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::press") > 1) {
    //   memcpy(press_p, press, local_grid_size*sizeof(CCTK_REAL));
    //   if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::press") > 2) {
    //     memcpy(press_p_p, press, local_grid_size*sizeof(CCTK_REAL));
    //   }
    // }
    CCTK_VINFO ("Pressure set on on reflevel %d.", GetRefinementLevel(cctkGH));

    return;
}



