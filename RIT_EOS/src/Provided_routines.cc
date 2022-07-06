/*******************************************************************************
    Functions provided to other thorns and listed in interface.ccl. As such,
    arguments are passed by value, not by reference, since Cactus automatically
    generates prototypes for them as functions whose arguments must be passed by
    value, not by reference.

    N.B.: Cactus forces input variables to be 'const CCTK_<something>' and
          output variables to be 'CCTK_<something> *restrict' pointers. 
********************************************************************************/

#include <cmath>  // For log () and exp()
#include <cctk.h>
#include <cctk_Parameters.h>
#include <RIT_EOS.hh>

using RIT_EOS_helpers::poly_Gamma_minus_1;
using RIT_EOS_helpers::poly_K_poly_Gamma;
using namespace RIT_EOS_function_callers;



// ************************* POLYTROPIC EOS ROUTINES ***************************

/* Get pressure from the mass density for a polytropic EOS. This routine is
   implemented because polytropic calls may be needed in some situations even
   when using a tabulated EOS.                                                  */
extern "C" void RIT_EOS_press_from_rho_polyEOS_impl(const CCTK_REAL     xrho,
                                                    CCTK_REAL *restrict out_press) {
    DECLARE_CCTK_PARAMETERS;

    // FIXME: is this really needed? Could slow down execution
    if (out_press == nullptr)  CCTK_ERROR("out_press is a null pointer.");

    // Set output
    out_press[0] = poly_K*pow(xrho, poly_Gamma);

    return;
}



/* Get pressure and specific internal energy from the mass density for a
   polytropic EOS. This routine is implemented because polytropic calls may be
   needed in some situations even when using a tabulated EOS.                   */
extern "C" void RIT_EOS_press_eps_from_rho_polyEOS_impl(const CCTK_REAL     xrho,
                                                        CCTK_REAL *restrict out_press,
                                                        CCTK_REAL *restrict out_eps) {
    DECLARE_CCTK_PARAMETERS;

    // FIXME: is this really needed? Could slow down execution
    if (out_press == nullptr)  CCTK_ERROR("out_press is a null pointer.");
    if (out_eps   == nullptr)  CCTK_ERROR("out_eps is a null pointer.");

    // Set output
    out_press[0] = poly_K*pow(xrho, poly_Gamma);
    out_eps[0]   = out_press[0]/(poly_Gamma_minus_1*xrho);

    return;
}



/* Get specific internal energy from mass density and pressure for a polytropic
   EOS. This routine is implemented because polytropic calls may be needed in
   some situations even when using a tabulated EOS.                             */
extern "C" void RIT_EOS_eps_from_rho_press_polyEOS_impl(const CCTK_REAL     xrho,
                                                        const CCTK_REAL     xpress, 
                                                        CCTK_REAL *restrict out_eps) {
    // FIXME: is this really needed? Could slow down execution
    if (out_eps   == nullptr)  CCTK_ERROR("out_eps is a null pointer.");

    // Set output
    out_eps[0] = xpress/(poly_Gamma_minus_1*xrho);

    return;
}



/* Get d(press)/d(rho) from rho for a polytropic EOS. This routine is
   implemented because polytropic calls may be needed in some situations even
   when using a tabulated EOS.                                                  */
extern "C" void RIT_EOS_dpdrhoe_from_rho_polyEOS_impl(const CCTK_REAL     xrho,
                                                      CCTK_REAL *restrict out_dpdrhoe) {
    // FIXME: is this really needed? Could slow down execution
    if (out_dpdrhoe == nullptr)  CCTK_ERROR("out_dpdrhoe is a null pointer.");

    // Set output
    out_dpdrhoe[0] = poly_K_poly_Gamma*pow(xrho, poly_Gamma_minus_1);

    return;
}

// *****************************************************************************





// ******************** ROUTINES WHERE TEMPERATURE IS KNOWN ********************

/* Retrieve pressure and specific internal energy from mass density, temperature
   and electron fraction. If using a polytropic EOS (EOS_key = 1), then just
   passdummy real arguments in place of xtemp and xye.                          */
extern "C" void RIT_EOS_press_eps_from_rho_temp_ye_impl(const CCTK_REAL     xrho,
                                                        const CCTK_REAL     xtemp,
                                                        const CCTK_REAL     xye,
                                                        CCTK_REAL *restrict out_press,
                                                        CCTK_REAL *restrict out_eps) {
    // FIXME: is this really needed? Could slow down execution
    if (out_press == nullptr)  CCTK_ERROR("out_press is a null pointer.");
    if (out_eps   == nullptr)  CCTK_ERROR("out_eps is a null pointer.");

    // Get pressure and specific internal energy
    press_eps_from_rho_temp_ye_caller(xrho, xtemp, xye,
                                      out_press[0], out_eps[0]);

    return;
}

// *****************************************************************************





// *************** ROUTINES WHERE TEMPERATURE IS UNKNOWN ***********************
/* FIXME: think of merging <3Dvar>_temp_from_rho_ye_eps_impl into one single
          routine taking the name of 3Dvar as input                             */
/* Retrieve temperature and pressure from mass density, electron fraction and
   specific internal energy
   N.B.: xeps is NOT const: inside press_temp_from_rho_ye_eps_caller(), the
         energy shift is re-added to it and the log is taken.                   */
extern "C" void RIT_EOS_press_temp_from_rho_ye_eps_impl(const CCTK_REAL     xrho,
                                                        const CCTK_REAL     xye,
                                                        CCTK_REAL           xeps,
                                                        CCTK_REAL *restrict out_press,
                                                        CCTK_REAL *restrict out_temp) {
    // FIXME: is this really needed? Could slow down execution
    if (out_temp  == nullptr)  CCTK_ERROR("out_temp is a null pointer.");
    if (out_press == nullptr)  CCTK_ERROR("out_press is a null pointer.");

    // Get pressure and temperature
    press_temp_from_rho_ye_eps_caller(xrho, xye, xeps,
                                      out_temp[0], out_press[0]);

    return;
}





/* Retrieve temperature and d(press)/d(rho) from mass density, electron fraction
   and specific internal energy
   N.B.: xeps is NOT const: inside dpdrhoe_temp_from_rho_ye_eps_caller(), the
         energy shift is re-added to it and the log is taken.                   */
/* FIXME: is this routine really what I need or should I implement a similar
          one, but where the temperature is known?                              */
extern "C" void RIT_EOS_dpdrhoe_temp_from_rho_ye_eps_impl(const CCTK_REAL     xrho,
                                                          const CCTK_REAL     xye,
                                                          CCTK_REAL           xeps,
                                                          CCTK_REAL *restrict out_dpdrhoe,
                                                          CCTK_REAL *restrict out_temp) {
    // FIXME: is this really needed? Could slow down execution
    if (out_temp    == nullptr)  CCTK_ERROR("out_temp is a null pointer.");
    if (out_dpdrhoe == nullptr)  CCTK_ERROR("out_dpdrhoe is a null pointer.");

    // Get d(press)/d(rho) and temperature
    dpdrhoe_temp_from_rho_ye_eps_caller(xrho, xye, xeps,
                                        out_dpdrhoe[0], out_temp[0]);
    return;
}





/* Retrieve temperature and d(press)/d(eps) from mass density, electron fraction
   and specific internal energy
   N.B.: xeps is NOT const: inside dpderho_temp_from_rho_ye_eps_caller(), the
         energy shift is re-added to it and the log is taken.                   */
/* FIXME: is this routine really what I need or should I implement a similar
          one, but where the temperature is known?                              */
extern "C" void RIT_EOS_dpderho_temp_from_rho_ye_eps_impl(const CCTK_REAL     xrho,
                                                          const CCTK_REAL     xye,
                                                          CCTK_REAL           xeps,
                                                          CCTK_REAL *restrict out_dpderho,
                                                          CCTK_REAL *restrict out_temp) {
    // FIXME: is this really needed? Could slow down execution
    if (out_temp    == nullptr)  CCTK_ERROR("out_temp is a null pointer.");
    if (out_dpderho == nullptr)  CCTK_ERROR("out_dpderho is a null pointer.");

    // Get d(press)/d(eps) and temperature
    dpderho_temp_from_rho_ye_eps_caller(xrho, xye, xeps,
                                        out_dpderho[0], out_temp[0]);
    return;
}





extern "C" void RIT_EOS_dpdrhoe_dpderho_temp_from_rho_ye_eps_impl(const CCTK_REAL     xrho,
                                                                  const CCTK_REAL     xye,
                                                                  CCTK_REAL           xeps,
                                                                  CCTK_REAL *restrict out_dpdrhoe,
                                                                  CCTK_REAL *restrict out_dpderho,
                                                                  CCTK_REAL *restrict out_temp) {
    // FIXME: is this really needed? Could slow down execution
    if (out_temp    == nullptr)  CCTK_ERROR("out_temp is a null pointer.");
    if (out_dpdrhoe == nullptr)  CCTK_ERROR("out_dpdrhoe is a null pointer.");
    if (out_dpderho == nullptr)  CCTK_ERROR("out_dpderho is a null pointer.");

    // Get d(press)/d(rho), d(press)/d(eps) and temperature
    dpdrhoe_dpderho_temp_from_rho_ye_eps_caller(xrho, xye, xeps,
                                                out_dpdrhoe[0], out_dpderho[0],
                                                out_temp[0]);
    return;
}

// *****************************************************************************
