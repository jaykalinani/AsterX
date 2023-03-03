#include <cmath>  // For pow()
#include <cctk.h>
#include <cctk_Parameters.h>
#include "RIT_EOS.hh"

using RIT_EOS_helpers::poly_Gamma_minus_1;
using RIT_EOS_helpers::poly_Gamma_times_poly_Gamma_minus_1;
using RIT_EOS_helpers::poly_K_poly_Gamma;



/* Retrieve pressure and specific internal energy from the mass density
   ***** DUMMY ARGUMENTS *****
   *        xtemp, xye       *
   ***************************                                                  */
void press_eps_from_rho_temp_ye_polyEOS(const CCTK_REAL &xrho, 
                                        const CCTK_REAL &xtemp,
                                        const CCTK_REAL &xye,
                                              CCTK_REAL &out_press,
                                              CCTK_REAL &out_eps) {
    DECLARE_CCTK_PARAMETERS;

    /* Sanity check
       N.B.: not worrying about xtemp, since it's just a dummy input            */
    if (xrho <= 0.)  CCTK_ERROR("xrho <= 0: not a good sign");

    // Set output
    out_press = poly_K*pow(xrho, poly_Gamma);
    out_eps   = out_press/(poly_Gamma_minus_1*xrho);

    return;
}





/* Retrieve temperature and pressure from mass density, electron fraction and
   specific internal energy.
   ***** DUMMY ARGUMENTS *****
   *    xtemp, xye, xeps     *
   ***************************
   N.B.: xeps is NOT const because it is not const in the "nucEOS" version of
         this routine and both routines must be assigned to the same function
         pointer ('function caller').
   N.B.: this function returns a temperature value of zero, which is *NONSENSE*
         (it's a polytropic EOS call) and must NOT be used. This is ONLY done
         because the "twin" routine press_temp_from_rho_ye_eps_nucEOS() exists
         (see NucEOS_routines.cc) and the two routines must return the same
         type and number of arguments (and must accept the same number and type
         of arguments, as well) in order for the function pointer pointing to
         either of the two to be well defined.
   N.B.: xeps is not used and this routine only gives the pressure in the end
         -> FIXME: ***** IS THIS THING OF HAVING A FUNCTION CALLER   *****
                   ***** (and so the 'polyEOS' and 'nucEOS' version  *****
                   ***** of each routine) REALLY NECESSARY?          *****
                   ***** In the end, it's just needed to avoid ONE   *****
                   ***** if condition (like if-polyEOS...//else...), *****
                   ***** but there are PLENTY of them in the acual   *****
                   ***** interpolators etc. and it makes the code a  *****
                   ***** mess.                                       *****      */
void press_temp_from_rho_ye_eps_polyEOS(const CCTK_REAL &xrho,
                                        const CCTK_REAL &xye,
                                              CCTK_REAL &xeps,
                                              CCTK_REAL &out_press,
                                              CCTK_REAL &out_temp) {
    DECLARE_CCTK_PARAMETERS;

    /* Sanity check
       N.B.: not worrying about xeps, since it's just a dummy input             */
    if (xrho <= 0.)  CCTK_ERROR("xrho <= 0: not a good sign");

    // Set output
    out_press = poly_K*pow(xrho, poly_Gamma);
    out_temp  = 0.;  // *JUNK* temp value

    return;
}





/* Retrieve temperature and d(press)/d(rho) from mass density, electron fraction
   and specific internal energy.
   ***** DUMMY ARGUMENTS *****
   *           xye           *
   ***************************
   N.B.: xeps is NOT const because it is not const in the "nucEOS" version of
         this routine and both routines must be assigned to the same function
         pointer ('function caller').
   N.B.: this function returns a temperature value of zero, which is *NONSENSE*
         (it's  polytropic EOS call) and must NOT be used. This is ONLY done
         because the "twin" routine dpdrhoe_temp_from_rho_ye_eps_nucEOS() exists
         (see NucEOS_routines.cc) and the two routines must return the same type
         and number of arguments (and must accept the same number and type of
         arguments, as well) in order for the function pointer pointing to
         either of the two to be well defined.                                  */
void dpdrhoe_temp_from_rho_ye_eps_polyEOS(const CCTK_REAL &xrho, 
                                          const CCTK_REAL &xye,
                                                CCTK_REAL &xeps,
                                                CCTK_REAL &out_dpdrhoe,
                                                CCTK_REAL &out_temp) {
    // Sanity checks
    if (xrho <= 0.)  CCTK_ERROR("xrho <= 0: not a good sign");
    if (xeps <= 0.)  CCTK_ERROR("xeps <= 0: not a good sign");

    // Set output
    out_dpdrhoe = poly_Gamma_times_poly_Gamma_minus_1*xeps;
    out_temp    = 0.;  // *JUNK* temp value

    return;
}





/* Retrieve temperature and d(press)/d(eps) from mass density, electron fraction
   and specific internal energy.
   ***** DUMMY ARGUMENTS *****
   *        xye, xeps        *
   ***************************
   N.B.: xeps is NOT const because it is not const in the "nucEOS" version of
         this routine and both routines must be assigned to the same function
         pointer ('function caller').
   N.B.: this function returns a temperature value of zero, which is *NONSENSE*
         (it's  polytropic EOS call) and must NOT be used. This is ONLY done
         because the "twin" routine dpderho_temp_from_rho_ye_eps_nucEOS() exists
         (see NucEOS_routines.cc) and the two routines must return the same type
         and number of arguments (and must accept the same number and type of
         arguments, as well) in order for the function pointer pointing to
         either of the two to be well defined.                                  */
void dpderho_temp_from_rho_ye_eps_polyEOS(const CCTK_REAL &xrho, 
                                          const CCTK_REAL &xye,
                                                CCTK_REAL &xeps,
                                                CCTK_REAL &out_dpderho,
                                                CCTK_REAL &out_temp) {
    /* Sanity check
       N.B.: not worrying about xeps, since it's just a dummy input             */
    if (xrho <= 0.)  CCTK_ERROR("xrho <= 0: not a good sign");

    // Set output
    out_dpderho = poly_Gamma_minus_1*xrho;
    out_temp    = 0.;  // *JUNK* temp value

    return;
}





/* Retrieve temperature, d(press)/d(rho) and d(press)/d(eps) from mass density,
   electron fraction and specific internal energy.
   ***** DUMMY ARGUMENTS *****
   *        xye, xeps        *
   ***************************
   N.B.: xeps is NOT const because it is not const in the "nucEOS" version of
         this routine and both routines must be assigned to the same function
         pointer ('function caller').
   N.B.: this function returns a temperature value of zero, which is *NONSENSE*
         (it's  polytropic EOS call) and must NOT be used. This is ONLY done
         because the "twin" routine
         dpdrhoe_dpderho_temp_from_rho_ye_eps_nucEOS() exists (see
         NucEOS_routines.cc) and the two routines must return the same type and
         number of arguments (and must accept the same number and type of
         arguments, as well) in order for the function pointer pointing to
         either of the two to be well defined.                                  */
void dpdrhoe_dpderho_temp_from_rho_ye_eps_polyEOS(const CCTK_REAL &xrho, 
                                                  const CCTK_REAL &xye,
                                                        CCTK_REAL &xeps,
                                                        CCTK_REAL &out_dpdrhoe,
                                                        CCTK_REAL &out_dpderho,
                                                        CCTK_REAL &out_temp) {
    /* Sanity check
       N.B.: not worrying about xeps, since it's just a dummy input             */
    if (xrho <= 0.)  CCTK_ERROR("xrho <= 0: not a good sign");

    // Set output
    out_dpdrhoe = poly_Gamma_times_poly_Gamma_minus_1*xeps;
    out_dpderho = poly_Gamma_minus_1*xrho;
    out_temp    = 0.;  // *JUNK* temp value

    return;
}
