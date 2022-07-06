#include <cmath>  // For exp()
#include <cctk.h>
#include <vector>
#include "RIT_EOS.hh"

using std::vector;
using namespace RIT_EOS_keys;



/* Retrieve pressure and specific internal energy from mass density, temperature
   and electron fraction                                                        */
void press_eps_from_rho_temp_ye_nucEOS(const CCTK_REAL &xrho, 
                                       const CCTK_REAL &xtemp,
                                       const CCTK_REAL &xye,
                                             CCTK_REAL &out_press,
                                             CCTK_REAL &out_eps) {
    // Sanity checks
    if (xrho  <= 0.)  CCTK_ERROR("xrho <= 0: log(xrho) is undefined and negative mass densities are anyway strange.");
    if (xtemp <= 0.)  CCTK_ERROR("xtemp <= 0: log(xtemp) is undefined and negative temperatures are anyway strange.");

    // Convert input to logarithms
    CCTK_REAL xlogrho  = log(xrho);
    CCTK_REAL xlogtemp = log(xtemp);

    // Choose the 3D quantities to interpolate
    const vector<key> which_vars3D{press_key, eps_key};

    /* Get pressure and specific internal energy from xlogrho, xlogtemp, xye.
       out_vars[0] = press
       out_vars[1] = eps                                                        */
    CCTK_REAL out_vars[2];
    vars3D_from_logrho_logtemp_ye(xlogrho, xlogtemp, xye,
                                  which_vars3D, out_vars);

    // Set output
    out_press = out_vars[0];
    out_eps   = out_vars[1];

    return;
}





/* Retrieve temperature and pressure from mass density, electron fraction and
   specific internal energy
   N.B.: xeps is NOT const: inside logtemp_from_logrho_ye_var3D(), the energy
         shift is re-added to it and the log is taken.                          */
void press_temp_from_rho_ye_eps_nucEOS(const CCTK_REAL &xrho,
                                       const CCTK_REAL &xye,
                                             CCTK_REAL &xeps,
                                             CCTK_REAL &out_press,
                                             CCTK_REAL &out_temp) {
    // Sanity checks
    if (xrho <= 0.)  CCTK_ERROR("xrho <= 0: log(xrho) is undefined (and anyway, not a good sign).");
    if (xeps <= 0.)  CCTK_ERROR("xeps <= 0: log(xeps) is undefined (and anyway, not a good sign).");

    // Convert input to logarithm
    CCTK_REAL xlogrho = log(xrho);

    // Get temperature from xlogrho, xye, xeps
    CCTK_REAL out_logtemp = logtemp_from_logrho_ye_var3D(xlogrho, xye,
                                                         xeps, eps_key);
    out_temp = exp(out_logtemp);

    /* Get pressure from xlogrho, out_logtemp, xye
       out_vars[0] = press                                                      */
    const vector<key> which_vars3D{press_key};
    vars3D_from_logrho_logtemp_ye(xlogrho, out_logtemp, xye,
                                  which_vars3D, &out_press);
    return;
}





/* Retrieve temperature and d(press)/d(rho) from mass density, electron fraction
   and specific internal energy
   N.B.: xeps is NOT const: inside logtemp_from_logrho_ye_var3D(), the energy
         shift is re-added to it and the log is taken.                          */
void dpdrhoe_temp_from_rho_ye_eps_nucEOS(const CCTK_REAL &xrho, 
                                         const CCTK_REAL &xye,
                                               CCTK_REAL &xeps,
                                               CCTK_REAL &out_dpdrhoe,
                                               CCTK_REAL &out_temp) {
    // Sanity checks
    if (xrho <= 0.)  CCTK_ERROR("xrho <= 0: log(xrho) is undefined (and anyway, not a good sign).");
    if (xeps <= 0.)  CCTK_ERROR("xeps <= 0: log(xeps) is undefined (and anyway, not a good sign).");

    // Convert input to logarithms
    CCTK_REAL xlogrho = log(xrho);

    // Get temperature from xlogrho, xye, xeps
    CCTK_REAL out_logtemp = logtemp_from_logrho_ye_var3D(xlogrho, xye,
                                                         xeps, eps_key);
    out_temp = exp(out_logtemp);

    // Get d(press)/d(rho) from xlogrho, out_logtemp, xye
    const vector<key> which_vars3D{dpdrhoe_key};
    vars3D_from_logrho_logtemp_ye(xlogrho, out_logtemp, xye,
                                  which_vars3D, &out_dpdrhoe);
    return;
}





/* Retrieve temperature and d(press)/d(eps) from mass density, electron fraction
   and specific internal energy
   N.B.: xeps is NOT const: inside logtemp_from_logrho_ye_var3D(), the energy
         shift is re-added to it and the log is taken.                          */
void dpderho_temp_from_rho_ye_eps_nucEOS(const CCTK_REAL &xrho, 
                                         const CCTK_REAL &xye,
                                               CCTK_REAL &xeps,
                                               CCTK_REAL &out_dpderho,
                                               CCTK_REAL &out_temp) {
    // Sanity checks
    if (xrho <= 0.)  CCTK_ERROR("xrho <= 0: log(xrho) is undefined (and anyway, not a good sign).");
    if (xeps <= 0.)  CCTK_ERROR("xeps <= 0: log(xeps) is undefined (and anyway, not a good sign).");

    // Convert input to logarithms
    CCTK_REAL xlogrho = log(xrho);

    // Get temperature from xlogrho, xye, xeps
    CCTK_REAL out_logtemp = logtemp_from_logrho_ye_var3D(xlogrho, xye,
                                                         xeps, eps_key);
    out_temp = exp(out_logtemp);

    // Get d(press)/d(rho) from xlogrho, out_logtemp, xye
    const vector<key> which_vars3D{dpderho_key};
    vars3D_from_logrho_logtemp_ye(xlogrho, out_logtemp, xye,
                                  which_vars3D, &out_dpderho);
    return;
}





/* Retrieve temperature, d(press)/d(rho) and d(press)/d(eps) from mass density,
   electron fraction and specific internal energy
   N.B.: xeps is NOT const: inside logtemp_from_logrho_ye_var3D(), the energy
         shift is re-added to it and the log is taken.                          */
void dpdrhoe_dpderho_temp_from_rho_ye_eps_nucEOS(const CCTK_REAL &xrho, 
                                                 const CCTK_REAL &xye,
                                                       CCTK_REAL &xeps,
                                                       CCTK_REAL &out_dpdrhoe,
                                                       CCTK_REAL &out_dpderho,
                                                       CCTK_REAL &out_temp) {
    // Sanity checks
    if (xrho <= 0.)  CCTK_ERROR("xrho <= 0: log(xrho) is undefined (and anyway, not a good sign).");
    if (xeps <= 0.)  CCTK_ERROR("xeps <= 0: log(xeps) is undefined (and anyway, not a good sign).");

    // Convert input to logarithms
    CCTK_REAL xlogrho = log(xrho);

    // Get temperature from xlogrho, xye, xeps
    CCTK_REAL out_logtemp = logtemp_from_logrho_ye_var3D(xlogrho, xye,
                                                         xeps, eps_key);
    out_temp = exp(out_logtemp);

    /* Get d(press)/d(rho) and d(press)/d(eps) from xlogrho, xlogtemp, xye.
       out_vars[0] = d(press)/d(rho)
       out_vars[1] = d(press)/d(eps)                                            */
    CCTK_REAL out_vars[2];

    const vector<key> which_vars3D{dpdrhoe_key, dpderho_key};
    vars3D_from_logrho_logtemp_ye(xlogrho, out_logtemp, xye,
                                  which_vars3D, out_vars);

    // Set output
    out_dpdrhoe = out_vars[0];
    out_dpderho = out_vars[1];

    return;
}
