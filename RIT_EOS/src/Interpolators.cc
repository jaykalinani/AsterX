/*******************************************************************************
    Routines to interpolate the EOS table
********************************************************************************/

#include <cmath>  // for log(), exp() and fabs()
#include <vector>
#include <cctk.h>
#include <cctk_Parameters.h>  // For bisection_maxit and bisection_eps
#include "RIT_EOS.hh"

using std::min;
using std::max;
using std::vector;
using namespace RIT_EOS_tablevars;
using namespace RIT_EOS_keys;



/* Routine to retrieve values from 3D quantities for values of logrho, logtemp,
   ye which are not tabulated; the cheapest way to do this is trilinear
   interpolation.                                                               */
void vars3D_from_logrho_logtemp_ye(
    const CCTK_REAL   &xlogrho,
    const CCTK_REAL   &xlogtemp,
    const CCTK_REAL   &xye,
    const vector<key> &which_vars3D,
          CCTK_REAL   *restrict out_vars3D) {

    // Declaration of variables
    CCTK_INT   count;
    CCTK_INT   ilow, iup, jlow, jup, klow, kup;
    CCTK_INT   klow_ntemp, kup_ntemp;
    CCTK_INT   klow_ntemp_p_jlow_nrho, klow_ntemp_p_jup_nrho;
    CCTK_INT   kup_ntemp_p_jlow_nrho,  kup_ntemp_p_jup_nrho;
    CCTK_REAL  deltalow, deltaup, deltatot_inv;
    CCTK_REAL  aux1, aux2, aux3, aux4, aux5, aux6;
    CCTK_REAL *var3D;


    // Sanity checks
    if (xlogrho < logrho[0] or xlogrho > logrho[nrho_m1])
        CCTK_VERROR("Desired interpolation point (xlogrho = %g) out of bounds [%g, %g] in 'logrho' direction",
                    xlogrho, logrho[0], logrho[nrho_m1]);

    if (xlogtemp < logtemp[0] or xlogtemp > logtemp[ntemp_m1])
        CCTK_VERROR("Desired interpolation point (xlogtemp = %g) out of bounds [%g, %g] in 'logtemp' direction",
                    xlogtemp, logtemp[0], logtemp[ntemp_m1]);

    if (xye < ye[0] or xye > ye[nye_m1])
        CCTK_VERROR("Desired interpolation point (xye = %g) out of bounds [%g, %g] in 'ye' direction",
                    xye, ye[0], ye[nye_m1]);

    if (out_vars3D == nullptr)
        CCTK_ERROR("out_vars3D is a nullptr.");


    // Select the 3D quantities to interpolate and perform the interpolations
    count = 0;

    for (const auto &this_key : which_vars3D) {
        var3D = this_key.var3D;

        // Get bounds for logrho, logtemp and ye
        GetBounds_var1D(xlogrho,  logrho,  nrho,  ilow, iup);
        GetBounds_var1D(xlogtemp, logtemp, ntemp, jlow, jup);
        GetBounds_var1D(xye,      ye,      nye,   klow, kup);

        // Set up helper variables
        klow_ntemp             = klow*ntemp;
        kup_ntemp              = kup*ntemp;
        klow_ntemp_p_jlow_nrho = (klow_ntemp + jlow)*nrho;
        klow_ntemp_p_jup_nrho  = (klow_ntemp + jup)*nrho;
        kup_ntemp_p_jlow_nrho  = (kup_ntemp  + jlow)*nrho;
        kup_ntemp_p_jup_nrho   = (kup_ntemp  + jup)*nrho;

        // Partial interpolations for fixed logrho
        deltalow     = xlogrho - logrho[ilow];
        deltaup      = xlogrho - logrho[iup];
        deltatot_inv = 1./(logrho[iup] - logrho[ilow]);
        aux1         = (var3D[klow_ntemp_p_jlow_nrho + iup]*deltalow  - var3D[klow_ntemp_p_jlow_nrho + ilow]*deltaup)*deltatot_inv;
        aux2         = (var3D[klow_ntemp_p_jup_nrho  + iup]*deltalow  - var3D[klow_ntemp_p_jup_nrho  + ilow]*deltaup)*deltatot_inv;
        aux3         = (var3D[kup_ntemp_p_jlow_nrho  + iup]*deltalow  - var3D[kup_ntemp_p_jlow_nrho  + ilow]*deltaup)*deltatot_inv;
        aux4         = (var3D[kup_ntemp_p_jup_nrho   + iup]*deltalow  - var3D[kup_ntemp_p_jup_nrho   + ilow]*deltaup)*deltatot_inv;

        // Partial interpolations for fixed logtemp
        deltalow     = xlogtemp - logtemp[jlow];
        deltaup      = xlogtemp - logtemp[jup];
        deltatot_inv = 1./(logtemp[jup] - logtemp[jlow]);
        aux5         = deltatot_inv*(deltalow*aux2 - deltaup*aux1);
        aux6         = deltatot_inv*(deltalow*aux4 - deltaup*aux3);

        // Final interpolation (fixed ye)
        out_vars3D[count] = (aux6*(xye - ye[klow]) - aux5*(xye - ye[kup]))/(ye[kup] - ye[klow]);

        /* Exponentiate and/or subtract the energy shift if needed.
           N.B.: energy_shift is subtracted AFTER exponentiation!               */
        if (this_key.do_exp_log)       out_vars3D[count]  = exp(out_vars3D[count]);
        if (this_key.do_energy_shift)  out_vars3D[count] -= energy_shift;

        ++count;
    }

    return;
}





/* Routine to retrieve the temperature knowing rho, ye, and one 3D variable.
   N.B.: xvar3D is NOT const: the energy shift is re-added to it if
         xvar3D = xeps and the log is taken is xvar3D = xpress or xvar3D = xeps.*/
CCTK_REAL logtemp_from_logrho_ye_var3D(
    const CCTK_REAL &xlogrho,
    const CCTK_REAL &xye,
          CCTK_REAL xvar3D,
    const key       &which_var3D) {

    DECLARE_CCTK_PARAMETERS;

    // Declaration of variables
    CCTK_REAL *var3D;
    CCTK_REAL xvar3D_as_input;
    const vector<key> which_vars3D{which_var3D};
    CCTK_INT  ilow, iup, klow, kup;
    CCTK_INT  klow_ntemp,                 kup_ntemp;
    CCTK_INT  klow_ntemp_nrho,            kup_ntemp_nrho;
    CCTK_INT  klow_ntemp_p_ntemp_m1_nrho, kup_ntemp_p_ntemp_m1_nrho;
    CCTK_INT  jlow1, jup1, jlow2, jup2, jlow3, jup3, jlow4, jup4;
    CCTK_INT  j, jlow, jup, n;
    CCTK_REAL out_logtemp, f;
    CCTK_REAL out_logtemp_bound_low, out_logtemp_bound_up;
    CCTK_REAL var3D_jlow, var3D_j;
    bool      is_var3D_jlow_LT_xvar3D_as_input;


    // Sanity checks
    if (xlogrho < logrho[0] or xlogrho > logrho[nrho_m1])
        CCTK_VERROR("Desired interpolation point (xlogrho = %g) out of bounds [%g, %g] in 'logrho' direction",
                    xlogrho, logrho[0], logrho[nrho_m1]);

    if (xye < ye[0] or xye > ye[nye_m1])
        CCTK_VERROR("Desired interpolation point (xye = %g) out of bounds [%g, %g] in 'ye' direction",
                    xye, ye[0], ye[nye_m1]);


    /* Select the 3D quantity to use to retrieve the temperature and re-add the
       energy shift and/or take the log if needed. Save a copy of xvar3D as it
       was input in the variable 'xvar3D_as_input'; this may be needed when
       retrieving the temperature (see below).
       N.B.: the energy shift is added BEFORE taking the log!                 */
    xvar3D_as_input = xvar3D;
    var3D           = which_var3D.var3D;

    if (which_var3D.do_energy_shift)  xvar3D += energy_shift;
    if (which_var3D.do_exp_log) {
        if (xvar3D > 0.) xvar3D = log(xvar3D);
        else CCTK_ERROR("xvar3D <= 0, can't take the log");
    }


    // Get bounds for logrho and ye
    GetBounds_var1D(xlogrho, logrho, nrho, ilow, iup);
    GetBounds_var1D(xye,     ye,     nye,  klow, kup);


    // Set up helper variables
    klow_ntemp = klow*ntemp;
    kup_ntemp  = kup*ntemp;
    klow_ntemp_nrho            = klow_ntemp*nrho;  // (klow_ntemp + 0)*nrho;
    klow_ntemp_p_ntemp_m1_nrho = (klow_ntemp + ntemp_m1)*nrho;
    kup_ntemp_nrho             = kup_ntemp*nrho;   // (kup_ntemp + 0)*nrho;
    kup_ntemp_p_ntemp_m1_nrho  = (kup_ntemp + ntemp_m1)*nrho;


    // Additional sanity checks
    if (xvar3D < var3D[klow_ntemp_nrho            + ilow] or  // [klow][0][ilow]
        xvar3D > var3D[klow_ntemp_p_ntemp_m1_nrho + ilow])    // [klow][ntemp_m1][ilow]
        CCTK_VERROR("Desired interpolation point (xvar3D = %g) out of bounds [%g, %g] in 'var3D' at [klow, ilow] along the 'logtemp' direction",
                    xvar3D,
                    var3D[klow_ntemp_nrho            + ilow],
                    var3D[klow_ntemp_p_ntemp_m1_nrho + ilow]);

    if (xvar3D < var3D[klow_ntemp_nrho            + iup] or   // [klow][0][iup]
        xvar3D > var3D[klow_ntemp_p_ntemp_m1_nrho + iup])     // [klow][ntemp_m1][iup]
        CCTK_VERROR("Desired interpolation point (xvar3D = %g) out of bounds [%g, %g] in 'var3D' at [klow, iup] along the 'logtemp' direction",
                    xvar3D,
                    var3D[klow_ntemp_nrho            + iup],
                    var3D[klow_ntemp_p_ntemp_m1_nrho + iup]);

    if (xvar3D < var3D[kup_ntemp_nrho             + ilow] or  // [kup][0][ilow]
        xvar3D > var3D[kup_ntemp_p_ntemp_m1_nrho  + ilow])    // [kup][ntemp_m1][ilow]
        CCTK_VERROR("Desired interpolation point (xvar3D = %g) out of bounds [%g, %g] in 'var3D' at [kup, ilow] along the 'logtemp' direction",
                    xvar3D,
                    var3D[kup_ntemp_nrho            + ilow],
                    var3D[kup_ntemp_p_ntemp_m1_nrho + ilow]);

    if (xvar3D < var3D[kup_ntemp_nrho             + iup] or   // [kup][0][iup]
        xvar3D > var3D[kup_ntemp_p_ntemp_m1_nrho  + iup])     // [kup][ntemp_m1][iup]
        CCTK_VERROR("Desired interpolation point (xvar3D = %g) out of bounds [%g, %g] in 'var3D' at [kup, iup] along the 'logtemp' direction",
                    xvar3D,
                    var3D[kup_ntemp_nrho            + iup],
                    var3D[kup_ntemp_p_ntemp_m1_nrho + iup]);


    /* Get bounds for var3D on four different lines of varying logtemp:
       1. for fixed (logrho[ilow], ye[klow])
       2. for fixed (logrho[ilow], ye[kup])
       3. for fixed (logrho[iup], ye[klow])
       4. for fixed (logrho[iup], ye[kup])                                      */
    GetBounds_var3D_logrho_ye_tabulated(xvar3D, var3D, ilow, klow, jlow1, jup1);
    GetBounds_var3D_logrho_ye_tabulated(xvar3D, var3D, ilow, kup,  jlow2, jup2);
    GetBounds_var3D_logrho_ye_tabulated(xvar3D, var3D, iup,  klow, jlow3, jup3);
    GetBounds_var3D_logrho_ye_tabulated(xvar3D, var3D, iup,  kup,  jlow4, jup4);


    // Set bounds to bracket out_logtemp
    jlow = min(min(jlow1, jlow2), min(jlow3, jlow4));
    jup  = max(max(jup1,  jup2),  max(jup3,  jup4));

    // If jup = jlow + 1, then out_logtemp is already bracketed
    if (jup == jlow + 1) {
        out_logtemp_bound_low = logtemp[jlow];
        out_logtemp_bound_up  = logtemp[jup];
    }

    /* If not, move the candidate upper bound on out_logtemp -- called j -- to
       the "right" until out_logtemp is bracketed                               */
    else {
        /* Get the value 'var3D_jlow' of the 3D variable corresponding to the
           temporary lower bound on out_logtemp -- jlow                         */ 
        vars3D_from_logrho_logtemp_ye(xlogrho, logtemp[jlow], xye,
                                      which_vars3D, &var3D_jlow);

        /* is_var3D_jlow_LT_xvar3D_as_input = true   if var3D_jlow < xvar3D_as_input
           is_var3D_jlow_LT_xvar3D_as_input = false  otherwise                  */
        is_var3D_jlow_LT_xvar3D_as_input = (var3D_jlow < xvar3D_as_input);

        /* Get the value 'var3D_j' of the 3D variable corresponding to the
           index just above jlow, the temporary lower bound on out_logtemp      */ 
        j = jlow + 1;
        vars3D_from_logrho_logtemp_ye(xlogrho, logtemp[j], xye,
                                      which_vars3D, &var3D_j);

        /* Keep on adding 1 to j (i.e., keep on moving the candidate upper bound
           on out_logtemp to the "right") as long as both var3D_jlow and
           var3D_j are smaller or larger than xvar3D (i.e., as long as
           out_logtemp is not bracketed)                                        */           
        // FIXME: maybe a do-while is better here, but not sure (probably not)
        while((var3D_j < xvar3D_as_input) == is_var3D_jlow_LT_xvar3D_as_input) {
            ++j;
            vars3D_from_logrho_logtemp_ye(xlogrho, logtemp[j], xye,
                                          which_vars3D, &var3D_j);
        }

        // At this point, the lower bound is for sure j - 1
        out_logtemp_bound_low = logtemp[j - 1];
        out_logtemp_bound_up  = logtemp[j];
    }


    /* Find out_logtemp such that var3D(xye, out_logtemp, xlogrho) = xvar3D
       using bisection                                                          */
    for (n = 0; n < bisection_maxit; ++n) {
        out_logtemp = 0.5*(out_logtemp_bound_low + out_logtemp_bound_up);
        if (fabs(out_logtemp_bound_up - out_logtemp_bound_low) < bisection_eps)
            return out_logtemp;

        vars3D_from_logrho_logtemp_ye(xlogrho, out_logtemp, xye,
                                      which_vars3D, &f);

        if      (f < xvar3D_as_input)  out_logtemp_bound_low = out_logtemp;
        else if (f > xvar3D_as_input)  out_logtemp_bound_up  = out_logtemp;
        else                           return out_logtemp;
    }

    CCTK_VERROR("Bisection failed (more than %d iterations done without converging to precision %g)",
                bisection_maxit, bisection_eps);
}
