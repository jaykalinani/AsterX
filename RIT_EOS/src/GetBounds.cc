/*******************************************************************************
  Routines to get the tabulated values of a some quantity which are closest to
  some chosen value for that quantity, which must be in the EOS table's range
  (otherwise execution is aborted)
********************************************************************************/

#include <cctk.h>
#include "RIT_EOS.hh"

using namespace RIT_EOS_tablevars;



/* Routine to get the indices of the two elements in vector var1D binding the
   value x                                                                      */
void GetBounds_var1D(
    const CCTK_REAL                 &x,      // xlogrho, xlogtemp or xye
          CCTK_REAL *const restrict &var1D,  // logrho,  logtemp  or ye
    const CCTK_INT                  &n,      // nrho,    ntemp    or nye
          CCTK_INT                  &ilow,   // Lower index
          CCTK_INT                  &iup) {  // Upper index
    /* FIXME: sanity check disabled for the sake of performance. This assumes
              that x is in the validity range of var1D, which should have been
              ensured by the interpolator routines BEFORE calling this function.*/
    // Sanity check
    /*if (x < var1D[0] or x > var1D[n - 1])
        CCTK_ERROR("Desired point out of bounds");*/

    // Declaration/setup of variables
    ilow = 0;
    iup  = n - 1;
    CCTK_INT ihalf;

    /* Locate x inside array var1D using binary search
       N.B.: this assumes the elements of var1D are put into ascending order.
             This is checked in ReadEOStable().                                 */
    // TODO: consider using std::binary_search()
    while (ilow < iup - 1) {
        ihalf = (ilow + iup)/2;  // ***** INTEGER division *****
        (x >= var1D[ihalf]) ? (ilow = ihalf) : (iup = ihalf);
    }

    return;
}





/* Routine to get the indices of the two elements in array var3D binding the
   value x in the temperature direction.
   // TODO: extend the list
   x     = xlogpress, xlogenergy or xentropy and
   var3D =  logpress,  logenergy or  entropy.                                   */
void GetBounds_var3D_logrho_ye_tabulated(
    const CCTK_REAL                 &x,      // xlogpress, xlogenergy or xentropy  // TODO: extend the list
          CCTK_REAL *const restrict &var3D,  // logpress,  logenergy  or  entropy  // TODO: extend the list
    const CCTK_INT                  &i,      // logrho index (fixed!)
    const CCTK_INT                  &k,      // ye index     (fixed!)
          CCTK_INT                  &jlow,   // Lower index for logtemp
          CCTK_INT                  &jup) {  // Upper index for logtemp

    // Sanity checks
    if (i < 0 or i >= nrho)
        CCTK_ERROR("Desired interpolation point out of bounds in 'logrho' direction");
    if (k < 0 or k >= nye)
        CCTK_VERROR("Desired interpolation point out of bounds in 'ye' direction");

    /* FIXME: sanity check disabled for the sake of performance. This assumes
              that x is in the validity range of var3D along the logtemp
              direction, which should have been ensured by the interpolator
              routines BEFORE calling this function.                            */
    /*if (x < var3D[k][0][i] or x > var3D[k][ntemp - 1][i])
        CCTK_ERROR("Desired interpolation point out of bounds along in the 3D variable in 'logtemp' direction");*/

    // Declaration/setup of variables
    jlow = 0;
    jup  = ntemp - 1;
    CCTK_INT jhalf;

    /* Locate x inside array var3D[index_ye][:][index_logrho] using binary
       search.
       N.B.: if var3D is pressure or specific internal energy, then the routine
             temp_from_rho_ye_var3D(), which calls this function, takes care of
             re-adding the energy shift (for eps only) and taking log10(x).     */
    while (jlow < jup - 1) {
        jhalf = (jlow + jup)/2;
        (x >= var3D[(k*ntemp + jhalf)*nrho + i]) ?
            (jlow = jhalf) : (jup = jhalf);
    }

    return;
}
