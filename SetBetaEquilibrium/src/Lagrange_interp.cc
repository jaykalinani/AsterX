/* Implmentation of the Lagrange interpolation routine used in
   SetBetaEquilibrium_Compose() and SetBetaEquilibrium_Lorene()                 */

#include <vector>
#include "SetBetaEquilibrium.hh"

using std::vector;



CCTK_REAL Lagrange_interp(const vector<CCTK_REAL> &xvec,
                          const vector<CCTK_REAL> &yvec,
                          const CCTK_REAL         &x,
                          const CCTK_INT          &n) {
    // Declaration of variables
          CCTK_INT i, j;
          CCTK_INT n_half = n/2;
    const CCTK_INT s      = xvec.size();

    CCTK_INT i_left  = 0;
    CCTK_INT i_right = s - 1;
    CCTK_INT i_mid;

    vector<CCTK_REAL> xpoints, ypoints;  // Interpolation windows along xvec and yvec

    CCTK_REAL y_interp = 0.;
    CCTK_REAL y_interp_tmp;



    // Sanity checks
    if (x < xvec.front() or x > xvec.back())
        CCTK_VERROR("Target interpolation point %f is out of bounds [%f, %f].",
                    x, xvec.front(), xvec.back());

    if (s != yvec.size())
        CCTK_ERROR("xvec and yvec have different sizes.");

    if (n < 1)
        CCTK_VERROR("Interpolation order set to %d < 1; makes no sense.", n);

    if (s < n)
        CCTK_VERROR("Too few points (%d) to perform a %d-point Lagrange interpolation",
            s, n);

    if (s < 1)
        CCTK_ERROR("xvec contains no elements.");

    for (i = 1; i < s; ++i)
        if (xvec.at(i) <= xvec.at(i - 1))
            CCTK_ERROR("xvec is not sorted into ascending order. Both binary search and Lagrange interpolation are not gonna work.");
            


    /* Find the indices corresponding to the elements of xvec which are closest
       to x from the left (i_left) and from the right (i_right) with binary
       search                                                                   */
    while (i_left < i_right - 1) {
        i_mid = (i_left + i_right)/2;
        (x >= xvec.at(i_mid)) ? (i_left = i_mid) : (i_right = i_mid);
    }



    /* Build the interpolation window along xvec
       N.B.: by construction, after binary search, iright = ileft + 1           */

    // If n is even
    if ((n & 1) == 0) {
        // If x is near the beginning of xvec
        if (i_right < n_half)
            for (j = 0; j < n; ++j) {
                xpoints.push_back(xvec.at(j));
                ypoints.push_back(yvec.at(j));
            }

        // If x is near the end of xvec
        else if (i_right > s - n_half)
            for (j = s - n; j < s; ++j) {
                xpoints.push_back(xvec.at(j));
                ypoints.push_back(yvec.at(j));
            }

        // If x is well inside xvec
        else
            for (j = i_right - n_half; j < i_right + n_half; ++j) {
                xpoints.push_back(xvec.at(j));
                ypoints.push_back(yvec.at(j));
            }
    }

    // If n is odd
    else {
        // If x is near the beginning of xvec
        if (i_left < n_half)
            for (i = 0; i < n; ++i) {
                xpoints.push_back(xvec.at(i));
                ypoints.push_back(yvec.at(i));
            }

        // If x is near the end of xvec
        else if (i_right >= s - n_half)
            for (i = s - n; i < s; ++i) {
                xpoints.push_back(xvec.at(i));
                ypoints.push_back(yvec.at(i));
            }

        // If x is well inside xvec
        else
            for (i = i_left - n_half; i <= i_left + n_half; ++i) {
                xpoints.push_back(xvec.at(i));
                ypoints.push_back(yvec.at(i));
            }
    }



    // Additional sanity check
    if (xpoints.size() != n)
        CCTK_VINFO("There are %d points in the interpolation window, but they should be %d, as the order of Lagrange interpolation.",
                   xpoints.size(), n);



    // Perform the n-points Lagrange interpolation
    for (i = 0; i < n; ++i) {
        y_interp_tmp = ypoints.at(i);

        for (j = 0; j < n; ++j)
            if (j != i)
                y_interp_tmp *= (x - xpoints.at(j))/(xpoints.at(i) - xpoints.at(j));

        y_interp += y_interp_tmp;
    }

    return y_interp;
}
