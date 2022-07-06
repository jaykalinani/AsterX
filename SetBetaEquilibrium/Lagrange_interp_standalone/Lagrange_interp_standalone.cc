#include <iostream>
#include <vector>
#include "Lagrange_interp_standalone.hh"

using namespace std;



// Lagrange interpolation
double Lagrange_interp(const vector<double> &xvec,
                       const vector<double> &yvec,
                       const double         &x,
                       const int            &n) {
    // Declaration of variables
    int i, j;
    int n_half = n/2;

    const int s = xvec.size();

    int i_left  = 0;
    int i_right = s - 1;
    int i_mid;

    vector<double> xpoints, ypoints;  // Interpolation windows along xvec and yvec

    double y_interp = 0.;
    double y_interp_tmp;



    // Sanity checks
    if (x < xvec.front() or x > xvec.back()) {
        cout << endl
             << "Target interpolation point " << x << " is out of bounds ["
             << xvec.front() << ", " << xvec.back() << "]."
             << endl << endl;
        exit(EXIT_FAILURE);
    }

    if (s != yvec.size()) {
        cout << endl
             << "xvec and yvec have different sizes."
             << endl << endl;
        exit(EXIT_FAILURE);
    }

    if (n < 1) {
        cout << endl
             << "Interpolation order set to " << n << " < 1; makes no sense."
             << endl << endl;
        exit(EXIT_FAILURE);
    }

    if (s < n) {
        cout << endl
             << "Too few points (" << s << ") to perform a " << n
             << "-point Lagrange interpolation"
             << endl << endl;
        exit(EXIT_FAILURE);
    }

    if (s < 1) {
        cout << endl
             << "xvec contains no elements."
             << endl << endl;
        exit(EXIT_FAILURE);
    }

    for (i = 1; i < s; ++i) {
        if (xvec.at(i) <= xvec.at(i - 1)) {
            cout << endl
                 << "xvec is not sorted into ascending order. Both binary search and Lagrange interpolation are not gonna work."
                 << endl << endl;
            exit(EXIT_FAILURE);
        }
    }



    /* Find the indices corresponding to the elements of xvec which are closest
       to x from the left (i_left) and from the right (i_right) with binary
       search                                                                   */
    while (i_left < i_right - 1) {
        i_mid = (i_left + i_right)/2;
        (x >= xvec.at(i_mid)) ? (i_left = i_mid) : (i_right = i_mid);
    }



    /* Build the interpolation window along xvec
       N.B.: by construction, after binary search, iright = ileft+1             */

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
    if (xpoints.size() != n) {
        cout << "There are " << xpoints.size() <<
                " points in the interpolation window, but they should be "
             << n << ", as the order of Lagrange interpolation."
             << endl << endl;
        exit(EXIT_FAILURE);
    }



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
