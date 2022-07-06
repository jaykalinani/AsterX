#include <iostream>
#include <cmath>
#include <iomanip>
#include "Lagrange_interp_standalone.hh"

using namespace std;


// ######################### USER-DEFINED PARAMETERS ###########################
#define XMIN         0.
#define XMAX         10.
#define N            100
#define INTERP_ORDER 5
#define X            0.113
// #############################################################################



// Function to interpolate
double f(const double &x) {
    return sin(x);
}



int main() {
    int            i;
    double         x, y_interp;
    double         delta = (XMAX - XMIN)/(1.*N);
    vector<double> v1, v2;

    for (i = 0; i <= N; ++i) {
        x = delta*i;
        v1.push_back(x);
        v2.push_back(f(x));
    }

    y_interp = Lagrange_interp(v1, v2, X, INTERP_ORDER);

    cout << endl << "Analytical   function value: " << setprecision(10) << f(X)
         << endl << "Interpolated function value: " << setprecision(10) << y_interp
         << endl << endl;

    return 0;
}
