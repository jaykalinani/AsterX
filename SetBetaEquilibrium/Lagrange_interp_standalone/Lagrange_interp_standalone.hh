#ifndef LAGRANGE_INTERP_STANDALONE
#define LAGRANGE_INTERP_STANDALONE

#include <vector>

using namespace std;

double Lagrange_interp(const vector<double> &xvec,
                       const vector<double> &yvec,
                       const double         &x,
                       const int            &n);

#endif // LAGRANGE_INTERP_STANDALONE
