#ifndef SHTOOLS_H
#define SHTOOLS_H

#include <cctk.h>

#ifdef __cplusplus
extern "C" {
#endif

#define C_SHExpandDH c_shexpanddh_
void C_SHExpandDH(const double *griddh, const int *n, double *cilm, int *lmax,
                  const int *norm, const int *sampling, const int *csphase,
                  const int *lmax_calc, int *exitstatus);

#define C_MakeGridDH c_makegriddh_
void C_MakeGridDH(double *griddh, int *n, const double *cilm, const int *lmax,
                  const int *norm, const int *sampling, const int *csphase,
                  const int *lmax_calc, int *exitstatus);

#ifdef __cplusplus
}
#endif

#endif // #ifndef SHTOOLS_H
