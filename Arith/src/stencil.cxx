#define NSIMD_AVX2
#define NSIMD_FMA
#include <nsimd/nsimd-all.hpp>

void nsimd_deriv(double *__restrict__ const du,
                 const double *__restrict__ const u, const int n) {
  typedef nsimd::pack<double> vdouble;
  typedef nsimd::packl<double> vbool;
  const int vsize = nsimd::len(vdouble{});
  for (int i = 2; i < n - 2; i += vsize) {
    vbool mask = nsimd::mask_for_loop_tail<vbool>(i, n);
    vdouble um2 = nsimd::maskz_loadu(mask, &u[i - 2]);
    vdouble um1 = nsimd::maskz_loadu(mask, &u[i - 1]);
    vdouble up1 = nsimd::maskz_loadu(mask, &u[i + 1]);
    vdouble up2 = nsimd::maskz_loadu(mask, &u[i + 2]);
    nsimd::mask_storeu(mask, &du[i], 3 * (um2 - up2) + 2 * (um1 - up1));
  }
}
