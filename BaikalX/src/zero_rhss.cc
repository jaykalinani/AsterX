
#include <loop.hxx>    

#include "cctk.h"
#include "cctk_Arguments_Checked.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

void BaikalX_zero_rhss(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_BaikalX_zero_rhss;
  DECLARE_CCTK_PARAMETERS;
Loop::GF3D<CCTK_REAL, 0, 0, 0> aDD00_rhsGF_(cctkGH, aDD00_rhsGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> aDD01_rhsGF_(cctkGH, aDD01_rhsGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> aDD02_rhsGF_(cctkGH, aDD02_rhsGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> aDD11_rhsGF_(cctkGH, aDD11_rhsGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> aDD12_rhsGF_(cctkGH, aDD12_rhsGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> aDD22_rhsGF_(cctkGH, aDD22_rhsGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> alpha_rhsGF_(cctkGH, alpha_rhsGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> betU0_rhsGF_(cctkGH, betU0_rhsGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> betU1_rhsGF_(cctkGH, betU1_rhsGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> betU2_rhsGF_(cctkGH, betU2_rhsGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> cf_rhsGF_(cctkGH, cf_rhsGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> hDD00_rhsGF_(cctkGH, hDD00_rhsGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> hDD01_rhsGF_(cctkGH, hDD01_rhsGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> hDD02_rhsGF_(cctkGH, hDD02_rhsGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> hDD11_rhsGF_(cctkGH, hDD11_rhsGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> hDD12_rhsGF_(cctkGH, hDD12_rhsGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> hDD22_rhsGF_(cctkGH, hDD22_rhsGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> lambdaU0_rhsGF_(cctkGH, lambdaU0_rhsGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> lambdaU1_rhsGF_(cctkGH, lambdaU1_rhsGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> lambdaU2_rhsGF_(cctkGH, lambdaU2_rhsGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> trK_rhsGF_(cctkGH, trK_rhsGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> vetU0_rhsGF_(cctkGH, vetU0_rhsGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> vetU1_rhsGF_(cctkGH, vetU1_rhsGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> vetU2_rhsGF_(cctkGH, vetU2_rhsGF);
Loop::loop_all<0,0,0>(cctkGH, [&](const Loop::PointDesc &p){
    const int i0 = p.i;
    const int i1 = p.j;
    const int i2 = p.k;
aDD00_rhsGF_(i0, i1, i2) = 0.0;
aDD01_rhsGF_(i0, i1, i2) = 0.0;
aDD02_rhsGF_(i0, i1, i2) = 0.0;
aDD11_rhsGF_(i0, i1, i2) = 0.0;
aDD12_rhsGF_(i0, i1, i2) = 0.0;
aDD22_rhsGF_(i0, i1, i2) = 0.0;
alpha_rhsGF_(i0, i1, i2) = 0.0;
betU0_rhsGF_(i0, i1, i2) = 0.0;
betU1_rhsGF_(i0, i1, i2) = 0.0;
betU2_rhsGF_(i0, i1, i2) = 0.0;
cf_rhsGF_(i0, i1, i2) = 0.0;
hDD00_rhsGF_(i0, i1, i2) = 0.0;
hDD01_rhsGF_(i0, i1, i2) = 0.0;
hDD02_rhsGF_(i0, i1, i2) = 0.0;
hDD11_rhsGF_(i0, i1, i2) = 0.0;
hDD12_rhsGF_(i0, i1, i2) = 0.0;
hDD22_rhsGF_(i0, i1, i2) = 0.0;
lambdaU0_rhsGF_(i0, i1, i2) = 0.0;
lambdaU1_rhsGF_(i0, i1, i2) = 0.0;
lambdaU2_rhsGF_(i0, i1, i2) = 0.0;
trK_rhsGF_(i0, i1, i2) = 0.0;
vetU0_rhsGF_(i0, i1, i2) = 0.0;
vetU1_rhsGF_(i0, i1, i2) = 0.0;
vetU2_rhsGF_(i0, i1, i2) = 0.0;

});
}
