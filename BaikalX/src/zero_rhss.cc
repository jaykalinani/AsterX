
#include <loop.hxx>    

#include "cctk.h"
#include "cctk_Arguments_Checked.h"
#include "cctk_Parameters.h"

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
Loop::GF3D<CCTK_REAL, 0, 0, 0> aDD00GF_(cctkGH, aDD00GF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> aDD01GF_(cctkGH, aDD01GF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> aDD02GF_(cctkGH, aDD02GF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> aDD11GF_(cctkGH, aDD11GF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> aDD12GF_(cctkGH, aDD12GF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> aDD22GF_(cctkGH, aDD22GF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> alphaGF_(cctkGH, alphaGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> betU0GF_(cctkGH, betU0GF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> betU1GF_(cctkGH, betU1GF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> betU2GF_(cctkGH, betU2GF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> cfGF_(cctkGH, cfGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> hDD00GF_(cctkGH, hDD00GF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> hDD01GF_(cctkGH, hDD01GF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> hDD02GF_(cctkGH, hDD02GF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> hDD11GF_(cctkGH, hDD11GF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> hDD12GF_(cctkGH, hDD12GF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> hDD22GF_(cctkGH, hDD22GF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> lambdaU0GF_(cctkGH, lambdaU0GF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> lambdaU1GF_(cctkGH, lambdaU1GF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> lambdaU2GF_(cctkGH, lambdaU2GF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> trKGF_(cctkGH, trKGF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> vetU0GF_(cctkGH, vetU0GF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> vetU1GF_(cctkGH, vetU1GF);
Loop::GF3D<CCTK_REAL, 0, 0, 0> vetU2GF_(cctkGH, vetU2GF);
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
Loop::loop_bnd<0,0,0>(cctkGH, [&](const Loop::PointDesc &p){
    const int i0 = p.i;
    const int i1 = p.j;
    const int i2 = p.k;
aDD00GF_(i0, i1, i2) = 0.0;
aDD01GF_(i0, i1, i2) = 0.0;
aDD02GF_(i0, i1, i2) = 0.0;
aDD11GF_(i0, i1, i2) = 0.0;
aDD12GF_(i0, i1, i2) = 0.0;
aDD22GF_(i0, i1, i2) = 0.0;
alphaGF_(i0, i1, i2) = 0.0;
betU0GF_(i0, i1, i2) = 0.0;
betU1GF_(i0, i1, i2) = 0.0;
betU2GF_(i0, i1, i2) = 0.0;
cfGF_(i0, i1, i2) = 0.0;
hDD00GF_(i0, i1, i2) = 0.0;
hDD01GF_(i0, i1, i2) = 0.0;
hDD02GF_(i0, i1, i2) = 0.0;
hDD11GF_(i0, i1, i2) = 0.0;
hDD12GF_(i0, i1, i2) = 0.0;
hDD22GF_(i0, i1, i2) = 0.0;
lambdaU0GF_(i0, i1, i2) = 0.0;
lambdaU1GF_(i0, i1, i2) = 0.0;
lambdaU2GF_(i0, i1, i2) = 0.0;
trKGF_(i0, i1, i2) = 0.0;
vetU0GF_(i0, i1, i2) = 0.0;
vetU1GF_(i0, i1, i2) = 0.0;
vetU2GF_(i0, i1, i2) = 0.0;

});
}
