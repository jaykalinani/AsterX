
#include <math.h>

#include <loop.hxx>

#include "cctk.h"
#include "cctk_Arguments_Checked.h"
#include "cctk_Parameters.h"

//void BSSN_RHSs()

#define UPWIND_ALG(UpwindVecU) UpwindVecU > 0.0 ? 1.0 : 0.0

void driver_pt2_BSSN_RHSs(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS_driver_pt2_BSSN_RHSs;
    DECLARE_CCTK_PARAMETERS;

    const CCTK_REAL invdx0 = 1.0/CCTK_DELTA_SPACE(0);
    const CCTK_REAL invdx1 = 1.0/CCTK_DELTA_SPACE(1);
    const CCTK_REAL invdx2 = 1.0/CCTK_DELTA_SPACE(2);
Loop::GF3D<std::remove_reference<decltype(*trKGF)>::type,0,0,0> trKGF_(cctkGH,trKGF);
Loop::GF3D<std::remove_reference<decltype(*trK_rhsGF)>::type,0,0,0> trK_rhsGF_(cctkGH,trK_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*alphaGF)>::type,0,0,0> alphaGF_(cctkGH,alphaGF);
Loop::GF3D<std::remove_reference<decltype(*alpha_rhsGF)>::type,0,0,0> alpha_rhsGF_(cctkGH,alpha_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*cfGF)>::type,0,0,0> cfGF_(cctkGH,cfGF);
Loop::GF3D<std::remove_reference<decltype(*cf_rhsGF)>::type,0,0,0> cf_rhsGF_(cctkGH,cf_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*lambdaU0GF)>::type,0,0,0> lambdaU0GF_(cctkGH,lambdaU0GF);
Loop::GF3D<std::remove_reference<decltype(*lambdaU0_rhsGF)>::type,0,0,0> lambdaU0_rhsGF_(cctkGH,lambdaU0_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*vetU0GF)>::type,0,0,0> vetU0GF_(cctkGH,vetU0GF);
Loop::GF3D<std::remove_reference<decltype(*vetU0_rhsGF)>::type,0,0,0> vetU0_rhsGF_(cctkGH,vetU0_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*betU0GF)>::type,0,0,0> betU0GF_(cctkGH,betU0GF);
Loop::GF3D<std::remove_reference<decltype(*betU0_rhsGF)>::type,0,0,0> betU0_rhsGF_(cctkGH,betU0_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*hDD00GF)>::type,0,0,0> hDD00GF_(cctkGH,hDD00GF);
Loop::GF3D<std::remove_reference<decltype(*aDD00GF)>::type,0,0,0> aDD00GF_(cctkGH,aDD00GF);
Loop::GF3D<std::remove_reference<decltype(*hDD00_rhsGF)>::type,0,0,0> hDD00_rhsGF_(cctkGH,hDD00_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*aDD00_rhsGF)>::type,0,0,0> aDD00_rhsGF_(cctkGH,aDD00_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*RbarDD00GF)>::type,0,0,0> RbarDD00GF_(cctkGH,RbarDD00GF);
Loop::GF3D<std::remove_reference<decltype(*hDD01GF)>::type,0,0,0> hDD01GF_(cctkGH,hDD01GF);
Loop::GF3D<std::remove_reference<decltype(*aDD01GF)>::type,0,0,0> aDD01GF_(cctkGH,aDD01GF);
Loop::GF3D<std::remove_reference<decltype(*hDD01_rhsGF)>::type,0,0,0> hDD01_rhsGF_(cctkGH,hDD01_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*aDD01_rhsGF)>::type,0,0,0> aDD01_rhsGF_(cctkGH,aDD01_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*RbarDD01GF)>::type,0,0,0> RbarDD01GF_(cctkGH,RbarDD01GF);
Loop::GF3D<std::remove_reference<decltype(*hDD02GF)>::type,0,0,0> hDD02GF_(cctkGH,hDD02GF);
Loop::GF3D<std::remove_reference<decltype(*aDD02GF)>::type,0,0,0> aDD02GF_(cctkGH,aDD02GF);
Loop::GF3D<std::remove_reference<decltype(*hDD02_rhsGF)>::type,0,0,0> hDD02_rhsGF_(cctkGH,hDD02_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*aDD02_rhsGF)>::type,0,0,0> aDD02_rhsGF_(cctkGH,aDD02_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*RbarDD02GF)>::type,0,0,0> RbarDD02GF_(cctkGH,RbarDD02GF);
Loop::GF3D<std::remove_reference<decltype(*lambdaU1GF)>::type,0,0,0> lambdaU1GF_(cctkGH,lambdaU1GF);
Loop::GF3D<std::remove_reference<decltype(*lambdaU1_rhsGF)>::type,0,0,0> lambdaU1_rhsGF_(cctkGH,lambdaU1_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*vetU1GF)>::type,0,0,0> vetU1GF_(cctkGH,vetU1GF);
Loop::GF3D<std::remove_reference<decltype(*vetU1_rhsGF)>::type,0,0,0> vetU1_rhsGF_(cctkGH,vetU1_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*betU1GF)>::type,0,0,0> betU1GF_(cctkGH,betU1GF);
Loop::GF3D<std::remove_reference<decltype(*betU1_rhsGF)>::type,0,0,0> betU1_rhsGF_(cctkGH,betU1_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*hDD11GF)>::type,0,0,0> hDD11GF_(cctkGH,hDD11GF);
Loop::GF3D<std::remove_reference<decltype(*aDD11GF)>::type,0,0,0> aDD11GF_(cctkGH,aDD11GF);
Loop::GF3D<std::remove_reference<decltype(*hDD11_rhsGF)>::type,0,0,0> hDD11_rhsGF_(cctkGH,hDD11_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*aDD11_rhsGF)>::type,0,0,0> aDD11_rhsGF_(cctkGH,aDD11_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*RbarDD11GF)>::type,0,0,0> RbarDD11GF_(cctkGH,RbarDD11GF);
Loop::GF3D<std::remove_reference<decltype(*hDD12GF)>::type,0,0,0> hDD12GF_(cctkGH,hDD12GF);
Loop::GF3D<std::remove_reference<decltype(*aDD12GF)>::type,0,0,0> aDD12GF_(cctkGH,aDD12GF);
Loop::GF3D<std::remove_reference<decltype(*hDD12_rhsGF)>::type,0,0,0> hDD12_rhsGF_(cctkGH,hDD12_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*aDD12_rhsGF)>::type,0,0,0> aDD12_rhsGF_(cctkGH,aDD12_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*RbarDD12GF)>::type,0,0,0> RbarDD12GF_(cctkGH,RbarDD12GF);
Loop::GF3D<std::remove_reference<decltype(*lambdaU2GF)>::type,0,0,0> lambdaU2GF_(cctkGH,lambdaU2GF);
Loop::GF3D<std::remove_reference<decltype(*lambdaU2_rhsGF)>::type,0,0,0> lambdaU2_rhsGF_(cctkGH,lambdaU2_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*vetU2GF)>::type,0,0,0> vetU2GF_(cctkGH,vetU2GF);
Loop::GF3D<std::remove_reference<decltype(*vetU2_rhsGF)>::type,0,0,0> vetU2_rhsGF_(cctkGH,vetU2_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*betU2GF)>::type,0,0,0> betU2GF_(cctkGH,betU2GF);
Loop::GF3D<std::remove_reference<decltype(*betU2_rhsGF)>::type,0,0,0> betU2_rhsGF_(cctkGH,betU2_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*hDD22GF)>::type,0,0,0> hDD22GF_(cctkGH,hDD22GF);
Loop::GF3D<std::remove_reference<decltype(*aDD22GF)>::type,0,0,0> aDD22GF_(cctkGH,aDD22GF);
Loop::GF3D<std::remove_reference<decltype(*hDD22_rhsGF)>::type,0,0,0> hDD22_rhsGF_(cctkGH,hDD22_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*aDD22_rhsGF)>::type,0,0,0> aDD22_rhsGF_(cctkGH,aDD22_rhsGF);
Loop::GF3D<std::remove_reference<decltype(*RbarDD22GF)>::type,0,0,0> RbarDD22GF_(cctkGH,RbarDD22GF);


#include "BSSN_RHSs.h"
}
