
#include <math.h>

#include <loop.hxx>

#include "cctk.h"
#include "cctk_Arguments_Checked.h"
#include "cctk_Parameters.h"

void BaikalX_ADM_to_BSSN_all_but_lambdaU(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS_BaikalX_ADM_to_BSSN_all_but_lambdaU;
    DECLARE_CCTK_PARAMETERS;
    
const CCTK_REAL *alphaSphorCartGF = alp;
    const CCTK_REAL *BSphorCartU0GF = dtbetax;
    const CCTK_REAL *BSphorCartU1GF = dtbetay;
    const CCTK_REAL *BSphorCartU2GF = dtbetaz;
    const CCTK_REAL *KSphorCartDD00GF = kxx;
    const CCTK_REAL *KSphorCartDD01GF = kxy;
    const CCTK_REAL *KSphorCartDD02GF = kxz;
    const CCTK_REAL *KSphorCartDD11GF = kyy;
    const CCTK_REAL *KSphorCartDD12GF = kyz;
    const CCTK_REAL *KSphorCartDD22GF = kzz;
    const CCTK_REAL *betaSphorCartU0GF = betax;
    const CCTK_REAL *betaSphorCartU1GF = betay;
    const CCTK_REAL *betaSphorCartU2GF = betaz;
    const CCTK_REAL *gammaSphorCartDD00GF = gxx;
    const CCTK_REAL *gammaSphorCartDD01GF = gxy;
    const CCTK_REAL *gammaSphorCartDD02GF = gxz;
    const CCTK_REAL *gammaSphorCartDD11GF = gyy;
    const CCTK_REAL *gammaSphorCartDD12GF = gyz;
    const CCTK_REAL *gammaSphorCartDD22GF = gzz;
Loop::GF3D<std::remove_reference<decltype(*alphaSphorCartGF)>::type,0,0,0> alphaSphorCartGF_(cctkGH,alphaSphorCartGF);
Loop::GF3D<std::remove_reference<decltype(*trKGF)>::type,0,0,0> trKGF_(cctkGH,trKGF);
Loop::GF3D<std::remove_reference<decltype(*alphaGF)>::type,0,0,0> alphaGF_(cctkGH,alphaGF);
Loop::GF3D<std::remove_reference<decltype(*cfGF)>::type,0,0,0> cfGF_(cctkGH,cfGF);
Loop::GF3D<std::remove_reference<decltype(*betaSphorCartU0GF)>::type,0,0,0> betaSphorCartU0GF_(cctkGH,betaSphorCartU0GF);
Loop::GF3D<std::remove_reference<decltype(*BSphorCartU0GF)>::type,0,0,0> BSphorCartU0GF_(cctkGH,BSphorCartU0GF);
Loop::GF3D<std::remove_reference<decltype(*vetU0GF)>::type,0,0,0> vetU0GF_(cctkGH,vetU0GF);
Loop::GF3D<std::remove_reference<decltype(*betU0GF)>::type,0,0,0> betU0GF_(cctkGH,betU0GF);
Loop::GF3D<std::remove_reference<decltype(*gammaSphorCartDD00GF)>::type,0,0,0> gammaSphorCartDD00GF_(cctkGH,gammaSphorCartDD00GF);
Loop::GF3D<std::remove_reference<decltype(*KSphorCartDD00GF)>::type,0,0,0> KSphorCartDD00GF_(cctkGH,KSphorCartDD00GF);
Loop::GF3D<std::remove_reference<decltype(*hDD00GF)>::type,0,0,0> hDD00GF_(cctkGH,hDD00GF);
Loop::GF3D<std::remove_reference<decltype(*aDD00GF)>::type,0,0,0> aDD00GF_(cctkGH,aDD00GF);
Loop::GF3D<std::remove_reference<decltype(*gammaSphorCartDD01GF)>::type,0,0,0> gammaSphorCartDD01GF_(cctkGH,gammaSphorCartDD01GF);
Loop::GF3D<std::remove_reference<decltype(*KSphorCartDD01GF)>::type,0,0,0> KSphorCartDD01GF_(cctkGH,KSphorCartDD01GF);
Loop::GF3D<std::remove_reference<decltype(*hDD01GF)>::type,0,0,0> hDD01GF_(cctkGH,hDD01GF);
Loop::GF3D<std::remove_reference<decltype(*aDD01GF)>::type,0,0,0> aDD01GF_(cctkGH,aDD01GF);
Loop::GF3D<std::remove_reference<decltype(*gammaSphorCartDD02GF)>::type,0,0,0> gammaSphorCartDD02GF_(cctkGH,gammaSphorCartDD02GF);
Loop::GF3D<std::remove_reference<decltype(*KSphorCartDD02GF)>::type,0,0,0> KSphorCartDD02GF_(cctkGH,KSphorCartDD02GF);
Loop::GF3D<std::remove_reference<decltype(*hDD02GF)>::type,0,0,0> hDD02GF_(cctkGH,hDD02GF);
Loop::GF3D<std::remove_reference<decltype(*aDD02GF)>::type,0,0,0> aDD02GF_(cctkGH,aDD02GF);
Loop::GF3D<std::remove_reference<decltype(*betaSphorCartU1GF)>::type,0,0,0> betaSphorCartU1GF_(cctkGH,betaSphorCartU1GF);
Loop::GF3D<std::remove_reference<decltype(*BSphorCartU1GF)>::type,0,0,0> BSphorCartU1GF_(cctkGH,BSphorCartU1GF);
Loop::GF3D<std::remove_reference<decltype(*vetU1GF)>::type,0,0,0> vetU1GF_(cctkGH,vetU1GF);
Loop::GF3D<std::remove_reference<decltype(*betU1GF)>::type,0,0,0> betU1GF_(cctkGH,betU1GF);
Loop::GF3D<std::remove_reference<decltype(*gammaSphorCartDD11GF)>::type,0,0,0> gammaSphorCartDD11GF_(cctkGH,gammaSphorCartDD11GF);
Loop::GF3D<std::remove_reference<decltype(*KSphorCartDD11GF)>::type,0,0,0> KSphorCartDD11GF_(cctkGH,KSphorCartDD11GF);
Loop::GF3D<std::remove_reference<decltype(*hDD11GF)>::type,0,0,0> hDD11GF_(cctkGH,hDD11GF);
Loop::GF3D<std::remove_reference<decltype(*aDD11GF)>::type,0,0,0> aDD11GF_(cctkGH,aDD11GF);
Loop::GF3D<std::remove_reference<decltype(*gammaSphorCartDD12GF)>::type,0,0,0> gammaSphorCartDD12GF_(cctkGH,gammaSphorCartDD12GF);
Loop::GF3D<std::remove_reference<decltype(*KSphorCartDD12GF)>::type,0,0,0> KSphorCartDD12GF_(cctkGH,KSphorCartDD12GF);
Loop::GF3D<std::remove_reference<decltype(*hDD12GF)>::type,0,0,0> hDD12GF_(cctkGH,hDD12GF);
Loop::GF3D<std::remove_reference<decltype(*aDD12GF)>::type,0,0,0> aDD12GF_(cctkGH,aDD12GF);
Loop::GF3D<std::remove_reference<decltype(*betaSphorCartU2GF)>::type,0,0,0> betaSphorCartU2GF_(cctkGH,betaSphorCartU2GF);
Loop::GF3D<std::remove_reference<decltype(*BSphorCartU2GF)>::type,0,0,0> BSphorCartU2GF_(cctkGH,BSphorCartU2GF);
Loop::GF3D<std::remove_reference<decltype(*vetU2GF)>::type,0,0,0> vetU2GF_(cctkGH,vetU2GF);
Loop::GF3D<std::remove_reference<decltype(*betU2GF)>::type,0,0,0> betU2GF_(cctkGH,betU2GF);
Loop::GF3D<std::remove_reference<decltype(*gammaSphorCartDD22GF)>::type,0,0,0> gammaSphorCartDD22GF_(cctkGH,gammaSphorCartDD22GF);
Loop::GF3D<std::remove_reference<decltype(*KSphorCartDD22GF)>::type,0,0,0> KSphorCartDD22GF_(cctkGH,KSphorCartDD22GF);
Loop::GF3D<std::remove_reference<decltype(*hDD22GF)>::type,0,0,0> hDD22GF_(cctkGH,hDD22GF);
Loop::GF3D<std::remove_reference<decltype(*aDD22GF)>::type,0,0,0> aDD22GF_(cctkGH,aDD22GF);


Loop::loop_all<0,0,0>(cctkGH, [&](const Loop::PointDesc &p){
    const int i0 = p.i;
    const int i1 = p.j;
    const int i2 = p.k;
   /* 
    * NRPy+ Finite Difference Code Generation, Step 1 of 1: Read from main memory and compute finite difference stencils:
    */
   const double alphaSphorCart = alphaSphorCartGF_.ptr[alphaSphorCartGF_.offset(i0,i1,i2)];
   const double betaSphorCartU0 = betaSphorCartU0GF_.ptr[betaSphorCartU0GF_.offset(i0,i1,i2)];
   const double betaSphorCartU1 = betaSphorCartU1GF_.ptr[betaSphorCartU1GF_.offset(i0,i1,i2)];
   const double betaSphorCartU2 = betaSphorCartU2GF_.ptr[betaSphorCartU2GF_.offset(i0,i1,i2)];
   const double BSphorCartU0 = BSphorCartU0GF_.ptr[BSphorCartU0GF_.offset(i0,i1,i2)];
   const double BSphorCartU1 = BSphorCartU1GF_.ptr[BSphorCartU1GF_.offset(i0,i1,i2)];
   const double BSphorCartU2 = BSphorCartU2GF_.ptr[BSphorCartU2GF_.offset(i0,i1,i2)];
   const double gammaSphorCartDD00 = gammaSphorCartDD00GF_.ptr[gammaSphorCartDD00GF_.offset(i0,i1,i2)];
   const double gammaSphorCartDD01 = gammaSphorCartDD01GF_.ptr[gammaSphorCartDD01GF_.offset(i0,i1,i2)];
   const double gammaSphorCartDD02 = gammaSphorCartDD02GF_.ptr[gammaSphorCartDD02GF_.offset(i0,i1,i2)];
   const double gammaSphorCartDD11 = gammaSphorCartDD11GF_.ptr[gammaSphorCartDD11GF_.offset(i0,i1,i2)];
   const double gammaSphorCartDD12 = gammaSphorCartDD12GF_.ptr[gammaSphorCartDD12GF_.offset(i0,i1,i2)];
   const double gammaSphorCartDD22 = gammaSphorCartDD22GF_.ptr[gammaSphorCartDD22GF_.offset(i0,i1,i2)];
   const double KSphorCartDD00 = KSphorCartDD00GF_.ptr[KSphorCartDD00GF_.offset(i0,i1,i2)];
   const double KSphorCartDD01 = KSphorCartDD01GF_.ptr[KSphorCartDD01GF_.offset(i0,i1,i2)];
   const double KSphorCartDD02 = KSphorCartDD02GF_.ptr[KSphorCartDD02GF_.offset(i0,i1,i2)];
   const double KSphorCartDD11 = KSphorCartDD11GF_.ptr[KSphorCartDD11GF_.offset(i0,i1,i2)];
   const double KSphorCartDD12 = KSphorCartDD12GF_.ptr[KSphorCartDD12GF_.offset(i0,i1,i2)];
   const double KSphorCartDD22 = KSphorCartDD22GF_.ptr[KSphorCartDD22GF_.offset(i0,i1,i2)];
   /* 
    * NRPy+ Finite Difference Code Generation, Step 2 of 1: Evaluate SymPy expressions and write to main memory:
    */
   const double tmp0 = gammaSphorCartDD11*gammaSphorCartDD22;
   const double tmp1 = gammaSphorCartDD00*tmp0;
   const double tmp2 = gammaSphorCartDD02*gammaSphorCartDD12;
   const double tmp3 = 2*gammaSphorCartDD01*tmp2;
   const double tmp4 = ((gammaSphorCartDD12)*(gammaSphorCartDD12));
   const double tmp5 = gammaSphorCartDD00*tmp4;
   const double tmp6 = ((gammaSphorCartDD01)*(gammaSphorCartDD01));
   const double tmp7 = gammaSphorCartDD22*tmp6;
   const double tmp8 = ((gammaSphorCartDD02)*(gammaSphorCartDD02));
   const double tmp9 = gammaSphorCartDD11*tmp8;
   const double tmp10 = tmp1 + tmp3 - tmp5 - tmp7 - tmp9;
   const double tmp11 = (1.0/(tmp10));
   const double tmp12 = cbrt(tmp11);
   const double tmp13 = 2*tmp11;
   const double tmp14 = KSphorCartDD00*tmp11*(tmp0 - tmp4) + KSphorCartDD01*tmp13*(-gammaSphorCartDD01*gammaSphorCartDD22 + tmp2) + KSphorCartDD02*tmp13*(gammaSphorCartDD01*gammaSphorCartDD12 - gammaSphorCartDD02*gammaSphorCartDD11) + KSphorCartDD11*tmp11*(gammaSphorCartDD00*gammaSphorCartDD22 - tmp8) + KSphorCartDD12*tmp13*(-gammaSphorCartDD00*gammaSphorCartDD12 + gammaSphorCartDD01*gammaSphorCartDD02) + KSphorCartDD22*tmp11*(gammaSphorCartDD00*gammaSphorCartDD11 - tmp6);
   const double tmp15 = (1.0/3.0)*tmp14;
   hDD00GF_.ptr[hDD00GF_.offset(i0, i1, i2)] = gammaSphorCartDD00*tmp12 - 1;
   hDD01GF_.ptr[hDD01GF_.offset(i0, i1, i2)] = gammaSphorCartDD01*tmp12;
   hDD02GF_.ptr[hDD02GF_.offset(i0, i1, i2)] = gammaSphorCartDD02*tmp12;
   hDD11GF_.ptr[hDD11GF_.offset(i0, i1, i2)] = gammaSphorCartDD11*tmp12 - 1;
   hDD12GF_.ptr[hDD12GF_.offset(i0, i1, i2)] = gammaSphorCartDD12*tmp12;
   hDD22GF_.ptr[hDD22GF_.offset(i0, i1, i2)] = gammaSphorCartDD22*tmp12 - 1;
   aDD00GF_.ptr[aDD00GF_.offset(i0, i1, i2)] = tmp12*(KSphorCartDD00 - gammaSphorCartDD00*tmp15);
   aDD01GF_.ptr[aDD01GF_.offset(i0, i1, i2)] = tmp12*(KSphorCartDD01 - gammaSphorCartDD01*tmp15);
   aDD02GF_.ptr[aDD02GF_.offset(i0, i1, i2)] = tmp12*(KSphorCartDD02 - gammaSphorCartDD02*tmp15);
   aDD11GF_.ptr[aDD11GF_.offset(i0, i1, i2)] = tmp12*(KSphorCartDD11 - gammaSphorCartDD11*tmp15);
   aDD12GF_.ptr[aDD12GF_.offset(i0, i1, i2)] = tmp12*(KSphorCartDD12 - gammaSphorCartDD12*tmp15);
   aDD22GF_.ptr[aDD22GF_.offset(i0, i1, i2)] = tmp12*(KSphorCartDD22 - gammaSphorCartDD22*tmp15);
   trKGF_.ptr[trKGF_.offset(i0, i1, i2)] = tmp14;
   vetU0GF_.ptr[vetU0GF_.offset(i0, i1, i2)] = betaSphorCartU0;
   vetU1GF_.ptr[vetU1GF_.offset(i0, i1, i2)] = betaSphorCartU1;
   vetU2GF_.ptr[vetU2GF_.offset(i0, i1, i2)] = betaSphorCartU2;
   betU0GF_.ptr[betU0GF_.offset(i0, i1, i2)] = BSphorCartU0;
   betU1GF_.ptr[betU1GF_.offset(i0, i1, i2)] = BSphorCartU1;
   betU2GF_.ptr[betU2GF_.offset(i0, i1, i2)] = BSphorCartU2;
   alphaGF_.ptr[alphaGF_.offset(i0, i1, i2)] = alphaSphorCart;
   cfGF_.ptr[cfGF_.offset(i0, i1, i2)] = pow(tmp10/(tmp1*tmp11 + tmp11*tmp3 - tmp11*tmp5 - tmp11*tmp7 - tmp11*tmp9), -1.0/6.0);


});

}

void BaikalX_ADM_to_BSSN_lambdaU(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS_BaikalX_ADM_to_BSSN_lambdaU;
    DECLARE_CCTK_PARAMETERS;
    

Loop::GF3D<std::remove_reference<decltype(*lambdaU0GF)>::type,0,0,0> lambdaU0GF_(cctkGH,lambdaU0GF);
Loop::GF3D<std::remove_reference<decltype(*hDD00GF)>::type,0,0,0> hDD00GF_(cctkGH,hDD00GF);
Loop::GF3D<std::remove_reference<decltype(*hDD01GF)>::type,0,0,0> hDD01GF_(cctkGH,hDD01GF);
Loop::GF3D<std::remove_reference<decltype(*hDD02GF)>::type,0,0,0> hDD02GF_(cctkGH,hDD02GF);
Loop::GF3D<std::remove_reference<decltype(*lambdaU1GF)>::type,0,0,0> lambdaU1GF_(cctkGH,lambdaU1GF);
Loop::GF3D<std::remove_reference<decltype(*hDD11GF)>::type,0,0,0> hDD11GF_(cctkGH,hDD11GF);
Loop::GF3D<std::remove_reference<decltype(*hDD12GF)>::type,0,0,0> hDD12GF_(cctkGH,hDD12GF);
Loop::GF3D<std::remove_reference<decltype(*lambdaU2GF)>::type,0,0,0> lambdaU2GF_(cctkGH,lambdaU2GF);
Loop::GF3D<std::remove_reference<decltype(*hDD22GF)>::type,0,0,0> hDD22GF_(cctkGH,hDD22GF);


    const CCTK_REAL invdx0 = 1.0/CCTK_DELTA_SPACE(0);
    const CCTK_REAL invdx1 = 1.0/CCTK_DELTA_SPACE(1);
    const CCTK_REAL invdx2 = 1.0/CCTK_DELTA_SPACE(2);
Loop::loop_int<0,0,0>(cctkGH, [&](const Loop::PointDesc &p){
    const int i0 = p.i;
    const int i1 = p.j;
    const int i2 = p.k;
   /* 
    * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
    */
   const double hDD00_i0_i1_i2m1 = hDD00GF_.ptr[hDD00GF_.offset(i0,i1,i2-1)];
   const double hDD00_i0_i1m1_i2 = hDD00GF_.ptr[hDD00GF_.offset(i0,i1-1,i2)];
   const double hDD00_i0m1_i1_i2 = hDD00GF_.ptr[hDD00GF_.offset(i0-1,i1,i2)];
   const double hDD00 = hDD00GF_.ptr[hDD00GF_.offset(i0,i1,i2)];
   const double hDD00_i0p1_i1_i2 = hDD00GF_.ptr[hDD00GF_.offset(i0+1,i1,i2)];
   const double hDD00_i0_i1p1_i2 = hDD00GF_.ptr[hDD00GF_.offset(i0,i1+1,i2)];
   const double hDD00_i0_i1_i2p1 = hDD00GF_.ptr[hDD00GF_.offset(i0,i1,i2+1)];
   const double hDD01_i0_i1_i2m1 = hDD01GF_.ptr[hDD01GF_.offset(i0,i1,i2-1)];
   const double hDD01_i0_i1m1_i2 = hDD01GF_.ptr[hDD01GF_.offset(i0,i1-1,i2)];
   const double hDD01_i0m1_i1_i2 = hDD01GF_.ptr[hDD01GF_.offset(i0-1,i1,i2)];
   const double hDD01 = hDD01GF_.ptr[hDD01GF_.offset(i0,i1,i2)];
   const double hDD01_i0p1_i1_i2 = hDD01GF_.ptr[hDD01GF_.offset(i0+1,i1,i2)];
   const double hDD01_i0_i1p1_i2 = hDD01GF_.ptr[hDD01GF_.offset(i0,i1+1,i2)];
   const double hDD01_i0_i1_i2p1 = hDD01GF_.ptr[hDD01GF_.offset(i0,i1,i2+1)];
   const double hDD02_i0_i1_i2m1 = hDD02GF_.ptr[hDD02GF_.offset(i0,i1,i2-1)];
   const double hDD02_i0_i1m1_i2 = hDD02GF_.ptr[hDD02GF_.offset(i0,i1-1,i2)];
   const double hDD02_i0m1_i1_i2 = hDD02GF_.ptr[hDD02GF_.offset(i0-1,i1,i2)];
   const double hDD02 = hDD02GF_.ptr[hDD02GF_.offset(i0,i1,i2)];
   const double hDD02_i0p1_i1_i2 = hDD02GF_.ptr[hDD02GF_.offset(i0+1,i1,i2)];
   const double hDD02_i0_i1p1_i2 = hDD02GF_.ptr[hDD02GF_.offset(i0,i1+1,i2)];
   const double hDD02_i0_i1_i2p1 = hDD02GF_.ptr[hDD02GF_.offset(i0,i1,i2+1)];
   const double hDD11_i0_i1_i2m1 = hDD11GF_.ptr[hDD11GF_.offset(i0,i1,i2-1)];
   const double hDD11_i0_i1m1_i2 = hDD11GF_.ptr[hDD11GF_.offset(i0,i1-1,i2)];
   const double hDD11_i0m1_i1_i2 = hDD11GF_.ptr[hDD11GF_.offset(i0-1,i1,i2)];
   const double hDD11 = hDD11GF_.ptr[hDD11GF_.offset(i0,i1,i2)];
   const double hDD11_i0p1_i1_i2 = hDD11GF_.ptr[hDD11GF_.offset(i0+1,i1,i2)];
   const double hDD11_i0_i1p1_i2 = hDD11GF_.ptr[hDD11GF_.offset(i0,i1+1,i2)];
   const double hDD11_i0_i1_i2p1 = hDD11GF_.ptr[hDD11GF_.offset(i0,i1,i2+1)];
   const double hDD12_i0_i1_i2m1 = hDD12GF_.ptr[hDD12GF_.offset(i0,i1,i2-1)];
   const double hDD12_i0_i1m1_i2 = hDD12GF_.ptr[hDD12GF_.offset(i0,i1-1,i2)];
   const double hDD12_i0m1_i1_i2 = hDD12GF_.ptr[hDD12GF_.offset(i0-1,i1,i2)];
   const double hDD12 = hDD12GF_.ptr[hDD12GF_.offset(i0,i1,i2)];
   const double hDD12_i0p1_i1_i2 = hDD12GF_.ptr[hDD12GF_.offset(i0+1,i1,i2)];
   const double hDD12_i0_i1p1_i2 = hDD12GF_.ptr[hDD12GF_.offset(i0,i1+1,i2)];
   const double hDD12_i0_i1_i2p1 = hDD12GF_.ptr[hDD12GF_.offset(i0,i1,i2+1)];
   const double hDD22_i0_i1_i2m1 = hDD22GF_.ptr[hDD22GF_.offset(i0,i1,i2-1)];
   const double hDD22_i0_i1m1_i2 = hDD22GF_.ptr[hDD22GF_.offset(i0,i1-1,i2)];
   const double hDD22_i0m1_i1_i2 = hDD22GF_.ptr[hDD22GF_.offset(i0-1,i1,i2)];
   const double hDD22 = hDD22GF_.ptr[hDD22GF_.offset(i0,i1,i2)];
   const double hDD22_i0p1_i1_i2 = hDD22GF_.ptr[hDD22GF_.offset(i0+1,i1,i2)];
   const double hDD22_i0_i1p1_i2 = hDD22GF_.ptr[hDD22GF_.offset(i0,i1+1,i2)];
   const double hDD22_i0_i1_i2p1 = hDD22GF_.ptr[hDD22GF_.offset(i0,i1,i2+1)];
         const double hDD_dD000 = invdx0*(-1.0/2.0*hDD00_i0m1_i1_i2 + (1.0/2.0)*hDD00_i0p1_i1_i2);
         const double hDD_dD001 = invdx1*(-1.0/2.0*hDD00_i0_i1m1_i2 + (1.0/2.0)*hDD00_i0_i1p1_i2);
         const double hDD_dD002 = invdx2*(-1.0/2.0*hDD00_i0_i1_i2m1 + (1.0/2.0)*hDD00_i0_i1_i2p1);
         const double hDD_dD010 = invdx0*(-1.0/2.0*hDD01_i0m1_i1_i2 + (1.0/2.0)*hDD01_i0p1_i1_i2);
         const double hDD_dD011 = invdx1*(-1.0/2.0*hDD01_i0_i1m1_i2 + (1.0/2.0)*hDD01_i0_i1p1_i2);
         const double hDD_dD012 = invdx2*(-1.0/2.0*hDD01_i0_i1_i2m1 + (1.0/2.0)*hDD01_i0_i1_i2p1);
         const double hDD_dD020 = invdx0*(-1.0/2.0*hDD02_i0m1_i1_i2 + (1.0/2.0)*hDD02_i0p1_i1_i2);
         const double hDD_dD021 = invdx1*(-1.0/2.0*hDD02_i0_i1m1_i2 + (1.0/2.0)*hDD02_i0_i1p1_i2);
         const double hDD_dD022 = invdx2*(-1.0/2.0*hDD02_i0_i1_i2m1 + (1.0/2.0)*hDD02_i0_i1_i2p1);
         const double hDD_dD110 = invdx0*(-1.0/2.0*hDD11_i0m1_i1_i2 + (1.0/2.0)*hDD11_i0p1_i1_i2);
         const double hDD_dD111 = invdx1*(-1.0/2.0*hDD11_i0_i1m1_i2 + (1.0/2.0)*hDD11_i0_i1p1_i2);
         const double hDD_dD112 = invdx2*(-1.0/2.0*hDD11_i0_i1_i2m1 + (1.0/2.0)*hDD11_i0_i1_i2p1);
         const double hDD_dD120 = invdx0*(-1.0/2.0*hDD12_i0m1_i1_i2 + (1.0/2.0)*hDD12_i0p1_i1_i2);
         const double hDD_dD121 = invdx1*(-1.0/2.0*hDD12_i0_i1m1_i2 + (1.0/2.0)*hDD12_i0_i1p1_i2);
         const double hDD_dD122 = invdx2*(-1.0/2.0*hDD12_i0_i1_i2m1 + (1.0/2.0)*hDD12_i0_i1_i2p1);
         const double hDD_dD220 = invdx0*(-1.0/2.0*hDD22_i0m1_i1_i2 + (1.0/2.0)*hDD22_i0p1_i1_i2);
         const double hDD_dD221 = invdx1*(-1.0/2.0*hDD22_i0_i1m1_i2 + (1.0/2.0)*hDD22_i0_i1p1_i2);
         const double hDD_dD222 = invdx2*(-1.0/2.0*hDD22_i0_i1_i2m1 + (1.0/2.0)*hDD22_i0_i1_i2p1);
   /* 
    * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
    */
   const double tmp0 = hDD22 + 1;
   const double tmp1 = -hDD01*tmp0 + hDD02*hDD12;
   const double tmp2 = hDD01*hDD12;
   const double tmp3 = ((hDD01)*(hDD01));
   const double tmp4 = ((hDD02)*(hDD02));
   const double tmp5 = hDD11 + 1;
   const double tmp6 = ((hDD12)*(hDD12));
   const double tmp7 = hDD00 + 1;
   const double tmp8 = tmp5*tmp7;
   const double tmp9 = (1.0/(2*hDD02*tmp2 - tmp0*tmp3 + tmp0*tmp8 - tmp4*tmp5 - tmp6*tmp7));
   const double tmp10 = (1.0/2.0)*tmp9;
   const double tmp11 = tmp1*tmp10;
   const double tmp12 = -hDD02*tmp5 + tmp2;
   const double tmp13 = tmp10*tmp12;
   const double tmp14 = hDD_dD012 + hDD_dD021 - hDD_dD120;
   const double tmp15 = tmp9*(tmp0*tmp5 - tmp6);
   const double tmp16 = (1.0/2.0)*tmp15;
   const double tmp17 = hDD01*hDD02 - hDD12*tmp7;
   const double tmp18 = 2*tmp9;
   const double tmp19 = tmp17*tmp18;
   const double tmp20 = hDD_dD012 - hDD_dD021 + hDD_dD120;
   const double tmp21 = tmp12*tmp18;
   const double tmp22 = -hDD_dD012 + hDD_dD021 + hDD_dD120;
   const double tmp23 = tmp1*tmp18;
   const double tmp24 = 2*hDD_dD122 - hDD_dD221;
   const double tmp25 = 2*hDD_dD022 - hDD_dD220;
   const double tmp26 = tmp9*(-tmp3 + tmp8);
   const double tmp27 = -hDD_dD112 + 2*hDD_dD121;
   const double tmp28 = 2*hDD_dD011 - hDD_dD110;
   const double tmp29 = tmp9*(tmp0*tmp7 - tmp4);
   const double tmp30 = -hDD_dD001 + 2*hDD_dD010;
   const double tmp31 = -hDD_dD002 + 2*hDD_dD020;
   const double tmp32 = tmp10*tmp17;
   const double tmp33 = (1.0/2.0)*tmp29;
   const double tmp34 = (1.0/2.0)*tmp26;
   lambdaU0GF_.ptr[lambdaU0GF_.offset(i0, i1, i2)] = tmp15*(hDD_dD000*tmp16 + tmp11*tmp30 + tmp13*tmp31) + tmp19*(hDD_dD112*tmp11 + hDD_dD221*tmp13 + tmp14*tmp16) + tmp21*(hDD_dD002*tmp16 + hDD_dD220*tmp13 + tmp11*tmp20) + tmp23*(hDD_dD001*tmp16 + hDD_dD110*tmp11 + tmp13*tmp22) + tmp26*(hDD_dD222*tmp13 + tmp11*tmp24 + tmp16*tmp25) + tmp29*(hDD_dD111*tmp11 + tmp13*tmp27 + tmp16*tmp28);
   lambdaU1GF_.ptr[lambdaU1GF_.offset(i0, i1, i2)] = tmp15*(hDD_dD000*tmp11 + tmp30*tmp33 + tmp31*tmp32) + tmp19*(hDD_dD112*tmp33 + hDD_dD221*tmp32 + tmp11*tmp14) + tmp21*(hDD_dD002*tmp11 + hDD_dD220*tmp32 + tmp20*tmp33) + tmp23*(hDD_dD001*tmp11 + hDD_dD110*tmp33 + tmp22*tmp32) + tmp26*(hDD_dD222*tmp32 + tmp11*tmp25 + tmp24*tmp33) + tmp29*(hDD_dD111*tmp33 + tmp11*tmp28 + tmp27*tmp32);
   lambdaU2GF_.ptr[lambdaU2GF_.offset(i0, i1, i2)] = tmp15*(hDD_dD000*tmp13 + tmp30*tmp32 + tmp31*tmp34) + tmp19*(hDD_dD112*tmp32 + hDD_dD221*tmp34 + tmp13*tmp14) + tmp21*(hDD_dD002*tmp13 + hDD_dD220*tmp34 + tmp20*tmp32) + tmp23*(hDD_dD001*tmp13 + hDD_dD110*tmp32 + tmp22*tmp34) + tmp26*(hDD_dD222*tmp34 + tmp13*tmp25 + tmp24*tmp32) + tmp29*(hDD_dD111*tmp32 + tmp13*tmp28 + tmp27*tmp34);


});
Loop::loop_bnd<0,0,0>(cctkGH, [&](const Loop::PointDesc &p){
    const int i0 = p.i;
    const int i1 = p.j;
    const int i2 = p.k;

    lambdaU0GF_(i0,i1,i2) = 0.;
    lambdaU1GF_(i0,i1,i2) = 0.;
    lambdaU2GF_(i0,i1,i2) = 0.;
    
});

    //ExtrapolateGammas(cctkGH,lambdaU0GF);
    //ExtrapolateGammas(cctkGH,lambdaU1GF);
    //ExtrapolateGammas(cctkGH,lambdaU2GF);
}
