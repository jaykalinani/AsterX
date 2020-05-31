
#include <loop.hxx>

#include <math.h>

#include "cctk.h"
#include "cctk_Arguments_Checked.h"
#include "cctk_Parameters.h"

void BaikalX_BSSN_to_ADM(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS_BaikalX_BSSN_to_ADM;
    DECLARE_CCTK_PARAMETERS;
    
Loop::GF3D<std::remove_reference<decltype(*trKGF)>::type,0,0,0> trKGF_(cctkGH,trKGF);
Loop::GF3D<std::remove_reference<decltype(*alphaGF)>::type,0,0,0> alphaGF_(cctkGH,alphaGF);
Loop::GF3D<std::remove_reference<decltype(*alp)>::type,0,0,0> alp_(cctkGH,alp);
Loop::GF3D<std::remove_reference<decltype(*dtalp)>::type,0,0,0> dtalp_(cctkGH,dtalp);
Loop::GF3D<std::remove_reference<decltype(*dt2alp)>::type,0,0,0> dt2alp_(cctkGH,dt2alp);
Loop::GF3D<std::remove_reference<decltype(*cfGF)>::type,0,0,0> cfGF_(cctkGH,cfGF);
Loop::GF3D<std::remove_reference<decltype(*betax)>::type,0,0,0> betax_(cctkGH,betax);
Loop::GF3D<std::remove_reference<decltype(*dtbetax)>::type,0,0,0> dtbetax_(cctkGH,dtbetax);
Loop::GF3D<std::remove_reference<decltype(*dt2betax)>::type,0,0,0> dt2betax_(cctkGH,dt2betax);
Loop::GF3D<std::remove_reference<decltype(*vetU0GF)>::type,0,0,0> vetU0GF_(cctkGH,vetU0GF);
Loop::GF3D<std::remove_reference<decltype(*betU0GF)>::type,0,0,0> betU0GF_(cctkGH,betU0GF);
Loop::GF3D<std::remove_reference<decltype(*gxx)>::type,0,0,0> gxx_(cctkGH,gxx);
Loop::GF3D<std::remove_reference<decltype(*kxx)>::type,0,0,0> kxx_(cctkGH,kxx);
Loop::GF3D<std::remove_reference<decltype(*dtkxx)>::type,0,0,0> dtkxx_(cctkGH,dtkxx);
Loop::GF3D<std::remove_reference<decltype(*hDD00GF)>::type,0,0,0> hDD00GF_(cctkGH,hDD00GF);
Loop::GF3D<std::remove_reference<decltype(*aDD00GF)>::type,0,0,0> aDD00GF_(cctkGH,aDD00GF);
Loop::GF3D<std::remove_reference<decltype(*gxy)>::type,0,0,0> gxy_(cctkGH,gxy);
Loop::GF3D<std::remove_reference<decltype(*kxy)>::type,0,0,0> kxy_(cctkGH,kxy);
Loop::GF3D<std::remove_reference<decltype(*dtkxy)>::type,0,0,0> dtkxy_(cctkGH,dtkxy);
Loop::GF3D<std::remove_reference<decltype(*hDD01GF)>::type,0,0,0> hDD01GF_(cctkGH,hDD01GF);
Loop::GF3D<std::remove_reference<decltype(*aDD01GF)>::type,0,0,0> aDD01GF_(cctkGH,aDD01GF);
Loop::GF3D<std::remove_reference<decltype(*gxz)>::type,0,0,0> gxz_(cctkGH,gxz);
Loop::GF3D<std::remove_reference<decltype(*kxz)>::type,0,0,0> kxz_(cctkGH,kxz);
Loop::GF3D<std::remove_reference<decltype(*dtkxz)>::type,0,0,0> dtkxz_(cctkGH,dtkxz);
Loop::GF3D<std::remove_reference<decltype(*hDD02GF)>::type,0,0,0> hDD02GF_(cctkGH,hDD02GF);
Loop::GF3D<std::remove_reference<decltype(*aDD02GF)>::type,0,0,0> aDD02GF_(cctkGH,aDD02GF);
Loop::GF3D<std::remove_reference<decltype(*betay)>::type,0,0,0> betay_(cctkGH,betay);
Loop::GF3D<std::remove_reference<decltype(*dtbetay)>::type,0,0,0> dtbetay_(cctkGH,dtbetay);
Loop::GF3D<std::remove_reference<decltype(*dt2betay)>::type,0,0,0> dt2betay_(cctkGH,dt2betay);
Loop::GF3D<std::remove_reference<decltype(*vetU1GF)>::type,0,0,0> vetU1GF_(cctkGH,vetU1GF);
Loop::GF3D<std::remove_reference<decltype(*betU1GF)>::type,0,0,0> betU1GF_(cctkGH,betU1GF);
Loop::GF3D<std::remove_reference<decltype(*gyy)>::type,0,0,0> gyy_(cctkGH,gyy);
Loop::GF3D<std::remove_reference<decltype(*kyy)>::type,0,0,0> kyy_(cctkGH,kyy);
Loop::GF3D<std::remove_reference<decltype(*dtkyy)>::type,0,0,0> dtkyy_(cctkGH,dtkyy);
Loop::GF3D<std::remove_reference<decltype(*hDD11GF)>::type,0,0,0> hDD11GF_(cctkGH,hDD11GF);
Loop::GF3D<std::remove_reference<decltype(*aDD11GF)>::type,0,0,0> aDD11GF_(cctkGH,aDD11GF);
Loop::GF3D<std::remove_reference<decltype(*gyz)>::type,0,0,0> gyz_(cctkGH,gyz);
Loop::GF3D<std::remove_reference<decltype(*kyz)>::type,0,0,0> kyz_(cctkGH,kyz);
Loop::GF3D<std::remove_reference<decltype(*dtkyz)>::type,0,0,0> dtkyz_(cctkGH,dtkyz);
Loop::GF3D<std::remove_reference<decltype(*hDD12GF)>::type,0,0,0> hDD12GF_(cctkGH,hDD12GF);
Loop::GF3D<std::remove_reference<decltype(*aDD12GF)>::type,0,0,0> aDD12GF_(cctkGH,aDD12GF);
Loop::GF3D<std::remove_reference<decltype(*betaz)>::type,0,0,0> betaz_(cctkGH,betaz);
Loop::GF3D<std::remove_reference<decltype(*dtbetaz)>::type,0,0,0> dtbetaz_(cctkGH,dtbetaz);
Loop::GF3D<std::remove_reference<decltype(*dt2betaz)>::type,0,0,0> dt2betaz_(cctkGH,dt2betaz);
Loop::GF3D<std::remove_reference<decltype(*vetU2GF)>::type,0,0,0> vetU2GF_(cctkGH,vetU2GF);
Loop::GF3D<std::remove_reference<decltype(*betU2GF)>::type,0,0,0> betU2GF_(cctkGH,betU2GF);
Loop::GF3D<std::remove_reference<decltype(*gzz)>::type,0,0,0> gzz_(cctkGH,gzz);
Loop::GF3D<std::remove_reference<decltype(*kzz)>::type,0,0,0> kzz_(cctkGH,kzz);
Loop::GF3D<std::remove_reference<decltype(*dtkzz)>::type,0,0,0> dtkzz_(cctkGH,dtkzz);
Loop::GF3D<std::remove_reference<decltype(*hDD22GF)>::type,0,0,0> hDD22GF_(cctkGH,hDD22GF);
Loop::GF3D<std::remove_reference<decltype(*aDD22GF)>::type,0,0,0> aDD22GF_(cctkGH,aDD22GF);

Loop::loop_all<0,0,0>(cctkGH, [&](const Loop::PointDesc &p){
    const int i0 = p.i;
    const int i1 = p.j;
    const int i2 = p.k;
   /* 
    * NRPy+ Finite Difference Code Generation, Step 1 of 1: Read from main memory and compute finite difference stencils:
    */
   const double hDD00 = hDD00GF_.ptr[hDD00GF_.offset(i0,i1,i2)];
   const double hDD01 = hDD01GF_.ptr[hDD01GF_.offset(i0,i1,i2)];
   const double hDD02 = hDD02GF_.ptr[hDD02GF_.offset(i0,i1,i2)];
   const double hDD11 = hDD11GF_.ptr[hDD11GF_.offset(i0,i1,i2)];
   const double hDD12 = hDD12GF_.ptr[hDD12GF_.offset(i0,i1,i2)];
   const double hDD22 = hDD22GF_.ptr[hDD22GF_.offset(i0,i1,i2)];
   const double aDD00 = aDD00GF_.ptr[aDD00GF_.offset(i0,i1,i2)];
   const double aDD01 = aDD01GF_.ptr[aDD01GF_.offset(i0,i1,i2)];
   const double aDD02 = aDD02GF_.ptr[aDD02GF_.offset(i0,i1,i2)];
   const double aDD11 = aDD11GF_.ptr[aDD11GF_.offset(i0,i1,i2)];
   const double aDD12 = aDD12GF_.ptr[aDD12GF_.offset(i0,i1,i2)];
   const double aDD22 = aDD22GF_.ptr[aDD22GF_.offset(i0,i1,i2)];
   const double vetU0 = vetU0GF_.ptr[vetU0GF_.offset(i0,i1,i2)];
   const double vetU1 = vetU1GF_.ptr[vetU1GF_.offset(i0,i1,i2)];
   const double vetU2 = vetU2GF_.ptr[vetU2GF_.offset(i0,i1,i2)];
   const double betU0 = betU0GF_.ptr[betU0GF_.offset(i0,i1,i2)];
   const double betU1 = betU1GF_.ptr[betU1GF_.offset(i0,i1,i2)];
   const double betU2 = betU2GF_.ptr[betU2GF_.offset(i0,i1,i2)];
   const double trK = trKGF_.ptr[trKGF_.offset(i0,i1,i2)];
   const double cf = cfGF_.ptr[cfGF_.offset(i0,i1,i2)];
   const double alpha = alphaGF_.ptr[alphaGF_.offset(i0,i1,i2)];
   /* 
    * NRPy+ Finite Difference Code Generation, Step 2 of 1: Evaluate SymPy expressions and write to main memory:
    */
   const double tmp0 = (1.0/((cf)*(cf)));
   const double tmp1 = tmp0*(hDD00 + 1);
   const double tmp2 = hDD01*tmp0;
   const double tmp3 = hDD02*tmp0;
   const double tmp4 = tmp0*(hDD11 + 1);
   const double tmp5 = hDD12*tmp0;
   const double tmp6 = tmp0*(hDD22 + 1);
   const double tmp7 = (1.0/3.0)*trK;
   gxx_(p.I) = tmp1;
   gxy_(p.I) = tmp2;
   gxz_(p.I) = tmp3;
   gyy_(p.I) = tmp4;
   gyz_(p.I) = tmp5;
   gzz_(p.I) = tmp6;
   kxx_(p.I) = aDD00*tmp0 + tmp1*tmp7;
   kxy_(p.I) = aDD01*tmp0 + tmp2*tmp7;
   kxz_(p.I) = aDD02*tmp0 + tmp3*tmp7;
   kyy_(p.I) = aDD11*tmp0 + tmp4*tmp7;
   kyz_(p.I) = aDD12*tmp0 + tmp5*tmp7;
   kzz_(p.I) = aDD22*tmp0 + tmp6*tmp7;
   dtkxx_(p.I) = 0;
   dtkxy_(p.I) = 0;
   dtkxz_(p.I) = 0;
   dtkyy_(p.I) = 0;
   dtkyz_(p.I) = 0;
   dtkzz_(p.I) = 0;
   alp_(p.I) = alpha;
   dtalp_(p.I) = 0;
   dt2alp_(p.I) = 0;
   betax_(p.I) = vetU0;
   betay_(p.I) = vetU1;
   betaz_(p.I) = vetU2;
   dtbetax_(p.I) = betU0;
   dtbetay_(p.I) = betU1;
   dtbetaz_(p.I) = betU2;
   dt2betax_(p.I) = 0;
   dt2betay_(p.I) = 0;
   dt2betaz_(p.I) = 0;


});
}
