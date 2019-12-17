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
/* 
 * NRPy+ Finite Difference Code Generation, Step 2 of 1: Evaluate SymPy expressions and write to main memory:
 */
const double tmp0 = hDD00 + 1;
const double tmp1 = hDD22 + 1;
const double tmp2 = hDD11 + 1;
const double tmp3 = cbrt(fabs(1)/(-((hDD01)*(hDD01))*tmp1 + 2*hDD01*hDD02*hDD12 - ((hDD02)*(hDD02))*tmp2 - ((hDD12)*(hDD12))*tmp0 + tmp0*tmp1*tmp2));
hDD00GF_.ptr[hDD00GF_.offset(i0, i1, i2)] = tmp0*tmp3 - 1;
hDD01GF_.ptr[hDD01GF_.offset(i0, i1, i2)] = hDD01*tmp3;
hDD02GF_.ptr[hDD02GF_.offset(i0, i1, i2)] = hDD02*tmp3;
hDD11GF_.ptr[hDD11GF_.offset(i0, i1, i2)] = tmp2*tmp3 - 1;
hDD12GF_.ptr[hDD12GF_.offset(i0, i1, i2)] = hDD12*tmp3;
hDD22GF_.ptr[hDD22GF_.offset(i0, i1, i2)] = tmp1*tmp3 - 1;


});
