/*
Enforce det(gammabar) = det(gammahat) constraint.
 */
void enforce_detgammabar_constraint(const paramstruct *restrict params, REAL *restrict xx[3], REAL *restrict in_gfs) {
#include "set_Cparameters.h"

#pragma omp parallel for
    for(int i2=0; i2<Nxx_plus_2NGHOSTS2; i2++) {
        const REAL xx2 = xx[2][i2];
        for(int i1=0; i1<Nxx_plus_2NGHOSTS1; i1++) {
            const REAL xx1 = xx[1][i1];
            for(int i0=0; i0<Nxx_plus_2NGHOSTS0; i0++) {
                const REAL xx0 = xx[0][i0];
                   /* 
                    * NRPy+ Finite Difference Code Generation, Step 1 of 1: Read from main memory and compute finite difference stencils:
                    */
                   const double hDD00 = in_gfs[IDX4S(HDD00GF, i0,i1,i2)];
                   const double hDD01 = in_gfs[IDX4S(HDD01GF, i0,i1,i2)];
                   const double hDD02 = in_gfs[IDX4S(HDD02GF, i0,i1,i2)];
                   const double hDD11 = in_gfs[IDX4S(HDD11GF, i0,i1,i2)];
                   const double hDD12 = in_gfs[IDX4S(HDD12GF, i0,i1,i2)];
                   const double hDD22 = in_gfs[IDX4S(HDD22GF, i0,i1,i2)];
                   /* 
                    * NRPy+ Finite Difference Code Generation, Step 2 of 1: Evaluate SymPy expressions and write to main memory:
                    */
                   const double tmp0 = hDD00 + 1;
                   const double tmp1 = ((sin(xx1))*(sin(xx1)));
                   const double tmp2 = tmp1*((xx0)*(xx0)*(xx0)*(xx0));
                   const double tmp3 = ((xx0)*(xx0));
                   const double tmp4 = hDD11*tmp3 + tmp3;
                   const double tmp5 = tmp1*tmp3;
                   const double tmp6 = hDD22*tmp5 + tmp5;
                   const double tmp7 = cbrt(fabs(tmp2)/(-((hDD01)*(hDD01))*tmp3*tmp6 + 2*hDD01*hDD02*hDD12*tmp2 - ((hDD02)*(hDD02))*tmp4*tmp5 - ((hDD12)*(hDD12))*tmp0*tmp2 + tmp0*tmp4*tmp6));
                   in_gfs[IDX4S(HDD00GF, i0, i1, i2)] = tmp0*tmp7 - 1;
                   in_gfs[IDX4S(HDD01GF, i0, i1, i2)] = hDD01*tmp7;
                   in_gfs[IDX4S(HDD02GF, i0, i1, i2)] = hDD02*tmp7;
                   in_gfs[IDX4S(HDD11GF, i0, i1, i2)] = tmp7*(hDD11 + 1) - 1;
                   in_gfs[IDX4S(HDD12GF, i0, i1, i2)] = hDD12*tmp7;
                   in_gfs[IDX4S(HDD22GF, i0, i1, i2)] = tmp7*(hDD22 + 1) - 1;
                
                
            } // END LOOP: for(int i0=0; i0<Nxx_plus_2NGHOSTS0; i0++)
        } // END LOOP: for(int i1=0; i1<Nxx_plus_2NGHOSTS1; i1++)
    } // END LOOP: for(int i2=0; i2<Nxx_plus_2NGHOSTS2; i2++)
}
