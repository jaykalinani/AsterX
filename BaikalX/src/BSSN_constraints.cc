
#include <math.h>

#include <loop.hxx>

#include "cctk.h"
#include "cctk_Arguments_Checked.h"
#include "cctk_Parameters.h"

void BaikalX_BSSN_constraints(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS_BaikalX_BSSN_constraints;
    DECLARE_CCTK_PARAMETERS;

    const CCTK_REAL invdx0 = 1.0/CCTK_DELTA_SPACE(0);
    const CCTK_REAL invdx1 = 1.0/CCTK_DELTA_SPACE(1);
    const CCTK_REAL invdx2 = 1.0/CCTK_DELTA_SPACE(2);
Loop::GF3D<CCTK_REAL,0,0,0> trKGF_(cctkGH,trKGF);
Loop::GF3D<CCTK_REAL,0,0,0> cfGF_(cctkGH,cfGF);
Loop::GF3D<CCTK_REAL,0,0,0> HGF_(cctkGH,HGF);
Loop::GF3D<CCTK_REAL,0,0,0> lambdaU0GF_(cctkGH,lambdaU0GF);
Loop::GF3D<CCTK_REAL,0,0,0> MU0GF_(cctkGH,MU0GF);
Loop::GF3D<CCTK_REAL,0,0,0> hDD00GF_(cctkGH,hDD00GF);
Loop::GF3D<CCTK_REAL,0,0,0> aDD00GF_(cctkGH,aDD00GF);
Loop::GF3D<CCTK_REAL,0,0,0> hDD01GF_(cctkGH,hDD01GF);
Loop::GF3D<CCTK_REAL,0,0,0> aDD01GF_(cctkGH,aDD01GF);
Loop::GF3D<CCTK_REAL,0,0,0> hDD02GF_(cctkGH,hDD02GF);
Loop::GF3D<CCTK_REAL,0,0,0> aDD02GF_(cctkGH,aDD02GF);
Loop::GF3D<CCTK_REAL,0,0,0> lambdaU1GF_(cctkGH,lambdaU1GF);
Loop::GF3D<CCTK_REAL,0,0,0> MU1GF_(cctkGH,MU1GF);
Loop::GF3D<CCTK_REAL,0,0,0> hDD11GF_(cctkGH,hDD11GF);
Loop::GF3D<CCTK_REAL,0,0,0> aDD11GF_(cctkGH,aDD11GF);
Loop::GF3D<CCTK_REAL,0,0,0> hDD12GF_(cctkGH,hDD12GF);
Loop::GF3D<CCTK_REAL,0,0,0> aDD12GF_(cctkGH,aDD12GF);
Loop::GF3D<CCTK_REAL,0,0,0> lambdaU2GF_(cctkGH,lambdaU2GF);
Loop::GF3D<CCTK_REAL,0,0,0> MU2GF_(cctkGH,MU2GF);
Loop::GF3D<CCTK_REAL,0,0,0> hDD22GF_(cctkGH,hDD22GF);
Loop::GF3D<CCTK_REAL,0,0,0> aDD22GF_(cctkGH,aDD22GF);


#include "BSSN_constraints.h"
}
