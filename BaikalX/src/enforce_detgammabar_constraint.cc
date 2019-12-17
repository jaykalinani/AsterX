
#include <math.h>

#include <loop.hxx>

#include "cctk.h"
#include "cctk_Arguments_Checked.h"
#include "cctk_Parameters.h"

void enforce_detgammabar_constraint(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS_enforce_detgammabar_constraint;
    DECLARE_CCTK_PARAMETERS;
Loop::GF3D<CCTK_REAL,0,0,0> hDD00GF_(cctkGH,hDD00GF);
Loop::GF3D<CCTK_REAL,0,0,0> hDD01GF_(cctkGH,hDD01GF);
Loop::GF3D<CCTK_REAL,0,0,0> hDD02GF_(cctkGH,hDD02GF);
Loop::GF3D<CCTK_REAL,0,0,0> hDD11GF_(cctkGH,hDD11GF);
Loop::GF3D<CCTK_REAL,0,0,0> hDD12GF_(cctkGH,hDD12GF);
Loop::GF3D<CCTK_REAL,0,0,0> hDD22GF_(cctkGH,hDD22GF);


#include "enforce_detgammabar_constraint.h"
}
