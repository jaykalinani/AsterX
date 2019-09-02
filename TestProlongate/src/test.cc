#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "loop.hxx"

#include <stddef.h>

void TestProlongate_Test(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  CCTK_VINFO("Initializing data of size %d %d %d", cctk_lsh[0], cctk_lsh[1], cctk_lsh[2]);
#pragma omp parallel
  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) 
  {
    ptrdiff_t idx = p.idx;
    CCTK_REAL xL = p.x;
    CCTK_REAL yL = p.y;
    CCTK_REAL zL = p.z;

    printf("x[%d]: %g\n", idx, xL);
    data[idx] = xL;
  });
}

void TestProlongate_Sync(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  return; // do nothing
}
