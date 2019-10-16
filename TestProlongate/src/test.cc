#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "loop.hxx"

#include <cstddef>
#include <cstdio>

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

    data[idx] = xL;
  });
}

void TestProlongate_Sync(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  return; // do nothing
}

void TestProlongate_Regrid(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(cctk_iteration < 1)
    return;

  printf("Setting grid\n");

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    if(p.x > -0.5 && p.x < 0.5 &&
       p.y > -0.5 && p.y < 0.5 &&
       p.z > -0.5 && p.z < 0.5) {
      fprintf(stderr,"Setting error %g %g %g\n", p.x, p.y, p.z);
      regrid_error[p.idx] = 1e3;
    } else {
      regrid_error[p.idx] = 0.;
    }
  });
  fprintf(stderr, "done\n");
}
