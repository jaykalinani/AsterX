#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <cassert>
#include <cmath>
#include <memory>
#include <iostream>
#include <string>
#include <utility>
#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <tuple>
#include <type_traits>
#include <vector>

extern "C" void TestArrayGroup_Initialize(cGH *cctkGH) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS_TestArrayGroup_Initialize;

  int imax = 5;
  int jmax = 6;
  int nmax = 4;

  // Initialize test array
  for(int n=0; n<nmax; n++) {
    for(int j=0; j<jmax; j++) {
      for(int i=0; i<imax; i++) {
        int index = i + j*imax + n*imax*jmax;
        test1[index] = 1+i*j*n;
        test2[index] = 1+7*i*j*n;
        test3[index] = 1+13*i*j*n;
      }
    }
  }

  // Initialize test grid function
  int gi = CCTK_GroupIndex("TestArrayGroup::test_gf");
  cGroup group;
  int ierr = CCTK_GroupData(gi, &group);
  assert(!ierr);
  std::vector<int> lsh(3);
  ierr = CCTK_GrouplshGI(cctkGH, group.dim, lsh.data(), gi);
  assert(!ierr);

 for(int k=0; k < lsh[2]; k++) {
   for(int j=0; j < lsh[1]; j++) {
     for(int i=0; i < lsh[0]; i++) {
       int index = CCTK_GFINDEX3D(cctkGH, i, j, k);
       test_gf[index] = i*j*k;
     }
   }
 }

  *test_scalar = 1;

 return;
}

extern "C" void TestArrayGroup_Compare(cGH *cctkGH) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS_TestArrayGroup_Compare;

  int imax = 5;
  int jmax = 6;
  int nmax = 4;

  int error_count[3] = {0,0,0};
  std::string vname;

  for(int n=0; n<nmax; n++) {
    vname = "TestArrayGroup::test1[" + std::to_string(n) + "]";
    CCTK_REAL *var1 = (CCTK_REAL *)CCTK_VarDataPtr(cctkGH,0,vname.c_str());

    vname = "TestArrayGroup::test2[" + std::to_string(n) + "]";
    CCTK_REAL *var2 = (CCTK_REAL *)CCTK_VarDataPtr(cctkGH,0,vname.c_str());

    vname = "TestArrayGroup::test3[" + std::to_string(n) + "]";
    CCTK_REAL *var3 = (CCTK_REAL *)CCTK_VarDataPtr(cctkGH,0,vname.c_str());
    for(int j=0; j<jmax; j++) {
      for(int i=0; i<imax; i++) {
        int index = i + j*imax;
        if (var1[index] != 1+i*j*n) error_count[0] += 1;
        if (var2[index] != 1+7*i*j*n) error_count[1] += 1;
        if (var3[index] != 1+13*i*j*n) error_count[2] += 1;
      }
    }
  }

  if (error_count[0] > 0) {
    int size = nmax*jmax*imax;
    std::string msg;
    msg = "TestArrayGroup: grid array test1 failed in " + std::to_string(error_count[0]) + 
          " of " + std::to_string(size) + " elements.";
    CCTK_ERROR(msg.c_str());
  }

  if (error_count[1] > 0) {
    int size = nmax*jmax*imax;
    std::string msg;
    msg = "TestArrayGroup: grid array test2 failed in " + std::to_string(error_count[1]) + 
          " of " + std::to_string(size) + " elements.";
    CCTK_ERROR(msg.c_str());
  }

  if (error_count[2] > 0) {
    int size = nmax*jmax*imax;
    std::string msg;
    msg = "TestArrayGroup: grid array test3 failed in " + std::to_string(error_count[2]) + 
          " of " + std::to_string(size) + " elements.";
    CCTK_ERROR(msg.c_str());
  }
 return;
}
