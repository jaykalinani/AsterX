#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include "utils.hxx"

namespace AsterX {
using namespace std;
using namespace Loop;


extern "C" void AsterX_SourceTerms(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AsterX_SourceTerms;
  DECLARE_CCTK_PARAMETERS;

}

} //namespace AsterX


