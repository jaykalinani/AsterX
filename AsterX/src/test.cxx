#include <cctk.h>
#include <cctk_Arguments.h>

#include "test.hxx"

extern "C" void AsterX_Test(CCTK_ARGUMENTS) {
  using namespace AsterXTests;

  DECLARE_CCTK_ARGUMENTS;

  test_contraction_smat_upvec();
  test_contraction_upvec_downvec();

  test_corss_product();

  test_minmod();
}
