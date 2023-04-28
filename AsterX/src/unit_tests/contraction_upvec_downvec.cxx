#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#include "../test.hxx"
#include "../utils.hxx"

#include <random>

void AsterXTests::test_contraction_upvec_downvec() {
  using namespace Arith;
  using namespace AsterX;

  using std::mt19937_64;
  using std::uniform_real_distribution;

  // Construct symmetric matrix of random real elements
  mt19937_64 engine{random_seed};
  uniform_real_distribution<CCTK_REAL> distrib{-1.0, 1.0};

  const vec<CCTK_REAL, 3> v_up{distrib(engine), distrib(engine),
                               distrib(engine)};

  const vec<CCTK_REAL, 3> v_down{distrib(engine), distrib(engine),
                                 distrib(engine)};

  CCTK_VINFO("Testing up_vec - down_vec contraction");
  {
    bool success{true};
    const auto actual_contraction{calc_contraction(v_up, v_down)};
    const auto expected_contraction{v_up(0) * v_down(0) + v_up(1) * v_down(1) +
                                    v_up(2) * v_down(2)};

    if (!isapprox(actual_contraction, expected_contraction)) {
      CCTK_VINFO("  \033[31;1mFAILED\033[0m. Reason: Contraction expected to "
                 "be %.16f, but instead got %.16f",
                 expected_contraction, actual_contraction);
      success &= false;
    }

    if (success) {
      CCTK_VINFO("  \033[32;1mPASSED\033[0m.");
    }
  }
}