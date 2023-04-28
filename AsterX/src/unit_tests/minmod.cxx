#include <cctk.h>
#include <cctk_Arguments.h>

#include "../test.hxx"
#include "../reconstruct.hxx"

#include <random>

void AsterXTests::test_minmod() {
  using namespace AsterX;
  using std::max;
  using std::min;

  using std::mt19937_64;
  using std::uniform_real_distribution;

  // Construct symmetric matrix of random real elements
  mt19937_64 engine{random_seed};
  uniform_real_distribution<CCTK_REAL> positive_distrib{0.0, 1.0};
  uniform_real_distribution<CCTK_REAL> negative_distrib{-1.0, 0.0};

  CCTK_VINFO("Testing minmod");
  {
    bool success{true};

    const auto p1{positive_distrib(engine)};
    const auto p2{positive_distrib(engine)};

    const auto n1{negative_distrib(engine)};
    const auto n2{negative_distrib(engine)};

    const auto positive_minmod{minmod(p1, p2)};
    const auto negative_minmod{minmod(n1, n2)};

    const auto mixed_minmod_1{minmod(n1, p1)};
    const auto mixed_minmod_2{minmod(n2, p2)};

    const auto min_p1p2{min(p1, p2)};
    const auto max_n1n2{max(n1, n2)};

    if (!isapprox(positive_minmod, min_p1p2)) {
      CCTK_VINFO("  \033[31;1mFAILED\033[0m. Reason: minmod with positive "
                 "numbers expected %.f16 and got %.f16",
                 min_p1p2, positive_minmod);
      success &= false;
    }

    if (!isapprox(negative_minmod, max_n1n2)) {
      CCTK_VINFO("  \033[31;1mFAILED\033[0m. Reason: minmod with negative "
                 "numbers expected %.f16 and got %.f16",
                 min_p1p2, negative_minmod);
      success &= false;
    }

    if (!isapprox(mixed_minmod_1, CCTK_REAL(0))) {
      CCTK_VINFO("  \033[31;1mFAILED\033[0m. Reason: minmod with mixed sign "
                 "numbers numbers expected 0.0 and got %.f16",
                 mixed_minmod_1);
      success &= false;
    }

    if (!isapprox(mixed_minmod_2, CCTK_REAL(0))) {
      CCTK_VINFO("  \033[31;1mFAILED\033[0m. Reason: minmod with mixed sign "
                 "numbers numbers expected 0.0 and got %.f16",
                 mixed_minmod_2);
      success &= false;
    }

    if (success) {
      CCTK_VINFO("  \033[32;1mPASSED\033[0m.");
    }
  }
}