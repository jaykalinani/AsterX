#include <cctk.h>
#include <cctk_Arguments.h>

#include "../test.hxx"
#include "../reconstruct.hxx"

#include <random>

void AsterXTests::test_minmod(std::mt19937_64 &engine, int repetitions) {
  using namespace ReconX;
  using std::abs;
  using std::max;
  using std::min;
  using std::uniform_real_distribution;

  uniform_real_distribution<CCTK_REAL> positive_distrib{0.0, 1.0};
  uniform_real_distribution<CCTK_REAL> negative_distrib{-1.0, 0.0};
  uniform_real_distribution<CCTK_REAL> real_distrib{-1.0, 1.0};

  for (int i = 0; i < repetitions; i++) {

    CCTK_VINFO("Testing minmod repetition %i", i);
    {
      bool success{true};

      const auto p1{positive_distrib(engine)};
      const auto p2{positive_distrib(engine)};

      const auto n1{negative_distrib(engine)};
      const auto n2{negative_distrib(engine)};

      const auto r1{real_distrib(engine)};
      const auto r2{real_distrib(engine)};

      const auto positive_minmod{minmod(p1, p2)};
      const auto negative_minmod{minmod(n1, n2)};

      const auto mixed_minmod_1{minmod(n1, p1)};
      const auto mixed_minmod_2{minmod(n2, p2)};

      const auto real_minmod(minmod(r1, r2));

      const auto min_p1p2{min(p1, p2)};
      const auto max_n1n2{max(n1, n2)};

      const auto max_abs{max(abs(r1), abs(r2))};

      if (!isapprox(positive_minmod, min_p1p2)) {
        CCTK_VINFO("  FAILED. Reason: minmod with positive "
                   "numbers expected %.f16 and got %.f16",
                   min_p1p2, positive_minmod);
        success &= false;
      }

      if (!isapprox(negative_minmod, max_n1n2)) {
        CCTK_VINFO("  FAILED. Reason: minmod with negative "
                   "numbers expected %.f16 and got %.f16",
                   min_p1p2, negative_minmod);
        success &= false;
      }

      if (!isapprox(mixed_minmod_1, CCTK_REAL(0))) {
        CCTK_VINFO("  FAILED. Reason: minmod with mixed sign "
                   "numbers numbers expected 0.0 and got %.f16",
                   mixed_minmod_1);
        success &= false;
      }

      if (!isapprox(mixed_minmod_2, CCTK_REAL(0))) {
        CCTK_VINFO("  FAILED. Reason: minmod with mixed sign "
                   "numbers numbers expected 0.0 and got %.f16",
                   mixed_minmod_2);
        success &= false;
      }

      if (!(real_minmod < max_abs || isapprox(real_minmod, max_abs))) {
        CCTK_VINFO("  FAILED. Reason: minmod with reals %.f16 and %.f16 "
                   "returned %.f, which is not bounded by %.f16",
                   r1, r2, real_minmod, max_abs);
        success &= false;
      }

      if (success) {
        CCTK_VINFO("  PASSED.");
      }
    }
  }
}
