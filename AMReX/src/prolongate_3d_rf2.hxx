#ifndef PROLONGATE_3D_RF2_HXX
#define PROLONGATE_3D_RF2_HXX

#include "driver.hxx"

#include <AMReX_Interpolater.H>

#include <cassert>

namespace AMReX {
using namespace amrex;
using namespace std;

template <int CENTI, int CENTJ, int CENTK, bool CONSI, bool CONSJ, bool CONSK,
          int ORDERI, int ORDERJ, int ORDERK>
class prolongate_3d_rf2 final : public Interpolater {

  // Centering must be vertex (0) or cell (1)
  static_assert(CENTI == 0 || CENTI == 1, "");
  static_assert(CENTJ == 0 || CENTJ == 1, "");
  static_assert(CENTK == 0 || CENTK == 1, "");

  // Order must be nonnegative
  static_assert(ORDERI >= 0, "");
  static_assert(ORDERJ >= 0, "");
  static_assert(ORDERK >= 0, "");

  static constexpr array<int, dim> indextype() { return {CENTI, CENTJ, CENTK}; }
  static constexpr array<bool, dim> conservative() {
    return {CONSI, CONSJ, CONSK};
  }
  static constexpr array<int, dim> order() { return {ORDERI, ORDERJ, ORDERK}; }

public:
  virtual ~prolongate_3d_rf2() override;

  virtual Box CoarseBox(const Box &fine, int ratio) override;
  virtual Box CoarseBox(const Box &fine, const IntVect &ratio) override;

  virtual void interp(const FArrayBox &crse, int crse_comp, FArrayBox &fine,
                      int fine_comp, int ncomp, const Box &fine_region,
                      const IntVect &ratio, const Geometry &crse_geom,
                      const Geometry &fine_geom, Vector<BCRec> const &bcr,
                      int actual_comp, int actual_state,
                      RunOn gpu_or_cpu) override;
};

extern prolongate_3d_rf2<0, 0, 0, false, false, false, 1, 1, 1>
    prolongate_3d_rf2_c000_o1;
extern prolongate_3d_rf2<0, 0, 1, false, false, false, 1, 1, 1>
    prolongate_3d_rf2_c001_o1;
extern prolongate_3d_rf2<0, 1, 0, false, false, false, 1, 1, 1>
    prolongate_3d_rf2_c010_o1;
extern prolongate_3d_rf2<0, 1, 1, false, false, false, 1, 1, 1>
    prolongate_3d_rf2_c011_o1;
extern prolongate_3d_rf2<1, 0, 0, false, false, false, 1, 1, 1>
    prolongate_3d_rf2_c100_o1;
extern prolongate_3d_rf2<1, 0, 1, false, false, false, 1, 1, 1>
    prolongate_3d_rf2_c101_o1;
extern prolongate_3d_rf2<1, 1, 0, false, false, false, 1, 1, 1>
    prolongate_3d_rf2_c110_o1;
extern prolongate_3d_rf2<1, 1, 1, false, false, false, 1, 1, 1>
    prolongate_3d_rf2_c111_o1;

extern prolongate_3d_rf2<0, 0, 0, true, true, true, 0, 0, 0>
    prolongate_cons_3d_rf2_c000_o0;
extern prolongate_3d_rf2<0, 0, 1, true, true, true, 0, 0, 0>
    prolongate_cons_3d_rf2_c001_o0;
extern prolongate_3d_rf2<0, 1, 0, true, true, true, 0, 0, 0>
    prolongate_cons_3d_rf2_c010_o0;
extern prolongate_3d_rf2<0, 1, 1, true, true, true, 0, 0, 0>
    prolongate_cons_3d_rf2_c011_o0;
extern prolongate_3d_rf2<1, 0, 0, true, true, true, 0, 0, 0>
    prolongate_cons_3d_rf2_c100_o0;
extern prolongate_3d_rf2<1, 0, 1, true, true, true, 0, 0, 0>
    prolongate_cons_3d_rf2_c101_o0;
extern prolongate_3d_rf2<1, 1, 0, true, true, true, 0, 0, 0>
    prolongate_cons_3d_rf2_c110_o0;
extern prolongate_3d_rf2<1, 1, 1, true, true, true, 0, 0, 0>
    prolongate_cons_3d_rf2_c111_o0;

extern prolongate_3d_rf2<0, 0, 0, true, true, true, 1, 1, 1>
    prolongate_cons_3d_rf2_c000_o1;
extern prolongate_3d_rf2<0, 0, 1, true, true, true, 1, 1, 2>
    prolongate_cons_3d_rf2_c001_o1;
extern prolongate_3d_rf2<0, 1, 0, true, true, true, 1, 2, 1>
    prolongate_cons_3d_rf2_c010_o1;
extern prolongate_3d_rf2<0, 1, 1, true, true, true, 1, 2, 2>
    prolongate_cons_3d_rf2_c011_o1;
extern prolongate_3d_rf2<1, 0, 0, true, true, true, 2, 1, 1>
    prolongate_cons_3d_rf2_c100_o1;
extern prolongate_3d_rf2<1, 0, 1, true, true, true, 2, 1, 2>
    prolongate_cons_3d_rf2_c101_o1;
extern prolongate_3d_rf2<1, 1, 0, true, true, true, 2, 2, 1>
    prolongate_cons_3d_rf2_c110_o1;
extern prolongate_3d_rf2<1, 1, 1, true, true, true, 2, 2, 2>
    prolongate_cons_3d_rf2_c111_o1;

} // namespace AMReX

#endif // #ifndef PROLONGATE_3D_RF2_HXX
