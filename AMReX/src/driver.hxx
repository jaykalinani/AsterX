#ifndef DRIVER_HXX
#define DRIVER_HXX

#include <cctk.h>
#undef copysign
#undef fpclassify
#undef isfinite
#undef isinf
#undef isnan
#undef isnormal
#undef signbit

#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFab.H>

#include <array>
#include <iostream>
#include <memory>
#include <type_traits>
#include <vector>

namespace AMReX {
using namespace amrex;
using namespace std;

constexpr int dim = 3;

static_assert(AMREX_SPACEDIM == dim,
              "AMReX's AMREX_SPACEDIM must be the same as Cactus's cctk_dim");

static_assert(is_same<Real, CCTK_REAL>::value,
              "AMReX's Real type must be the same as Cactus's CCTK_REAL");

class CactusAmrCore final : public AmrCore {
public:
  CactusAmrCore();
  CactusAmrCore(const RealBox *rb, int max_level_in,
                const Vector<int> &n_cell_in, int coord = -1,
                Vector<IntVect> ref_ratios = Vector<IntVect>(),
                const int *is_per = nullptr);
  CactusAmrCore(const RealBox &rb, int max_level_in,
                const Vector<int> &n_cell_in, int coord,
                Vector<IntVect> const &ref_ratios,
                Array<int, AMREX_SPACEDIM> const &is_per);
  CactusAmrCore(const AmrCore &rhs) = delete;
  CactusAmrCore &operator=(const AmrCore &rhs) = delete;

  virtual ~CactusAmrCore() override;

  virtual void ErrorEst(int level, TagBoxArray &tags, Real time,
                        int ngrow) override;
  virtual void MakeNewLevelFromScratch(int level, Real time, const BoxArray &ba,
                                       const DistributionMapping &dm) override;
  virtual void MakeNewLevelFromCoarse(int level, Real time, const BoxArray &ba,
                                      const DistributionMapping &dm) override;
  virtual void RemakeLevel(int level, Real time, const BoxArray &ba,
                           const DistributionMapping &dm) override;
  virtual void ClearLevel(int level) override;
};

struct valid_t {
  bool valid_int, valid_bnd;
  valid_t() : valid_int(false), valid_bnd(false) {}

  friend bool operator<(const valid_t &x, const valid_t &y) {
    if (x.valid_int < y.valid_int)
      return true;
    return x.valid_bnd < y.valid_bnd;
  }
  friend ostream &operator<<(ostream &os, const valid_t v) {
    auto str = [](bool v) { return v ? "VAL" : "INV"; };
    return os << "[int:" << str(v.valid_int) << ",bnd:" << str(v.valid_bnd)
              << "]";
  }
  operator string() const {
    ostringstream buf;
    buf << *this;
    return buf.str();
  }
};

// Cactus grid hierarchy extension
struct GHExt {

  // AMReX grid structure
  // TODO: Remove unique_ptr once AmrCore has move constructors
  unique_ptr<CactusAmrCore> amrcore;

  struct LevelData {
    int level;
    // Empty MultiFab holding a cell-centred BoxArray for iterating
    // over grid functions.
    // TODO: Can we store the BoxArray directly?
    // The data in mfab0 is not valid. It's a dummy
    // variable to get the distribution mapping, the
    // grid size, ghost zones, etc. Only needed for
    // cctkGH.
    unique_ptr<MultiFab> mfab0;

    struct GroupData {
      int groupindex;
      int firstvarindex;
      int numvars;
      array<int, dim> indextype;

      // each MultiFab has numvars components
      vector<unique_ptr<MultiFab> > mfab; // [time level]
      vector<vector<valid_t> > valid;     // [time level][var index]

      // flux register between this and the next coarser level
      unique_ptr<FluxRegister> freg;
      // associated flux group indices
      array<int, dim> fluxes; // [dir]
    };
    vector<GroupData> groupdata;
  };
  vector<LevelData> leveldata;
};

extern unique_ptr<GHExt> ghext;

Interpolater *get_interpolator(const array<int, dim> indextype);

} // namespace AMReX

#endif // #ifndef DRIVER_HXX
