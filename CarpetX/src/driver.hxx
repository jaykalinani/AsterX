#ifndef DRIVER_HXX
#define DRIVER_HXX

#include "loop.hxx"
#include "valid.hxx"

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
#include <AMReX_FillPatchUtil.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFab.H>

#include <algorithm>
#include <array>
#include <functional>
#include <memory>
#include <ostream>
#include <tuple>
#include <type_traits>
#include <vector>

namespace CarpetX {
using namespace amrex;
using namespace std;

using Loop::dim;

// Taken from
// <https://stackoverflow.com/questions/27440953/stdunique-ptr-for-c-functions-that-need-free>
struct free_deleter {
  template <typename T> void operator()(T *p) const {
    std::free(const_cast<std::remove_const_t<T> *>(p));
  }
};
template <typename T> using unique_C_ptr = std::unique_ptr<T, free_deleter>;
static_assert(sizeof(char *) == sizeof(unique_C_ptr<char>),
              ""); // ensure no overhead

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

// Cactus grid hierarchy extension
struct GHExt {

  // AMReX grid structure
  // TODO: Remove unique_ptr once AmrCore has move constructors
  unique_ptr<CactusAmrCore> amrcore;

  struct CommonGroupData {
    int groupindex;
    int firstvarindex;
    int numvars;

    bool do_checkpoint; // whether to checkpoint
    bool do_restrict;   // whether to restrict

    vector<vector<why_valid_t> > valid; // [time level][var index]

    // TODO: add poison_invalid and check_valid functions
  };

  struct GlobalData {
    // all data that exists on all levels

    struct ScalarGroupData : public CommonGroupData {
      vector<vector<CCTK_REAL> > data; // [time level][var index]
    };
    // TODO: right now this is sized for the total number of groups
    vector<unique_ptr<ScalarGroupData> > scalargroupdata; // [group index]
  };
  GlobalData globaldata;

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

    struct GroupData : public CommonGroupData {
      array<int, dim> indextype;
      array<int, dim> nghostzones;
      vector<array<int, dim> > parities;

      // each MultiFab has numvars components
      vector<unique_ptr<MultiFab> > mfab; // [time level]

      // flux register between this and the next coarser level
      unique_ptr<FluxRegister> freg;
      // associated flux group indices
      array<int, dim> fluxes; // [dir]
    };
    // TODO: right now this is sized for the total number of groups
    vector<unique_ptr<GroupData> > groupdata; // [group index]
  };
  vector<LevelData> leveldata; // [reflevel]
};

extern unique_ptr<GHExt> ghext;

Interpolater *get_interpolator(const array<int, dim> indextype);

typedef void apply_physbcs_t(const Box &, const FArrayBox &, int, int,
                             const Geometry &, CCTK_REAL, const Vector<BCRec> &,
                             int, int);
typedef PhysBCFunct<apply_physbcs_t *> CarpetXPhysBCFunct;
tuple<CarpetXPhysBCFunct, Vector<BCRec> >
get_boundaries(const GHExt::LevelData &leveldata,
               const GHExt::LevelData::GroupData &groupdata);

} // namespace CarpetX

#endif // #ifndef DRIVER_HXX
