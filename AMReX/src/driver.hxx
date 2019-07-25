#ifndef DRIVER_HXX
#define DRIVER_HXX

#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_MultiFab.H>

#include <cctk.h>

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

class CactusAmrCore : public AmrCore {
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

  virtual ~CactusAmrCore();

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

  struct LevelData {
    int level;

    struct GroupData {
      int firstvarindex;
      int numvars;
      // each MultiFab has numvars components
      vector<unique_ptr<MultiFab> > mfab; // [time level]
    };
    vector<GroupData> groupdata;
  };
  vector<LevelData> leveldata;
};

extern unique_ptr<GHExt> ghext;

void CreateRefinedGrid(int level);

} // namespace AMReX

#endif // #ifndef DRIVER_HXX
