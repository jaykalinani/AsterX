#ifndef DRIVER_HXX
#define DRIVER_HXX

#include <AMReX.H>
#include <AMReX_AmrMesh.H>
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

class CactusAmrMesh : public AmrMesh {
public:
  CactusAmrMesh();
  CactusAmrMesh(const RealBox *rb, int max_level_in,
                const Vector<int> &n_cell_in, int coord = -1,
                Vector<IntVect> ref_ratios = Vector<IntVect>(),
                const int *is_per = nullptr);
  CactusAmrMesh(const RealBox &rb, int max_level_in,
                const Vector<int> &n_cell_in, int coord,
                Vector<IntVect> const &ref_ratios,
                Array<int, AMREX_SPACEDIM> const &is_per);
  CactusAmrMesh(const AmrMesh &rhs) = delete;
  CactusAmrMesh &operator=(const AmrMesh &rhs) = delete;

  virtual ~CactusAmrMesh();

  // virtual void MakeNewLevelFromScratch(int lev, Real time, const BoxArray
  // &ba,
  //                                      const DistributionMapping &dm)
  //                                      override;
  virtual void ErrorEst(int lev, TagBoxArray &tags, Real time,
                        int ngrow) override;
  // virtual void ManualTagsPlacement(int lev, TagBoxArray &tags,
  //                                  const Vector<IntVect> &bf_lev) override;
  // virtual BoxArray GetAreaNotToTag(int lev) override;
};

struct GHExt {

  // AMReX grid structure
  unique_ptr<CactusAmrMesh> amrmesh;

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

} // namespace AMReX

#endif // #ifndef DRIVER_HXX
