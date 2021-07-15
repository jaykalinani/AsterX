#include "pdesolvers.hxx"

// TODO: Don't include files from other thorn; create a proper interface
#include "../../CarpetX/src/driver.hxx"
#include "../../CarpetX/src/reduction.hxx"

#include <loop.hxx>
#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <mpi.h>

#include <petscdm.h>
#include <petscsnes.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <optional>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

namespace PDESolvers {

const int tl = 0;

// TODO: Generalize this
const Arith::vect<int, 3> indextype{0, 0, 0}; // vertex centred

// Tuning parameters for matrix preallocation
// const int dnz = 7;
// const int onz = 4;
const int dnz = 13;
const int onz = 8;

////////////////////////////////////////////////////////////////////////////////

std::optional<jacobian_t> jacobian;

////////////////////////////////////////////////////////////////////////////////

// Erik Schnetter the blocks for vertex centred grids overlap by one
// point at the block boundaries (see
// https://amrex-codes.github.io/amrex/docs_html/Basics.html#id4).
// does FillPatchSingleLevel ensure these points are consistent? or is
// the domain double-valued there, and only the ghost points (not
// shown in this figure) are set?

// weiqun If the original source data are consistent, the results of
// FillPatchSingleLevel are also consistent.  Otherwise, it's undefined which
// original values will used for which boxes.  It depends on how boxes are
// distributed.

// If you want to the data consistent, you can use MultiFab's OverrideSync
// functions.
//     void OverrideSync (const Periodicity& period =
//         Periodicity::NonPeriodic());
//     void OverrideSync (const iMultiFab& msk,
//         const Periodicity& period = Periodicity::NonPeriodic());

// The firs version will use data from the grid with the lowest grid number to
// override others.

// In the second version, you can pass your own mask.

// Erik Schnetter thanks!

// weiqun Even if you want to use the first version's mask, you might still want
// to call
//     //! Owner is the grid with the lowest grid number containing the data.
//     std::unique_ptr<iMultiFab> OwnerMask (const Periodicity& period =
//         Periodicity::NonPeriodic()) const;
// and then use the second function.

// If you need to call OverrideSync repeatedly, you could save the mask.

// Erik Schnetter thanks!

////////////////////////////////////////////////////////////////////////////////

template <typename F>
void loop_over_points(const int level, const int block,
                      const Arith::vect<int, 3> &lo,
                      const Arith::vect<int, 3> &hi, const F &point_kernel) {
  for (int k = lo[2]; k < hi[2]; ++k) {
    for (int j = lo[1]; j < hi[1]; ++j) {
      for (int i = lo[0]; i < hi[0]; ++i) {
        point_kernel(i, j, k);
      }
    }
  }
}

template <typename F, typename = std::invoke_result_t<
                          F, int, int, const Arith::vect<int, 3> &,
                          const Arith::vect<int, 3> &> >
void loop_over_blocks(const F &block_kernel) {
  // Decode Cactus variables
  static const int vn = CCTK_VarIndex("PDESolvers::idx");
  assert(vn >= 0);
  static const int gi = CCTK_GroupIndexFromVarI(vn);
  assert(gi >= 0);
  static const int v0 = CCTK_FirstVarIndexI(gi);
  assert(v0 >= 0);
  static const int vi = vn - v0;
  assert(vi >= 0);

  std::vector<std::function<void()> > tasks;

  // Set Cactus index vector
  for (const auto &leveldata : CarpetX::ghext->leveldata) {
    const int level = leveldata.level;
    amrex::MultiFab &mfab = *leveldata.groupdata.at(gi)->mfab.at(tl);
    for (int d = 0; d < 3; ++d)
      assert(mfab.ixType()[d] == (indextype[d]
                                      ? amrex::IndexType::CellIndex::CELL
                                      : amrex::IndexType::CellIndex::NODE));
    const int nblocks = mfab.local_size();
    for (int block = 0; block < nblocks; ++block) {
      tasks.emplace_back([&block_kernel, &mfab, level, block]() {
        const int global_block = mfab.IndexArray().at(block);
        const amrex::Box &box = mfab.box(global_block);
        const Arith::vect<int, 3> lo{box.smallEnd()[0], box.smallEnd()[1],
                                     box.smallEnd()[2]};
        const Arith::vect<int, 3> hi{box.bigEnd()[0] + 1, box.bigEnd()[1] + 1,
                                     box.bigEnd()[2] + 1};
        block_kernel(level, block, lo, hi);
      });
    } // for block
  }   // for level

#pragma omp parallel for schedule(dynamic, 1)
  for (int n = 0; n < int(tasks.size()); ++n)
    tasks.at(n)();
}

////////////////////////////////////////////////////////////////////////////////

void count_points(const std::vector<int> &varinds, int &npoints_local,
                  std::vector<std::vector<int> > &block_offsets) {
  int myproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &myproc);

  // Decode Cactus variables
  static const int vn = CCTK_VarIndex("PDESolvers::idx");
  assert(vn >= 0);
  static const int gi = CCTK_GroupIndexFromVarI(vn);
  assert(gi >= 0);
  static const int v0 = CCTK_FirstVarIndexI(gi);
  assert(v0 >= 0);
  static const int vi = vn - v0;
  assert(vi >= 0);

  const int nvars = varinds.size();

  npoints_local = 0;
  block_offsets.resize(CarpetX::ghext->leveldata.size());
  for (const auto &leveldata : CarpetX::ghext->leveldata) {
    const int level = leveldata.level;
    amrex::MultiFab &mfab = *leveldata.groupdata.at(gi)->mfab.at(tl);
    for (int d = 0; d < 3; ++d)
      assert(mfab.ixType()[d] == (indextype[d]
                                      ? amrex::IndexType::CellIndex::CELL
                                      : amrex::IndexType::CellIndex::NODE));
    const int nblocks = mfab.local_size();
    block_offsets.at(level).resize(nblocks);
    for (int block = 0; block < nblocks; ++block) {
      block_offsets.at(level).at(block) = npoints_local;
      const int global_block = mfab.IndexArray().at(block);
      assert(mfab.DistributionMap().ProcessorMap().at(global_block) == myproc);
      const amrex::Box &box = mfab.box(global_block);
      const Arith::vect<int, 3> lo{box.smallEnd()[0], box.smallEnd()[1],
                                   box.smallEnd()[2]};
      const Arith::vect<int, 3> hi{box.bigEnd()[0] + 1, box.bigEnd()[1] + 1,
                                   box.bigEnd()[2] + 1};
      npoints_local += nvars * prod(hi - lo);
    }
  }
}

void set_global_index(const std::vector<int> &varinds, const int npoints_local,
                      const int global_offset,
                      const std::vector<std::vector<int> > &block_offsets) {
  int myproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &myproc);

  // Decode Cactus variables
  const int nvars = varinds.size();
  const int vn = CCTK_VarIndex("PDESolvers::idx");
  assert(vn >= 0);
  const int gi = CCTK_GroupIndexFromVarI(vn);
  assert(gi >= 0);
  const int v0 = CCTK_FirstVarIndexI(gi);
  assert(v0 >= 0);
  const int vi = vn - v0;
  assert(vi >= 0);

  // Set Cactus index vector
  loop_over_blocks([&](const int level, const int block,
                       const Arith::vect<int, 3> &lo,
                       const Arith::vect<int, 3> &hi) {
    const int block_offset = global_offset + block_offsets.at(level).at(block);
    amrex::MultiFab &mfab =
        *CarpetX::ghext->leveldata.at(level).groupdata.at(gi)->mfab.at(tl);
    amrex::Array4<CCTK_REAL> cactus_arr = mfab.atLocalIdx(block).array(vi);
    const Arith::vect<int, 3> sh = hi - lo;
    const Arith::vect<int, 3> di{nvars, nvars * sh[0], nvars * sh[0] * sh[1]};
    const int offset =
        block_offset - di[0] * lo[0] - di[1] * lo[1] - di[2] * lo[2];
    loop_over_points(
        level, block, lo, hi, [&](const int i, const int j, const int k) {
          cactus_arr(i, j, k) = offset + di[0] * i + di[1] * j + di[2] * k;
        });
  });
}

void copy_Cactus_to_PETSc(Vec vec, const std::vector<int> &varinds,
                          const std::vector<std::vector<int> > &block_offsets) {
  PetscErrorCode ierr;

  int myproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &myproc);

  // Decode Cactus variables
  const int nvars = varinds.size();
  std::vector<int> gis, vis;
  gis.reserve(nvars);
  vis.reserve(nvars);
  for (int vn : varinds) {
    assert(vn >= 0);
    const int gi = CCTK_GroupIndexFromVarI(vn);
    assert(gi >= 0);
    gis.push_back(gi);
    const int v0 = CCTK_FirstVarIndexI(gi);
    assert(v0 >= 0);
    const int vi = vn - v0;
    assert(vi >= 0);
    vis.push_back(vi);
  }

  // Get vector from PETSc
  CCTK_REAL *vec_ptr;
  ierr = VecGetArray(vec, &vec_ptr);
  assert(!ierr);
  CCTK_REAL *restrict const petsc_ptr = vec_ptr;

  // Copy vector to Cactus
  loop_over_blocks([&](const int level, const int block,
                       const Arith::vect<int, 3> &lo,
                       const Arith::vect<int, 3> &hi) {
    const int block_offset = block_offsets.at(level).at(block);
    std::vector<amrex::Array4<const CCTK_REAL> > cactus_arrs;
    cactus_arrs.reserve(nvars);
    for (int n = 0; n < nvars; ++n) {
      const amrex::MultiFab &mfab =
          *CarpetX::ghext->leveldata.at(level).groupdata.at(gis.at(n))->mfab.at(
              tl);
      cactus_arrs.push_back(mfab.atLocalIdx(block).array(vis.at(n)));
    }
    const Arith::vect<int, 3> sh = hi - lo;
    const Arith::vect<int, 3> di{nvars, nvars * sh[0], nvars * sh[0] * sh[1]};
    const int offset =
        block_offset - di[0] * lo[0] - di[1] * lo[1] - di[2] * lo[2];
    loop_over_points(
        level, block, lo, hi, [&](const int i, const int j, const int k) {
          for (int n = 0; n < nvars; ++n)
            petsc_ptr[offset + n + di[0] * i + di[1] * j + di[2] * k] =
                cactus_arrs.at(n)(i, j, k);
        });
  });

  ierr = VecRestoreArray(vec, &vec_ptr);
  assert(!ierr);
}

void copy_PETSc_to_Cactus(Vec vec, const std::vector<int> &varinds,
                          const std::vector<std::vector<int> > &block_offsets) {
  PetscErrorCode ierr;

  int myproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &myproc);

  // Decode Cactus variables
  const int nvars = varinds.size();
  std::vector<int> gis, vis;
  gis.reserve(nvars);
  vis.reserve(nvars);
  for (int vn : varinds) {
    assert(vn >= 0);
    const int gi = CCTK_GroupIndexFromVarI(vn);
    assert(gi >= 0);
    gis.push_back(gi);
    const int v0 = CCTK_FirstVarIndexI(gi);
    assert(v0 >= 0);
    const int vi = vn - v0;
    assert(vi >= 0);
    vis.push_back(vi);
  }

  // Get vector from PETSc
  const CCTK_REAL *vec_ptr;
  ierr = VecGetArrayRead(vec, &vec_ptr);
  assert(!ierr);
  const CCTK_REAL *restrict const petsc_ptr = vec_ptr;

  // Copy vector to Cactus
  loop_over_blocks([&](const int level, const int block,
                       const Arith::vect<int, 3> &lo,
                       const Arith::vect<int, 3> &hi) {
    const int block_offset = block_offsets.at(level).at(block);
    std::vector<amrex::Array4<CCTK_REAL> > cactus_arrs;
    cactus_arrs.reserve(nvars);
    for (int n = 0; n < nvars; ++n) {
      amrex::MultiFab &mfab =
          *CarpetX::ghext->leveldata.at(level).groupdata.at(gis.at(n))->mfab.at(
              tl);
      cactus_arrs.push_back(mfab.atLocalIdx(block).array(vis.at(n)));
    }
    const Arith::vect<int, 3> sh = hi - lo;
    const Arith::vect<int, 3> di{nvars, nvars * sh[0], nvars * sh[0] * sh[1]};
    const int offset =
        block_offset - di[0] * lo[0] - di[1] * lo[1] - di[2] * lo[2];
    loop_over_points(
        level, block, lo, hi, [&](const int i, const int j, const int k) {
          for (int n = 0; n < nvars; ++n)
            cactus_arrs.at(n)(i, j, k) =
                petsc_ptr[offset + n + di[0] * i + di[1] * j + di[2] * k];
        });
  });

  ierr = VecRestoreArrayRead(vec, &vec_ptr);
  assert(!ierr);
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void PDESolvers_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_PDESolvers_Setup;
  DECLARE_CCTK_PARAMETERS;

  /*static*/ std::vector<std::string> args;
  args.push_back("cactus");
  std::istringstream buf(petsc_options);
  std::copy(std::istream_iterator<std::string>(buf),
            std::istream_iterator<std::string>(), std::back_inserter(args));

  /*static*/ std::vector<char *> argptrs;
  for (auto &arg : args)
    argptrs.push_back(arg.data());
  argptrs.push_back(nullptr);

  /*static*/ int argc = args.size();
  /*static*/ char **argv = argptrs.data();
  PetscErrorCode ierr = PetscInitialize(&argc, &argv, NULL, NULL);
  assert(!ierr);
}

extern "C" void PDESolvers_IdxInit(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_PDESolvers_IdxInit;

  const int dim = Loop::dim;

  const std::array<int, dim> nghostzones = {
      cctk_nghostzones[0], cctk_nghostzones[1], cctk_nghostzones[2]};
  Arith::vect<int, dim> imin, imax;
  Loop::GridDescBase(cctkGH).box_int<0, 0, 0>(nghostzones, imin, imax);
  const Loop::GF3D2layout layout1(cctkGH, indextype);

  const Loop::GF3D2<CCTK_REAL> gf_idx(layout1, idx);

  const Loop::GridDescBase grid(cctkGH);
  grid.loop_int<0, 0, 0>(grid.nghostzones,
                         [=](const Loop::PointDesc &p) ARITH_INLINE {
                           const Loop::GF3D2index index1(layout1, p.I);
                           using std::lrint;
                           gf_idx(index1) = -1;
                         });
}

extern "C" void PDESolvers_IdxBoundary(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_PDESolvers_IdxBoundary;

  const int dim = Loop::dim;

  const std::array<int, dim> nghostzones = {
      cctk_nghostzones[0], cctk_nghostzones[1], cctk_nghostzones[2]};
  Arith::vect<int, dim> imin, imax;
  Loop::GridDescBase(cctkGH).box_int<0, 0, 0>(nghostzones, imin, imax);
  const Loop::GF3D2layout layout1(cctkGH, indextype);

  const Loop::GF3D2<CCTK_REAL> gf_idx(layout1, idx);

  const Loop::GridDescBase grid(cctkGH);
  grid.loop_bnd<0, 0, 0>(grid.nghostzones,
                         [=](const Loop::PointDesc &p) ARITH_INLINE {
                           const Loop::GF3D2index index1(layout1, p.I);
                           using std::lrint;
                           gf_idx(index1) = -1;
                         });
}

namespace {
PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void *fun) {
  return (
      *static_cast<std::function<PetscErrorCode(SNES snes, Vec x, Vec f)> *>(
          fun))(snes, x, f);
}

PetscErrorCode FormJacobian(SNES snes, Vec x, Mat J, Mat B, void *fun) {
  return (*static_cast<
          std::function<PetscErrorCode(SNES snes, Vec x, Mat J, Mat B)> *>(
      fun))(snes, x, J, B);
}
} // namespace

extern "C" void PDESolvers_Solve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_PDESolvers_Solve;
  DECLARE_CCTK_PARAMETERS;

  PetscErrorCode ierr;

  // Grid layout

  const std::vector<int> solinds{CCTK_VarIndex("Poisson2::sol")};
  const std::vector<int> resinds{CCTK_VarIndex("Poisson2::res")};

  int npoints_local;
  std::vector<std::vector<int> > block_offsets;
  count_points(solinds, npoints_local, block_offsets);
  {
    int npoints_local1;
    std::vector<std::vector<int> > block_offsets1;
    count_points(resinds, npoints_local1, block_offsets1);
    assert(npoints_local1 == npoints_local);
    assert(block_offsets1 == block_offsets);
  }
  int npoints_offset = 0;
  MPI_Exscan(&npoints_local, &npoints_offset, 1, MPI_INT, MPI_SUM,
             MPI_COMM_WORLD);
  int npoints_global;
  MPI_Allreduce(&npoints_local, &npoints_global, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  assert(npoints_offset >= 0);
  assert(npoints_offset + npoints_local <= npoints_global);
  // We don't want this call
  CallScheduleGroup(cctkGH, "PDESolvers_IdxInitGroup");
  set_global_index(solinds, npoints_local, npoints_offset, block_offsets);
  CallScheduleGroup(cctkGH, "PDESolvers_IdxBoundaryGroup");

  // Nonlinear solver

  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD, &snes);
  assert(!ierr);
  ierr = SNESSetFromOptions(snes);
  assert(!ierr);
  ierr = SNESSetAlwaysComputesFinalResidual(snes, PETSC_TRUE);
  assert(!ierr);

  // State vector and evaluation function

  Vec r;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, npoints_local, npoints_global, &r);
  assert(!ierr);
  ierr = VecSetFromOptions(r);
  assert(!ierr);

  std::function<PetscErrorCode(SNES snes, Vec x, Vec f)> evalf =
      [&](SNES snes, Vec x, Vec f) {
        copy_PETSc_to_Cactus(x, solinds, block_offsets);
        CallScheduleGroup(cctkGH, "PDESolvers_Residual");
        copy_Cactus_to_PETSc(f, resinds, block_offsets);
        return 0;
      };
  ierr = SNESSetFunction(snes, r, FormFunction, &evalf);
  assert(!ierr);

  // Matrix and Jacobian evaluation function

  Mat J;
  ierr = MatCreateAIJ(PETSC_COMM_WORLD, npoints_local, npoints_local,
                      npoints_global, npoints_global, dnz, NULL, onz, NULL, &J);
  assert(!ierr);
  ierr = MatSetFromOptions(J);
  assert(!ierr);
  jacobian = std::make_optional<jacobian_t>(J);

  std::function<PetscErrorCode(SNES snes, Vec x, Mat J, Mat B)> evalJ =
      [&](SNES snes, Vec x, Mat J, Mat B) {
        ierr = MatZeroEntries(J);
        assert(!ierr);
        copy_PETSc_to_Cactus(x, solinds, block_offsets);
        CallScheduleGroup(cctkGH, "PDESolvers_Jacobian");
        ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
        assert(!ierr);
        ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
        assert(!ierr);
        if (false) {
          MatInfo info;
          ierr = MatGetInfo(J, MAT_GLOBAL_SUM, &info);
          assert(!ierr);
          CCTK_VINFO("Jacobian info: nz_allocated=%g nz_used=%g nz_unneeded=%g "
                     "memory=%g",
                     double(info.nz_allocated), double(info.nz_used),
                     double(info.nz_unneeded), double(info.memory));
        }
        return 0;
      };
#if 1
  ierr = SNESSetJacobian(snes, J, J, FormJacobian, &evalJ);
  assert(!ierr);
#else
  ierr = SNESSetJacobian(snes, J, J, NULL, NULL);
  assert(!ierr);
#endif

  {
    double atol, rtol, stol;
    int maxit, maxf;
    ierr = SNESGetTolerances(snes, &atol, &rtol, &stol, &maxit, &maxf);
    assert(!ierr);
    CCTK_VINFO("Settings: atol=%g, rtol=%g, stol=%g, maxit=%d, maxf=%d", atol,
               rtol, stol, maxit, maxf);
  }

  // Initialize solution

  Vec x;
  ierr = VecDuplicate(r, &x);
  assert(!ierr);
  copy_Cactus_to_PETSc(x, solinds, block_offsets);

#if 0
  ierr = FormJacobian(snes, x, J, J, &evalJ);
  assert(!ierr);
#endif

  // Solve

  CCTK_INFO("Solve");
  ierr = SNESSolve(snes, NULL, x);
  assert(!ierr);
  {
    int iters;
    ierr = SNESGetIterationNumber(snes, &iters);
    assert(!ierr);
    int liters;
    ierr = SNESGetLinearSolveIterations(snes, &liters);
    assert(!ierr);
    CCTK_VINFO("Iterations: %d snes, %d linear", iters, liters);
  }

  // Extract solution

  copy_PETSc_to_Cactus(x, solinds, block_offsets);
  CallScheduleGroup(cctkGH, "PDESolvers_Residual");
  {
    const int nvars = resinds.size();
    std::vector<int> gis, vis;
    gis.reserve(nvars);
    vis.reserve(nvars);
    for (int vn : resinds) {
      assert(vn >= 0);
      const int gi = CCTK_GroupIndexFromVarI(vn);
      assert(gi >= 0);
      gis.push_back(gi);
      const int v0 = CCTK_FirstVarIndexI(gi);
      assert(v0 >= 0);
      const int vi = vn - v0;
      assert(vi >= 0);
      vis.push_back(vi);
    }
    const CarpetX::reduction<CCTK_REAL, 3> res_norm =
        CarpetX::reduce(gis.at(0), vis.at(0), tl);
    std::ostringstream buf;
    buf << res_norm;
    CCTK_VINFO("Residual norm: %s", buf.str().c_str());
  }

  // Clean up

  jacobian.reset();
  VecDestroy(&x);
  VecDestroy(&r);
  MatDestroy(&J);
  SNESDestroy(&snes);
  CCTK_INFO("Done.");
}

extern "C" void PDESolvers_Shutdown(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_PDESolvers_Shutdown;
  PetscFinalize();
}

} // namespace PDESolvers
