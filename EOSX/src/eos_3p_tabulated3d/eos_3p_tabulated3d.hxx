// =============================================================================
//  EOS: 3‑parameter tabulated (rho, T, Ye) – modern, GPU‑friendly implementation
// -----------------------------------------------------------------------------
//  • Completely self‑contained header (inline definitions) – just include & use.
//  • Works with both serial and parallel HDF5 (broadcasts data when HDF5‑MPI is
//    not available).
//  • Consistent variable ordering with Stellar‑Collapse style tables (19 fields)
//    so every dataset lands at the expected index.
//  • Robust bounds‑checking & clamping for ρ, T, ε and Y_e before every lookup.
//  • Root‑finding in log T uses Brent with automatic bracketing inside table
//    limits.
//  • SIMD‑/GPU‑friendly math (avoid branches where possible, always operate in
//    log‑space inside the table).
// -----------------------------------------------------------------------------
//  Jay Kalinani – 18 Jun 2025
// =============================================================================
#ifndef EOS_3P_TABULATED3D_HXX
#define EOS_3P_TABULATED3D_HXX

// ===== Standard / external ===================================================
#include <cassert>
#include <cmath>
#include <limits>
#include <string>
#include <array>
#include <mpi.h>
#include <hdf5.h>

// ===== Cactus / project headers =============================================
#include <cctk.h>
#include "../eos_3p.hxx"
#include "../utils/eos_brent.hxx"               // zero_brent
#include "../utils/eos_linear_interp_ND.hxx"    // linear_interp_uniform_ND_t

// --------------------------------------------------------------------------------
//  Table layout — keep this in sync with the HDF5 file produced by stellarcollapse
// --------------------------------------------------------------------------------
#define NTABLES 19

namespace EOSX {
using namespace eos_constants;   // for unit‑conversion factors

// -----------------------------------------------------------------------------
//  Helper: clamp x to the closed interval [lo, hi]
// -----------------------------------------------------------------------------
template <typename T> CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE  
static T clamp(const T x, const T lo, const T hi) {
  return (x < lo ? lo : (x > hi ? hi : x));
}

// =============================================================================
//  Class declaration / definition
// =============================================================================
class eos_3p_tabulated3d : public eos_3p {
public:
  // ----- Dataset order inside the packed table --------------------------------
  enum EV : int {
    PRESS = 0,   // log(P)
    EPS,         // log(ε)
    S,           // entropy (k_B / baryon)
    MUNU,        // μ_ν (not used by hydrodynamics but supplied)
    CS2,         // c_s^2  (physical units)
    DEDT,        // ∂ε/∂T |ρ,Ye
    DPDRHOE,     // ∂P/∂ρ |ε,Ye
    DPDERHO,     // ∂P/∂ε |ρ,Ye
    MUHAT,       // μ̂
    MU_E,
    MU_P,
    MU_N,
    XA,
    XH,
    XN,
    XP,
    ABAR,
    ZBAR,
    GAMMA,       // adiabatic index Γ (optional — not used internally)
    NUM_VARS
  };

  // ---------------------------------------------------------------------------
  //  Public interface: init() must be called once on the host (rank 0 is fine).
  // ---------------------------------------------------------------------------
  CCTK_REAL gamma;
  CCTK_REAL       *energy_shift {nullptr};
  range            rgeps;        // global ε range (after energy‑shift removal)
  linear_interp_uniform_ND_t<CCTK_REAL, 3, NTABLES> *interptable {nullptr};

  // ---------------------------------------------------------------------------
  //  init  –  read the HDF5 EOS file, convert units, build interpolation table
  // ---------------------------------------------------------------------------
  CCTK_HOST void init(const std::string &fname,
                      range &out_rgeps,
                      const range &rgrho_in,
                      const range &rgye_in) {
    // -------------------------------------------------------
    // 0.  Read grid sizes (done on rank 0 if serial I/O)
    // -------------------------------------------------------
    int rank;  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);  assert(fapl >= 0);
#ifdef H5_HAVE_PARALLEL
    CHECK_ERROR( H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL) );
    CHECK_ERROR( H5Pset_all_coll_metadata_ops(fapl, true) );
#endif

    hid_t fid  = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, fapl);
    assert(fid >= 0);

    int nrho=0, ntemp=0, nye=0;
    if (rank == 0) {
      get_hdf5_int_dset(fid, "pointsrho",  1, &nrho);
      get_hdf5_int_dset(fid, "pointstemp", 1, &ntemp);
      get_hdf5_int_dset(fid, "pointsye",   1, &nye);
    }
    MPI_Bcast(&nrho,  1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ntemp, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nye,   1, MPI_INT, 0, MPI_COMM_WORLD);

    const size_t npoints = static_cast<size_t>(nrho) * ntemp * nye;

    // -------------------------------------------------------
    // 1.  Allocate managed (device‑visible) memory
    // -------------------------------------------------------
    auto *logrho   = (CCTK_REAL*) amrex::The_Managed_Arena()->alloc(nrho  * sizeof(CCTK_REAL));
    auto *logtemp  = (CCTK_REAL*) amrex::The_Managed_Arena()->alloc(ntemp * sizeof(CCTK_REAL));
    auto *yes      = (CCTK_REAL*) amrex::The_Managed_Arena()->alloc(nye   * sizeof(CCTK_REAL));
    auto *alltabl  = (CCTK_REAL*) amrex::The_Managed_Arena()->alloc(npoints * NTABLES * sizeof(CCTK_REAL));
    energy_shift   = (CCTK_REAL*) amrex::The_Managed_Arena()->alloc(sizeof(CCTK_REAL));

    static const char *dset_names[NTABLES] = {
      "logpress", "logenergy", "entropy", "munu",   "cs2",
      "dedt",     "dpdrhoe",  "dpderho", "muhat",  "mu_e",
      "mu_p",     "mu_n",     "Xa",      "Xh",     "Xn",
      "Xp",       "Abar",     "Zbar",     "gamma" };

    // -------------------------------------------------------
    // 2.  Read datasets (MPI‑IO if available else serial‑IO + Bcast)
    // -------------------------------------------------------
#ifdef H5_HAVE_PARALLEL
    // — Parallel path: every rank reads its own data (collective)
    get_hdf5_real_dset(fid, "logrho",     nrho,  logrho);
    get_hdf5_real_dset(fid, "logtemp",    ntemp, logtemp);
    get_hdf5_real_dset(fid, "ye",         nye,   yes);
    get_hdf5_real_dset(fid, "energy_shift", 1,   energy_shift);
    for (int v=0; v<NTABLES; ++v)
      get_hdf5_real_dset(fid, dset_names[v], npoints, &alltabl[v*npoints]);
#else
    if (rank == 0) {
      get_hdf5_real_dset(fid, "logrho",     nrho,  logrho);
      get_hdf5_real_dset(fid, "logtemp",    ntemp, logtemp);
      get_hdf5_real_dset(fid, "ye",         nye,   yes);
      get_hdf5_real_dset(fid, "energy_shift", 1,   energy_shift);
      for (int v=0; v<NTABLES; ++v)
        get_hdf5_real_dset(fid, dset_names[v], npoints, &alltabl[v*npoints]);
    }
    // broadcast to rest of ranks so everyone owns the table (needed for GPU)
    MPI_Bcast(logrho,  nrho,              MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(logtemp, ntemp,             MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(yes,     nye,               MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(energy_shift, 1,            MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(alltabl, npoints*NTABLES,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

    // Everyone closes file / property list
    CHECK_ERROR( H5Fclose(fid) );
    CHECK_ERROR( H5Pclose(fapl) );

    // -------------------------------------------------------
    // 3.  Unit conversions & log‑space preparation
    // -------------------------------------------------------
    const CCTK_REAL ln10 = log(10.0);
    *energy_shift *= EPSGF;                     // erg/g → code

    // logrho, logT : base‑10 → natural log, then scale to code units
#pragma omp parallel for simd if(nrho>256)
    for (int i=0; i<nrho; ++i)
      logrho[i] = logrho[i]*ln10 + log(RHOGF);
#pragma omp parallel for simd if(ntemp>256)
    for (int i=0; i<ntemp; ++i)
      logtemp[i] *= ln10;                       // T already in MeV (no scale)

    // Table: loop over all grid points once, convert in‑place
    const CCTK_REAL cs2_scale = LENGTHGF*LENGTHGF / (TIMEGF*TIMEGF);
#pragma omp parallel for if(npoints>1024)
    for (size_t idx=0; idx<npoints; ++idx) {
      size_t b = idx*NTABLES;
      alltabl[b + PRESS ] = alltabl[b + PRESS ] * ln10 + log(PRESSGF);
      alltabl[b + EPS   ] = alltabl[b + EPS   ] * ln10 + log(EPSGF);
      alltabl[b + CS2   ] = alltabl[b + CS2   ] * cs2_scale;
      alltabl[b + DEDT  ] = alltabl[b + DEDT  ] * EPSGF;
      alltabl[b + DPDRHOE] = alltabl[b + DPDRHOE] * PRESSGF / RHOGF;
      alltabl[b + DPDERHO] = alltabl[b + DPDERHO] * PRESSGF / EPSGF;
      // remaining datasets are dimensionless → no scaling needed
    }

    // -------------------------------------------------------
    // 4.  Build interpolator (managed memory → accessible on GPU)
    // -------------------------------------------------------
    std::array<size_t,3> dims{static_cast<size_t>(nrho),
                              static_cast<size_t>(ntemp),
                              static_cast<size_t>(nye)};
    interptable = reinterpret_cast<decltype(interptable)>(
        amrex::The_Managed_Arena()->alloc(sizeof(*interptable)));
    new (interptable) linear_interp_uniform_ND_t<CCTK_REAL,3,NTABLES>(
            alltabl, dims, logrho, logtemp, yes);

    // -------------------------------------------------------
    // 5.  Physical ranges
    // -------------------------------------------------------
    set_range_rho ( range{ exp(interptable->xmin<0>()), exp(interptable->xmax<0>()) } );
    set_range_temp( range{ exp(interptable->xmin<1>()), exp(interptable->xmax<1>()) } );
    set_range_ye  ( range{        interptable->xmin<2>(),        interptable->xmax<2>() } );

    rgeps = compute_eps_range_full_table();
    out_rgeps = rgeps;
  } // init

  // =============================================================================
  //  Device‑/host‑callable lookup helpers (always inlined)
  // =============================================================================

  // ---------------------------------------------------------------------------
  //  logT(ρ, ε, Ye)  –  invert table in ε to obtain T (all log‑space)
  // ---------------------------------------------------------------------------
  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  logtemp_from_valid_rho_eps_ye(CCTK_REAL rho, CCTK_REAL &eps, CCTK_REAL ye) const {
    rho = clamp(rho,  rgrho.min,  rgrho.max);
    ye  = clamp(ye,   rgye.min,   rgye.max);

    // Shift ε (allows negative physical ε)
    eps = clamp(eps, rgeps.min, rgeps.max);
    const CCTK_REAL lrho = log(rho);
    const CCTK_REAL leps = log(eps + *energy_shift);

    // Quick exit at table edges
    const CCTK_REAL eps_min = interptable->interpolate<EV::EPS>(lrho, interptable->xmin<1>(), ye)[0];
    const CCTK_REAL eps_max = interptable->interpolate<EV::EPS>(lrho, interptable->xmax<1>(), ye)[0];
    if (leps <= eps_min) { eps = exp(eps_min) - *energy_shift; return interptable->xmin<1>(); }
    if (leps >= eps_max) { eps = exp(eps_max) - *energy_shift; return interptable->xmax<1>(); }

    // Monotonic root: f(lt) = leps_target - logε_table(lt)
    auto f = [&](CCTK_REAL &lt){ return leps - interptable->interpolate<EV::EPS>(lrho, lt, ye)[0]; };
    return zero_brent(interptable->xmin<1>(), interptable->xmax<1>(), 1e-14, f);
  }

  // ---------------------------------------------------------------------------
  //  P(ρ,T,Ye)  –  forward lookup
  // ---------------------------------------------------------------------------
  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  press_from_valid_rho_temp_ye(CCTK_REAL rho, CCTK_REAL temp, CCTK_REAL ye) const {
    rho = clamp(rho,  rgrho.min,  rgrho.max);
    temp= clamp(temp, rgtemp.min, rgtemp.max);
    ye  = clamp(ye,   rgye.min,   rgye.max);
    const CCTK_REAL v = interptable->interpolate<EV::PRESS>(log(rho), log(temp), ye)[0];
    return exp(v);
  }

  // ---------------------------------------------------------------------------
  //  P(ρ,ε,Ye)  –  via inversion ε→T then forward
  // ---------------------------------------------------------------------------
  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  press_from_valid_rho_eps_ye(CCTK_REAL rho, CCTK_REAL &eps, CCTK_REAL ye) const {
    rho = clamp(rho, rgrho.min, rgrho.max);
    ye  = clamp(ye,  rgye.min,  rgye.max);
    const CCTK_REAL lt = logtemp_from_valid_rho_eps_ye(rho, eps, ye);
    const CCTK_REAL v  = interptable->interpolate<EV::PRESS>(log(rho), lt, ye)[0];
    return exp(v);
  }

  // ---------------------------------------------------------------------------
  //  ε(ρ,T,Ye)
  // ---------------------------------------------------------------------------
  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  eps_from_valid_rho_temp_ye(CCTK_REAL rho, CCTK_REAL temp, CCTK_REAL ye) const {
    rho  = clamp(rho,  rgrho.min,  rgrho.max);
    temp = clamp(temp, rgtemp.min, rgtemp.max);
    ye   = clamp(ye,   rgye.min,   rgye.max);
    const CCTK_REAL v = interptable->interpolate<EV::EPS>(log(rho), log(temp), ye)[0];
    return exp(v) - *energy_shift;
  }

  // ---------------------------------------------------------------------------
  //  c_s(ρ,T,Ye)   &   c_s(ρ,ε,Ye)
  // ---------------------------------------------------------------------------
  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  csnd_from_valid_rho_temp_ye(CCTK_REAL rho, CCTK_REAL temp, CCTK_REAL ye) const {
    rho  = clamp(rho,  rgrho.min,  rgrho.max);
    temp = clamp(temp, rgtemp.min, rgtemp.max);
    ye   = clamp(ye,   rgye.min,   rgye.max);
    const CCTK_REAL v = interptable->interpolate<EV::CS2>(log(rho), log(temp), ye)[0];
    assert(v>=0);
    return sqrt(v);
  }

  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  csnd_from_valid_rho_eps_ye(CCTK_REAL rho, CCTK_REAL &eps, CCTK_REAL ye) const {
    const CCTK_REAL lt = logtemp_from_valid_rho_eps_ye(rho, eps, ye);
    const CCTK_REAL v  = interptable->interpolate<EV::CS2>(log(clamp(rho,rgrho.min,rgrho.max)), lt, clamp(ye,rgye.min,rgye.max))[0];
    assert(v>=0);
    return sqrt(v);
  }

  // ---------------------------------------------------------------------------
  //  T(ρ,ε,Ye)
  // ---------------------------------------------------------------------------
  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  temp_from_valid_rho_eps_ye(CCTK_REAL rho, CCTK_REAL &eps, CCTK_REAL ye) const {
    return exp( logtemp_from_valid_rho_eps_ye(rho, eps, ye) );
  }

  // ---------------------------------------------------------------------------
  //  Entropy s(ρ,T,Ye)
  // ---------------------------------------------------------------------------
  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  entropy_from_valid_rho_temp_ye(CCTK_REAL rho, CCTK_REAL temp, CCTK_REAL ye) const {
    rho  = clamp(rho,  rgrho.min,  rgrho.max);
    temp = clamp(temp, rgtemp.min, rgtemp.max);
    ye   = clamp(ye,   rgye.min,   rgye.max);
    return interptable->interpolate<EV::S>(log(rho), log(temp), ye)[0];
  }

  // ---------------------------------------------------------------------------
  //  κ(ρ,ε,Ye) — for tabulated EOS κ ≡ s (entropy)
  // ---------------------------------------------------------------------------
  CCTK_HOST CCTK_DEVICE inline CCTK_REAL
  kappa_from_valid_rho_eps_ye(CCTK_REAL rho, CCTK_REAL &eps, CCTK_REAL ye) const {
    return entropy_from_valid_rho_temp_ye(rho, temp_from_valid_rho_eps_ye(rho, eps, ye), ye);
  }

  // ---------------------------------------------------------------------------
  //  ε range at fixed (ρ,Ye)     — rarely needed, but handy for fallback logic
  // ---------------------------------------------------------------------------
  CCTK_HOST CCTK_DEVICE inline range
  range_eps_from_valid_rho_ye(CCTK_REAL rho, CCTK_REAL ye) const {
    rho = clamp(rho, rgrho.min, rgrho.max);
    ye  = clamp(ye,  rgye.min,  rgye.max);
    const CCTK_REAL vmin = interptable->interpolate<EV::EPS>(log(rho), interptable->xmin<1>(), ye)[0];
    const CCTK_REAL vmax = interptable->interpolate<EV::EPS>(log(rho), interptable->xmax<1>(), ye)[0];
    return { exp(vmin) - *energy_shift, exp(vmax) - *energy_shift };
  }

  // ---------------------------------------------------------------------------
  //  Global ε range across full table (used during init)
  // ---------------------------------------------------------------------------
  CCTK_HOST CCTK_DEVICE inline range compute_eps_range_full_table() const {
    const size_t n0 = interptable->num_points[0];
    const size_t n1 = interptable->num_points[1];
    const size_t n2 = interptable->num_points[2];
    const size_t total = n0*n1*n2;
    const CCTK_REAL *eps_log = interptable->y + EV::EPS*total;

    CCTK_REAL eps_min =  std::numeric_limits<CCTK_REAL>::max();
    CCTK_REAL eps_max = -std::numeric_limits<CCTK_REAL>::max();
#pragma omp parallel for reduction(min:eps_min) reduction(max:eps_max)
    for (size_t i=0; i<total; ++i) {
      CCTK_REAL val = exp(eps_log[i]) - *energy_shift;
      eps_min = (val < eps_min ? val : eps_min);
      eps_max = (val > eps_max ? val : eps_max);
    }
    return { eps_min, eps_max };
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_valid_rho_press_ye(const CCTK_REAL rho, const CCTK_REAL press,
                              const CCTK_REAL ye) const {

    printf(
        "This routine should not be used. There is no monotonicity condition "
        "to enforce a succesfull inversion from eps(press). So you better "
        "rewrite your code to not require this call. \n");
    assert(false);
    return 0;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  press_derivs_from_valid_rho_eps_ye(CCTK_REAL &press, CCTK_REAL &dpdrho,
                                     CCTK_REAL &dpdeps, const CCTK_REAL rho,
                                     const CCTK_REAL eps,
                                     const CCTK_REAL ye) const {
    printf("press_derivs_from_valid_rho_eps_ye is not supported for now! \n");
    assert(false);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  press_from_valid_rho_kappa_ye(const CCTK_REAL rho,
                                const CCTK_REAL kappa, // kappa=entropy
                                const CCTK_REAL ye) const {
    printf(
        "press_from_valid_rho_kappa_ye is not supported for tabulated EOS! \n");
    assert(false);
    return 0.0;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_valid_rho_kappa_ye(const CCTK_REAL rho,
                              const CCTK_REAL kappa, // kappa=entropy
                              const CCTK_REAL ye) const {
    printf(
        "eps_from_valid_rho_kappa_ye is not supported for tabulated EOS! \n");
    assert(false);
    return 0.0;
  }
  
}; // class eos_3p_tabulated3d

} // namespace EOSX

#endif // EOS_3P_TABULATED3D_HXX

