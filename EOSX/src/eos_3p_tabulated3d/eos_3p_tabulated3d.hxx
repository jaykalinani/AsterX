#ifndef EOS_3P_TABULATED3D_HXX
#define EOS_3P_TABULATED3D_HXX

#include <cctk.h>
#include <cmath>
#include <hdf5.h>
#include <mpi.h>
#include <string>
#include "../eos_3p.hxx"
#include "../utils/eos_brent.hxx"
#include "../utils/eos_linear_interp_ND.hxx"

#define NTABLES 19

namespace EOSX {
using namespace std;
using namespace eos_constants;

class eos_3p_tabulated3d : public eos_3p {

public:
  enum errors {
    NO_ERRORS = 0,
    RHO_TOO_HIGH,
    RHO_TOO_LOW,
    YE_TOO_HIGH,
    YE_TOO_LOW,
    TEMP_TOO_HIGH,
    TEMP_TOO_LOW,
    num_errors
  };

  enum EV {
    PRESS = 0,
    EPS,
    S,
    CS2,
    MUE,
    MUP,
    MUN,
    XA,
    XH,
    XN,
    XP,
    ABAR,
    ZBAR,
    NUM_VARS
  };

  CCTK_REAL gamma;
  range rgeps;
  CCTK_REAL *energy_shift;
  linear_interp_uniform_ND_t<CCTK_REAL, 3, NTABLES> *interptable;

  CCTK_HOST void
  init(const string &filename, range &rgeps_, const range &rgrho_, const range &rgye_) {
    eos_readtable_scollapse(filename);
    set_range_rho(rgrho_);
    set_range_ye(rgye_);
    rgeps_ = range_eps_from_valid_rho_ye(rgrho_.min, rgye_.min);
    // TODO: first compute temp as a function of rho, ye, and eps, and then
    // initialize its range For now, as dummy, we pass range of eps as range of
    // temp
    set_range_temp(rgeps_);
  }

/*** Table readers ***/

// Routine reading the EOS table and filling the corresponding object
CCTK_HOST void
eos_readtable_scollapse(const string &filename) {

  CCTK_INT ntemp, nrho, nye;
  CCTK_VINFO("Reading Stellar Collapse EOS table '%s'", filename.c_str());

  auto fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  assert(fapl_id >= 0);
  hid_t file_id = 0;
  int rank_id;
  CHECK_ERROR(MPI_Comm_rank(MPI_COMM_WORLD, &rank_id));

#ifdef H5_HAVE_PARALLEL
  CHECK_ERROR(H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL));
  CHECK_ERROR(H5Pset_all_coll_metadata_ops(fapl_id, true));
  file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl_id);
#else
  fapl_id = H5P_DEFAULT;
  if (rank_id == 0) {
    file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl_id);
  }
// TODO: bcast table info to all ranks
#endif

  assert(file_id >= 0);

  // Get number of points
  get_hdf5_int_dset(file_id, "pointstemp", 1, &ntemp);
  get_hdf5_int_dset(file_id, "pointsrho", 1, &nrho);
  get_hdf5_int_dset(file_id, "pointsye", 1, &nye);

  const int npoints = ntemp * nrho * nye;

  CCTK_VINFO("EOS table dimensions: ntemp = %d, nrho = %d, nye = %d", ntemp,
             nrho, nye);

  energy_shift = (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(sizeof(CCTK_REAL));
  assert(energy_shift);

  // Allocate memory for tables

  CCTK_REAL *logrho, *logtemp, *yes;
  CCTK_REAL *epstable;
  CCTK_REAL *alltables;

  if (!(alltables = (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(
            npoints * NTABLES * sizeof(CCTK_REAL)))) {
    printf("Cannot allocate memory for alltables");
  }
  if (!(logrho = (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(nrho *
                                                         sizeof(CCTK_REAL)))) {
    printf("Cannot allocate memory for logrho");
  }
  if (!(logtemp = (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(ntemp *
                                                          sizeof(CCTK_REAL)))) {
    printf("Cannot allocate memory for logtemp");
  }
  if (!(yes =
            (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(nye * sizeof(CCTK_REAL)))) {
    printf("Cannot allocate memory for yes");
  }

  // allocate epstable; a linear-scale eps table
  // that allows us to extrapolate to negative eps
  if (!(epstable = (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(
            npoints * sizeof(CCTK_REAL)))) {
    printf("Cannot allocate memory for epstable");
  }

  // Prepare HDF5 to read hyperslabs into alltables_temp
  hsize_t table_dims[2] = {NTABLES, (hsize_t)npoints};
  hsize_t var3[2] = {1, (hsize_t)npoints};
  hid_t mem3 = H5Screate_simple(2, table_dims, NULL);

  // hydro (and munu)
  get_hdf5_real_dset(file_id, "logpress", npoints,
                     &alltables[0 * npoints]);
  get_hdf5_real_dset(file_id, "logenergy", npoints,
                     &alltables[1 * npoints]);
  get_hdf5_real_dset(file_id, "entropy", npoints, &alltables[2 * npoints]);
  get_hdf5_real_dset(file_id, "munu", npoints, &alltables[3 * npoints]);
  get_hdf5_real_dset(file_id, "cs2", npoints, &alltables[4 * npoints]);
  get_hdf5_real_dset(file_id, "dedt", npoints, &alltables[5 * npoints]);
  get_hdf5_real_dset(file_id, "dpdrhoe", npoints, &alltables[6 * npoints]);
  get_hdf5_real_dset(file_id, "dpderho", npoints, &alltables[7 * npoints]);

  // chemical potentials
  get_hdf5_real_dset(file_id, "muhat", npoints, &alltables[8 * npoints]);
  get_hdf5_real_dset(file_id, "mu_e", npoints, &alltables[9 * npoints]);
  get_hdf5_real_dset(file_id, "mu_p", npoints, &alltables[10 * npoints]);
  get_hdf5_real_dset(file_id, "mu_n", npoints, &alltables[11 * npoints]);

  // compositions
  get_hdf5_real_dset(file_id, "Xa", npoints, &alltables[12 * npoints]);
  get_hdf5_real_dset(file_id, "Xh", npoints, &alltables[13 * npoints]);
  get_hdf5_real_dset(file_id, "Xn", npoints, &alltables[14 * npoints]);
  get_hdf5_real_dset(file_id, "Xp", npoints, &alltables[15 * npoints]);

  // average nucleus
  get_hdf5_real_dset(file_id, "Abar", npoints, &alltables[16 * npoints]);
  get_hdf5_real_dset(file_id, "Zbar", npoints, &alltables[17 * npoints]);

  // Gamma
  get_hdf5_real_dset(file_id, "gamma", npoints, &alltables[18 * npoints]);

  // Read additional tables and variables
  get_hdf5_real_dset(file_id, "logrho", nrho, logrho);
  get_hdf5_real_dset(file_id, "logtemp", ntemp, logtemp);
  get_hdf5_real_dset(file_id, "ye", nye, yes);
  get_hdf5_real_dset(file_id, "energy_shift", 1, energy_shift);

  CHECK_ERROR(H5Pclose(fapl_id));
  CHECK_ERROR(H5Sclose(mem3));

#ifdef H5_HAVE_PARALLEL
  CHECK_ERROR(H5Fclose(file_id));
#else
  if (rank_id == 0) {
    CHECK_ERROR(H5Fclose(file_id));
  }
#endif

  *energy_shift = *energy_shift * EPSGF;

  for (int i = 0; i < nrho; i++) {
    // rewrite:
    // logrho[i] = log(pow(10.0,logrho[i]) * RHOGF);
    // by using log(a^b*c) = b*log(a)+log(c)
    logrho[i] = logrho[i] * log(10.) + log(RHOGF);
  }

  for (int i = 0; i < ntemp; i++) {
    // logtemp[i] = log(pow(10.0,logtemp[i]));
    logtemp[i] = logtemp[i] * log(10.0);
  }

  // convert units
  for (int i = 0; i < npoints; i++) {

    { // pressure
      int idx = 0 + NTABLES * i;
      alltables[idx] = alltables[idx] * log(10.0) + log(PRESSGF);
    }

    { // eps
      int idx = 1 + NTABLES * i;
      alltables[idx] = alltables[idx] * log(10.0) + log(EPSGF);
      epstable[i] = exp(alltables[idx]);
    }

    { // cs2
      int idx = 4 + NTABLES * i;
      alltables[idx] *= LENGTHGF * LENGTHGF / TIMEGF / TIMEGF;
    }

    { // dedT
      int idx = 5 + NTABLES * i;
      alltables[idx] *= EPSGF;
    }

    { // dpdrhoe
      int idx = 6 + NTABLES * i;
      alltables[idx] *= PRESSGF / RHOGF;
    }

    { // dpderho
      int idx = 7 + NTABLES * i;
      alltables[idx] *= PRESSGF / EPSGF;
    }
  }

  auto num_points =
      std::array<size_t, 3>{size_t(nrho), size_t(ntemp), size_t(nye)};

  interptable = (linear_interp_uniform_ND_t<CCTK_REAL, 3, NTABLES> *)amrex::The_Managed_Arena()->alloc(
        sizeof *interptable);
  assert(interptable);
  new (interptable) linear_interp_uniform_ND_t<CCTK_REAL, 3, NTABLES> (alltables, num_points, logrho, logtemp, yes);
  // set up steps, mins, maxes here?
  return;
}

  
/*** Function calls to compute hydro variables ***/

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  logtemp_from_eps(const CCTK_REAL rho, CCTK_REAL &eps,
                   const CCTK_REAL ye) const {
    const auto lrho = log(rho);
    const auto leps = log(eps + *energy_shift);

    // Get eps ranges
    const auto varsmin =
        interptable->interpolate<EV::EPS>(lrho, interptable->xmin<1>(), ye);
    const auto varsmax =
        interptable->interpolate<EV::EPS>(lrho, interptable->xmax<1>(), ye);

    if (leps <= varsmin[0]) {
      eps = exp(varsmin[0]) - *energy_shift;
      return interptable->xmin<1>();
    }
    if (leps >= varsmax[0]) {
      eps = exp(varsmax[0]) - *energy_shift;
      return interptable->xmax<1>();
    }

    // Root finding interface closure
    const auto func = [&](CCTK_REAL &lt) {
      const auto vars = interptable->interpolate<EV::EPS>(lrho, lt, ye);
      return leps - vars[0];
    };

    return zero_brent(interptable->xmin<1>(), interptable->xmax<1>(), 1.e-14,
                      func);
  }
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  press_from_valid_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                               const CCTK_REAL ye) const {
    // error = checkbounds<true>(rho, temp, ye);
    const auto lrho = log(rho);
    const auto ltemp = log(temp);
    const auto vars = interptable->interpolate<EV::PRESS>(lrho, ltemp, ye);

    return exp(vars[0]);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  press_from_valid_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                              const CCTK_REAL ye) const {
    const CCTK_REAL lrho = log(rho);
    const CCTK_REAL ltemp = logtemp_from_eps(lrho, eps, ye);
    const auto vars = interptable->interpolate<EV::PRESS>(lrho, ltemp, ye);

    return exp(vars[0]);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_valid_rho_press_ye(const CCTK_REAL rho, const CCTK_REAL press,
                              const CCTK_REAL ye) const {

    assert(
        !"This routine should not be used. There is no monotonicity condition "
         "to enforce a succesfull inversion from eps(press). So you better "
         "rewrite your code to not require this call...");

    return 0;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_valid_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                             const CCTK_REAL ye) const {

    const CCTK_REAL lrho = log(rho);
    const CCTK_REAL ltemp = log(temp);
    auto const vars = interptable->interpolate<EV::EPS>(lrho, ltemp, ye);
    const auto eps = exp(vars[0]) - *energy_shift;
    return eps;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  csnd_from_valid_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                              const CCTK_REAL ye) const {
    const auto lrho = log(rho);
    const auto ltemp = log(temp);
    const auto vars = interptable->interpolate<EV::CS2>(lrho, ltemp, ye);
    assert(vars[0] >= 0); // Soundspeed^2 should never ever be negative

    return sqrt(vars[0]);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  csnd_from_valid_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                             const CCTK_REAL ye) const {
    const CCTK_REAL lrho = log(rho);
    const CCTK_REAL ltemp = logtemp_from_eps(lrho, eps, ye);
    const auto vars = interptable->interpolate<EV::CS2>(lrho, ltemp, ye);
    assert(vars[0] >= 0); // Soundspeed^2 should never ever be negative

    return sqrt(vars[0]);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  temp_from_valid_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                             const CCTK_REAL ye) const {
    const CCTK_REAL lrho = log(rho);
    const CCTK_REAL ltemp = logtemp_from_eps(lrho, eps, ye);

    return exp(ltemp);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  press_derivs_from_valid_rho_eps_ye(CCTK_REAL &press, CCTK_REAL &dpdrho,
                                     CCTK_REAL &dpdeps, const CCTK_REAL rho,
                                     const CCTK_REAL eps,
                                     const CCTK_REAL ye) const {
    printf("press_derivs_from_valid_rho_eps_ye is not supported for now!");
    exit(EXIT_FAILURE);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  entropy_from_valid_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                                 const CCTK_REAL ye) const {
    const CCTK_REAL lrho = log(rho);
    const CCTK_REAL ltemp = log(temp);
    const auto vars = interptable->interpolate<EV::S>(lrho, ltemp, ye);

    return vars[0];
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline range
  range_eps_from_valid_rho_ye(const CCTK_REAL rho, const CCTK_REAL ye) const {
    //    const CCTK_REAL lrho = log(rho);
    //    rgeps.min = interptable->interpolate<EV::EPS>(lrho,
    //    interptable->xmin<1>(), ye); rgeps.max =
    //    interptable->interpolate<EV::EPS>(lrho, interptable->xmax<1>(), ye);

    return rgeps;
  }
};
} // namespace EOSX

#endif
