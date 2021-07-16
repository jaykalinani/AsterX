#include "pdesolvers.hxx"

#ifdef _OPENMP
#include <omp.h>
#else
static inline int omp_get_max_threads() { return 1; }
static inline int omp_get_num_threads() { return 1; }
static inline int omp_get_thread_num() { return 0; }
#endif

#include <cassert>
#include <tuple>

namespace PDESolvers {

void jacobian_t::clear() { entries.clear(); }

void jacobian_t::set_matrix_entries(Mat J) const {
  for (const auto &e : entries)
    MatSetValue(J, std::get<0>(e), std::get<1>(e), std::get<2>(e), ADD_VALUES);
}

////////////////////////////////////////////////////////////////////////////////

jacobians_t::jacobians_t() : jacobians(omp_get_max_threads()) {}

void jacobians_t::clear() {
  for (auto &j : jacobians)
    j.clear();
}

jacobian_t &jacobians_t::get_local() {
  return jacobians.at(omp_get_thread_num());
}

void jacobians_t::define_matrix(Mat J) const {
  PetscErrorCode ierr;
  ierr = MatZeroEntries(J);
  assert(!ierr);
  for (const auto &j : jacobians)
    j.set_matrix_entries(J);
  ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  assert(!ierr);
  ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
  assert(!ierr);
}

////////////////////////////////////////////////////////////////////////////////

std::optional<jacobians_t> jacobians;

} // namespace PDESolvers
