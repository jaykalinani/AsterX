#include "pdesolvers.hxx"

#ifdef _OPENMP
#include <omp.h>
#else
static inline int omp_get_max_threads() { return 1; }
static inline int omp_get_thread_num() { return 0; }
#endif

#include <cassert>

namespace PDESolvers {

bool csr_t::invariant() const {
  assert(m >= 0);
  assert(n >= 0);
  assert(int(rowptrs.size()) == m + 1);
  assert(int(colvals.size()) == rowptrs.back());
  assert(int(nzvals.size()) == rowptrs.back());
  assert(rowptrs.at(0) == 0);
  for (int i = 0; i < m; ++i) {
    assert(rowptrs.at(i) <= rowptrs.at(i + 1));
    for (int jp = rowptrs.at(i); jp < rowptrs.at(i + 1); ++jp) {
      const int j = colvals.at(jp);
      assert(j >= 0 && j < n);
      // We could check that j are sorted
    }
  }
  return true;
}

csr_t::csr_t(
    const int m, const int n,
    const std::vector<
        std::vector<std::vector<std::tuple<int, int, CCTK_REAL> > > > &values)
    : m(m), n(n) {
  for (const auto &values1 : values) {
    for (const auto &values2 : values1) {
      for (const auto &ijv : values2) {
        const int i = std::get<0>(ijv);
        const int j = std::get<1>(ijv);
        const CCTK_REAL v = std::get<2>(ijv);
        assert(i >= 0 && i < m);
        assert(j >= 0 && j < n);
        assert(i >= int(rowptrs.size()) - 1);
        while (int(rowptrs.size()) <= i)
          rowptrs.push_back(colvals.size());
        colvals.push_back(j);
        nzvals.push_back(v);
      }
    }
  }
  assert(int(rowptrs.size()) <= m);
  while (int(rowptrs.size()) <= m)
    rowptrs.push_back(colvals.size());
  assert(invariant());
}

////////////////////////////////////////////////////////////////////////////////

void jacobian_t::clear() { entries.clear(); }

void jacobian_t::set_matrix_entries(const csr_t &Jp, Mat J) const {
  std::vector<CCTK_REAL> values;
  for (const auto &e : entries) {
    const auto i = std::get<0>(e);
    const auto j = std::get<1>(e);
    const auto v = std::get<2>(e);
    assert(i >= 0);
    if (j < 0) {
      // ignore this point
    } else if (j >= 0 && j < prolongation_index_offset) {
      // regular point
      MatSetValue(J, i, j, v, ADD_VALUES);
    } else {
      // prolongated point
      const int row = j - prolongation_index_offset;
      assert(0 <= row && row < Jp.m);
      const int rowptr0 = Jp.rowptrs.at(row);
      const int rowptr1 = Jp.rowptrs.at(row + 1);
      values.resize(rowptr1 - rowptr0);
      for (int colindp = 0; colindp < rowptr1 - rowptr0; ++colindp)
        values.at(colindp) = v * Jp.nzvals.at(rowptr0 + colindp);
      MatSetValues(J, 1, &i, rowptr1 - rowptr0, &Jp.colvals.at(rowptr0),
                   values.data(), ADD_VALUES);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

jacobians_t::jacobians_t() : jacobians(omp_get_max_threads()) {}

jacobian_t &jacobians_t::get_local() {
  return jacobians.at(omp_get_thread_num());
}

void jacobians_t::clear() {
  for (auto &j : jacobians)
    j.clear();
}

void jacobians_t::define_matrix(const csr_t &Jp, Mat J) const {
  PetscErrorCode ierr;
  // TODO: Call MatMPIAIJSetPreallocation again, with good estimates.
  // Or better call MatMPIAIJSetPreallocationCSR?
  ierr = MatZeroEntries(J);
  assert(!ierr);
  for (const auto &j : jacobians)
    j.set_matrix_entries(Jp, J);
  ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  assert(!ierr);
  ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
  assert(!ierr);
}

////////////////////////////////////////////////////////////////////////////////

std::optional<jacobians_t> jacobians;

} // namespace PDESolvers
