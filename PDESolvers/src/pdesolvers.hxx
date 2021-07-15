#ifndef PDESOLVERS_HXX
#define PDESOLVERS_HXX

#include <fixmath.hxx>
#include <cctk.h>

#include <petsc.h>

#include <cassert>
#include <optional>

namespace PDESolvers {

class jacobian_t {
  Mat J;

public:
  jacobian_t(const jacobian_t &) = default;
  jacobian_t(jacobian_t &&) = default;
  jacobian_t &operator=(const jacobian_t &) = default;
  jacobian_t &operator=(jacobian_t &&) = default;

  jacobian_t(Mat J_) : J(J_) {}
  void add_value(int i, int j, CCTK_REAL v) const {
#pragma omp critical(PDESolvers_add_value)
    MatSetValue(J, i, j, v, ADD_VALUES);
  }
};

extern std::optional<jacobian_t> jacobian;
} // namespace PDESolvers

#endif // #ifndef PDESOLVERS_HXX
