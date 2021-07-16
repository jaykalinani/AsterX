#ifndef PDESOLVERS_HXX
#define PDESOLVERS_HXX

#include <fixmath.hxx>
#include <cctk.h>

#include <petsc.h>

#include <optional>
#include <vector>

namespace PDESolvers {

class jacobian_t {
  std::vector<std::tuple<int, int, CCTK_REAL> > entries;

public:
  jacobian_t(const jacobian_t &) = default;
  jacobian_t(jacobian_t &&) = default;
  jacobian_t &operator=(const jacobian_t &) = default;
  jacobian_t &operator=(jacobian_t &&) = default;

  jacobian_t() = default;
  void clear();
  void add_value(int i, int j, CCTK_REAL v) { entries.emplace_back(i, j, v); }
  void set_matrix_entries(Mat J) const;
};

////////////////////////////////////////////////////////////////////////////////

class jacobians_t {
  std::vector<jacobian_t> jacobians;

public:
  jacobians_t(const jacobians_t &) = delete;
  jacobians_t(jacobians_t &&) = default;
  jacobians_t &operator=(const jacobians_t &) = delete;
  jacobians_t &operator=(jacobians_t &&) = default;

  jacobians_t();
  void clear();
  jacobian_t &get_local();
  void define_matrix(Mat J) const;
};

////////////////////////////////////////////////////////////////////////////////

extern std::optional<jacobians_t> jacobians;

} // namespace PDESolvers

#endif // #ifndef PDESOLVERS_HXX
