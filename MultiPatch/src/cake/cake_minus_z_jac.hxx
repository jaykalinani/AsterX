
#ifndef MULTIPATCH_CAKE_MINUS_Z_JAC_HXX
#define MULTIPATCH_CAKE_MINUS_Z_JAC_HXX

#include "cake.hxx"

namespace MultiPatch {
namespace Cake {

template <typename T>
inline std::tuple<jac_t, djac_t>
cake_minus_z_jac(const PatchTransformations &pt, const svec_u &local_vars,
                 const jacobian_data<T> &data) {
  jac_t J{};
  djac_t dJ{};

  const auto a = local_vars(0);
  const auto b = local_vars(1);

  const auto core = data.core;
  const auto d_core_da = data.d_core_da;
  const auto d_core_db = data.d_core_db;
  const auto d_core_dc = data.d_core_dc;
  const auto d2_core_da_da = data.d2_core_da_da;
  const auto d2_core_da_db = data.d2_core_da_db;
  const auto d2_core_da_dc = data.d2_core_da_dc;
  const auto d2_core_db_db = data.d2_core_db_db;
  const auto d2_core_db_dc = data.d2_core_db_dc;
  const auto d2_core_dc_dc = data.d2_core_dc_dc;

  J(0)(0) = core + a * d_core_da;
  J(0)(1) = a * d_core_db;
  J(0)(2) = a * d_core_dc;
  J(1)(0) = b * d_core_da;
  J(1)(1) = core + b * d_core_db;
  J(1)(2) = b * d_core_dc;
  J(2)(0) = -d_core_da;
  J(2)(1) = -d_core_db;
  J(2)(2) = -d_core_dc;
  dJ(0)(0, 0) = 2 * d_core_da + a * d2_core_da_da;
  dJ(0)(0, 1) = d_core_db + a * d2_core_da_db;
  dJ(0)(0, 2) = d_core_dc + a * d2_core_da_dc;
  dJ(0)(1, 0) = d_core_db + a * d2_core_da_db;
  dJ(0)(1, 1) = a * d2_core_db_db;
  dJ(0)(1, 2) = a * d2_core_db_dc;
  dJ(0)(2, 0) = d_core_dc + a * d2_core_da_dc;
  dJ(0)(2, 1) = a * d2_core_db_dc;
  dJ(0)(2, 2) = a * d2_core_dc_dc;
  dJ(1)(0, 0) = b * d2_core_da_da;
  dJ(1)(0, 1) = d_core_da + b * d2_core_da_db;
  dJ(1)(0, 2) = b * d2_core_da_dc;
  dJ(1)(1, 0) = d_core_da + b * d2_core_da_db;
  dJ(1)(1, 1) = 2 * d_core_db + b * d2_core_db_db;
  dJ(1)(1, 2) = d_core_dc + b * d2_core_db_dc;
  dJ(1)(2, 0) = b * d2_core_da_dc;
  dJ(1)(2, 1) = d_core_dc + b * d2_core_db_dc;
  dJ(1)(2, 2) = b * d2_core_dc_dc;
  dJ(2)(0, 0) = -d2_core_da_da;
  dJ(2)(0, 1) = -d2_core_da_db;
  dJ(2)(0, 2) = -d2_core_da_dc;
  dJ(2)(1, 0) = -d2_core_da_db;
  dJ(2)(1, 1) = -d2_core_db_db;
  dJ(2)(1, 2) = -d2_core_db_dc;
  dJ(2)(2, 0) = -d2_core_da_dc;
  dJ(2)(2, 1) = -d2_core_db_dc;
  dJ(2)(2, 2) = -d2_core_dc_dc;

  return std::make_tuple(J, dJ);
}

} // namespace Cake
} // namespace MultiPatch

#endif // MULTIPATCH_CAKE_MINUS_Z_JAC_HXX
