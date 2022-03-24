
#ifndef MULTIPATCH_CAKE_CARTESIAN_JAC_HXX
#define MULTIPATCH_CAKE_CARTESIAN_JAC_HXX

#include "cake.hxx"

namespace MultiPatch {
namespace Cake {

template <typename T>
inline std::tuple<jac_t, djac_t>
cake_cartesian_jac(const PatchTransformations &, const svec_u &,
                   const jacobian_data<T> &) {
  jac_t J = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
  djac_t dJ{};

  return std::make_tuple(J, dJ);
}

} // namespace Cake
} // namespace MultiPatch

#endif // MULTIPATCH_CAKE_CARTESIAN_JAC_HXX
