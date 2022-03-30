
#ifndef MULTIPATCH_CAKE_PLUS_Y_JAC_HXX
#define MULTIPATCH_CAKE_PLUS_Y_JAC_HXX

#include "cake.hxx"

namespace MultiPatch {
namespace Cake {

inline std::tuple<jac_t, djac_t> cake_plus_y_jac(const PatchTransformations &pt,
                                                 const svec_u &local_vars) {
  jac_t J{};
  djac_t dJ{};

  const auto r0 = pt.cake_inner_boundary_radius;
  const auto r1 = pt.cake_outer_boundary_radius;

  const auto global_vars =
      pt.local2global(pt, static_cast<int>(patch_piece::plus_y), local_vars);
  const auto x = global_vars(0);
  const auto y = global_vars(1);
  const auto z = global_vars(2);

  J(0)(0) = 0;
  J(0)(1) = -(z / Power(y, 2));
  J(0)(2) = 1 / y;
  J(1)(0) = -(1 / y);
  J(1)(1) = x / Power(y, 2);
  J(1)(2) = 0;
  J(2)
  (0) = (2 * x *
         (1 + (2 * r0 * (r0 - r1) + Power(x, 2) + Power(z, 2)) /
                  Sqrt(4 * Power(r1, 2) * Power(y, 2) +
                       Power(Power(x, 2) + Power(z, 2), 2) +
                       4 * Power(r0, 2) *
                           (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
                       4 * r0 * r1 *
                           (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))))) /
        Power(r0 - r1, 2);
  J(2)
  (1) = (4 * y) /
        Sqrt(4 * Power(r1, 2) * Power(y, 2) +
             Power(Power(x, 2) + Power(z, 2), 2) +
             4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
             4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2)));
  J(2)
  (2) = (2 * z *
         (1 + (2 * r0 * (r0 - r1) + Power(x, 2) + Power(z, 2)) /
                  Sqrt(4 * Power(r1, 2) * Power(y, 2) +
                       Power(Power(x, 2) + Power(z, 2), 2) +
                       4 * Power(r0, 2) *
                           (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
                       4 * r0 * r1 *
                           (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))))) /
        Power(r0 - r1, 2);
  dJ(0)(0, 0) = 0;
  dJ(0)(0, 1) = 0;
  dJ(0)(0, 2) = 0;
  dJ(0)(1, 0) = 0;
  dJ(0)(1, 1) = (2 * z) / Power(y, 3);
  dJ(0)(1, 2) = -Power(y, -2);
  dJ(0)(2, 0) = 0;
  dJ(0)(2, 1) = -Power(y, -2);
  dJ(0)(2, 2) = 0;
  dJ(1)(0, 0) = 0;
  dJ(1)(0, 1) = Power(y, -2);
  dJ(1)(0, 2) = 0;
  dJ(1)(1, 0) = Power(y, -2);
  dJ(1)(1, 1) = (-2 * x) / Power(y, 3);
  dJ(1)(1, 2) = 0;
  dJ(1)(2, 0) = 0;
  dJ(1)(2, 1) = 0;
  dJ(1)(2, 2) = 0;
  dJ(2)(0, 0) =
      (2 +
       (2 * (2 * r0 * (r0 - r1) + 3 * Power(x, 2) + Power(z, 2))) /
           Sqrt(4 * Power(r1, 2) * Power(y, 2) +
                Power(Power(x, 2) + Power(z, 2), 2) +
                4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
                4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))) -
       (4 * Power(x, 2) *
        Power(2 * r0 * (r0 - r1) + Power(x, 2) + Power(z, 2), 2)) /
           (Sqrt(4 * Power(r1, 2) * Power(y, 2) +
                 Power(Power(x, 2) + Power(z, 2), 2) +
                 4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
                 4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))) *
            (4 * Power(r1, 2) * Power(y, 2) +
             Power(Power(x, 2) + Power(z, 2), 2) +
             4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
             4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))))) /
      Power(r0 - r1, 2);
  dJ(2)(0, 1) =
      (-8 * x * y * (2 * r0 * (r0 - r1) + Power(x, 2) + Power(z, 2))) /
      (Sqrt(4 * Power(r1, 2) * Power(y, 2) +
            Power(Power(x, 2) + Power(z, 2), 2) +
            4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
            4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))) *
       (4 * Power(r1, 2) * Power(y, 2) + Power(Power(x, 2) + Power(z, 2), 2) +
        4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
        4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))));
  dJ(2)(0, 2) =
      (4 * x * z *
       (-Power(2 * r0 * (r0 - r1) + Power(x, 2) + Power(z, 2), 2) +
        (4 * Power(r1, 2) * Power(y, 2) + Power(Power(x, 2) + Power(z, 2), 2) +
         4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
         4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))))) /
      (Power(r0 - r1, 2) *
       Sqrt(4 * Power(r1, 2) * Power(y, 2) +
            Power(Power(x, 2) + Power(z, 2), 2) +
            4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
            4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))) *
       (4 * Power(r1, 2) * Power(y, 2) + Power(Power(x, 2) + Power(z, 2), 2) +
        4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
        4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))));
  dJ(2)(1, 0) =
      (-8 * x * y * (2 * r0 * (r0 - r1) + Power(x, 2) + Power(z, 2))) /
      (Sqrt(4 * Power(r1, 2) * Power(y, 2) +
            Power(Power(x, 2) + Power(z, 2), 2) +
            4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
            4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))) *
       (4 * Power(r1, 2) * Power(y, 2) + Power(Power(x, 2) + Power(z, 2), 2) +
        4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
        4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))));
  dJ(2)(1, 1) =
      (4 *
       (-4 * Power(r0 - r1, 2) * Power(y, 2) +
        (4 * Power(r1, 2) * Power(y, 2) + Power(Power(x, 2) + Power(z, 2), 2) +
         4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
         4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))))) /
      (Sqrt(4 * Power(r1, 2) * Power(y, 2) +
            Power(Power(x, 2) + Power(z, 2), 2) +
            4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
            4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))) *
       (4 * Power(r1, 2) * Power(y, 2) + Power(Power(x, 2) + Power(z, 2), 2) +
        4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
        4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))));
  dJ(2)(1, 2) =
      (-8 * y * z * (2 * r0 * (r0 - r1) + Power(x, 2) + Power(z, 2))) /
      (Sqrt(4 * Power(r1, 2) * Power(y, 2) +
            Power(Power(x, 2) + Power(z, 2), 2) +
            4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
            4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))) *
       (4 * Power(r1, 2) * Power(y, 2) + Power(Power(x, 2) + Power(z, 2), 2) +
        4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
        4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))));
  dJ(2)(2, 0) =
      (4 * x * z *
       (-Power(2 * r0 * (r0 - r1) + Power(x, 2) + Power(z, 2), 2) +
        (4 * Power(r1, 2) * Power(y, 2) + Power(Power(x, 2) + Power(z, 2), 2) +
         4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
         4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))))) /
      (Power(r0 - r1, 2) *
       Sqrt(4 * Power(r1, 2) * Power(y, 2) +
            Power(Power(x, 2) + Power(z, 2), 2) +
            4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
            4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))) *
       (4 * Power(r1, 2) * Power(y, 2) + Power(Power(x, 2) + Power(z, 2), 2) +
        4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
        4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))));
  dJ(2)(2, 1) =
      (-8 * y * z * (2 * r0 * (r0 - r1) + Power(x, 2) + Power(z, 2))) /
      (Sqrt(4 * Power(r1, 2) * Power(y, 2) +
            Power(Power(x, 2) + Power(z, 2), 2) +
            4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
            4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))) *
       (4 * Power(r1, 2) * Power(y, 2) + Power(Power(x, 2) + Power(z, 2), 2) +
        4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
        4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))));
  dJ(2)(2, 2) =
      (2 +
       (2 * (2 * r0 * (r0 - r1) + Power(x, 2) + 3 * Power(z, 2))) /
           Sqrt(4 * Power(r1, 2) * Power(y, 2) +
                Power(Power(x, 2) + Power(z, 2), 2) +
                4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
                4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))) -
       (4 * Power(z, 2) *
        Power(2 * r0 * (r0 - r1) + Power(x, 2) + Power(z, 2), 2)) /
           (Sqrt(4 * Power(r1, 2) * Power(y, 2) +
                 Power(Power(x, 2) + Power(z, 2), 2) +
                 4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
                 4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))) *
            (4 * Power(r1, 2) * Power(y, 2) +
             Power(Power(x, 2) + Power(z, 2), 2) +
             4 * Power(r0, 2) * (Power(x, 2) + Power(y, 2) + Power(z, 2)) -
             4 * r0 * r1 * (Power(x, 2) + 2 * Power(y, 2) + Power(z, 2))))) /
      Power(r0 - r1, 2);

  return std::make_tuple(J, dJ);
}

} // namespace Cake
} // namespace MultiPatch

#endif // MULTIPATCH_CAKE_PLUS_Y_JAC_HXX
