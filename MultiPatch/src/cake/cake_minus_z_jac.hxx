
#ifndef MULTIPATCH_CAKE_MINUS_Z_JAC_HXX
#define MULTIPATCH_CAKE_MINUS_Z_JAC_HXX

#include "cake.hxx"

namespace MultiPatch {
namespace Cake {

inline std::tuple<jac_t, djac_t>
cake_minus_z_jac(const PatchTransformations &pt, const svec_u &local_vars) {
  jac_t J{};
  djac_t dJ{};

  const auto r0 = pt.cake_inner_boundary_radius;
  const auto r1 = pt.cake_outer_boundary_radius;

  const auto global_vars =
      pt.local2global(pt, static_cast<int>(patch_piece::minus_z), local_vars);
  const auto x = global_vars(0);
  const auto y = global_vars(1);
  const auto z = global_vars(2);

  J(0)(0) = -(1 / z);
  J(0)(1) = 0;
  J(0)(2) = x / Power(z, 2);
  J(1)(0) = 0;
  J(1)(1) = -(1 / z);
  J(1)(2) = y / Power(z, 2);
  J(2)
  (0) = (2 * x *
         (1 + (2 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) /
                  Sqrt((Power(x, 2) + Power(y, 2)) *
                           (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
                       4 * Power(r0 - r1, 2) * Power(z, 2)))) /
        Power(r0 - r1, 2);
  J(2)
  (1) = (2 * y *
         (1 + (2 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) /
                  Sqrt((Power(x, 2) + Power(y, 2)) *
                           (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
                       4 * Power(r0 - r1, 2) * Power(z, 2)))) /
        Power(r0 - r1, 2);
  J(2)
  (2) = (4 * z) / Sqrt((Power(x, 2) + Power(y, 2)) *
                           (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
                       4 * Power(r0 - r1, 2) * Power(z, 2));
  dJ(0)(0, 0) = 0;
  dJ(0)(0, 1) = 0;
  dJ(0)(0, 2) = Power(z, -2);
  dJ(0)(1, 0) = 0;
  dJ(0)(1, 1) = 0;
  dJ(0)(1, 2) = 0;
  dJ(0)(2, 0) = Power(z, -2);
  dJ(0)(2, 1) = 0;
  dJ(0)(2, 2) = (-2 * x) / Power(z, 3);
  dJ(1)(0, 0) = 0;
  dJ(1)(0, 1) = 0;
  dJ(1)(0, 2) = 0;
  dJ(1)(1, 0) = 0;
  dJ(1)(1, 1) = 0;
  dJ(1)(1, 2) = Power(z, -2);
  dJ(1)(2, 0) = 0;
  dJ(1)(2, 1) = Power(z, -2);
  dJ(1)(2, 2) = (-2 * y) / Power(z, 3);
  dJ(2)(0, 0) =
      (2 +
       (2 * (2 * r0 * (r0 - r1) + 3 * Power(x, 2) + Power(y, 2))) /
           Sqrt((Power(x, 2) + Power(y, 2)) *
                    (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
                4 * Power(r0 - r1, 2) * Power(z, 2)) -
       (4 * Power(x, 2) *
        Power(2 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2), 2)) /
           (Sqrt((Power(x, 2) + Power(y, 2)) *
                     (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
                 4 * Power(r0 - r1, 2) * Power(z, 2)) *
            ((Power(x, 2) + Power(y, 2)) *
                 (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
             4 * Power(r0 - r1, 2) * Power(z, 2)))) /
      Power(r0 - r1, 2);
  dJ(2)(0, 1) = (4 * x * y *
                 (-Power(2 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2), 2) +
                  ((Power(x, 2) + Power(y, 2)) *
                       (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
                   4 * Power(r0 - r1, 2) * Power(z, 2)))) /
                (Power(r0 - r1, 2) *
                 Sqrt((Power(x, 2) + Power(y, 2)) *
                          (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
                      4 * Power(r0 - r1, 2) * Power(z, 2)) *
                 ((Power(x, 2) + Power(y, 2)) *
                      (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
                  4 * Power(r0 - r1, 2) * Power(z, 2)));
  dJ(2)(0, 2) =
      (-8 * x * (2 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) * z) /
      (Sqrt((Power(x, 2) + Power(y, 2)) *
                (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
            4 * Power(r0 - r1, 2) * Power(z, 2)) *
       ((Power(x, 2) + Power(y, 2)) *
            (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
        4 * Power(r0 - r1, 2) * Power(z, 2)));
  dJ(2)(1, 0) = (4 * x * y *
                 (-Power(2 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2), 2) +
                  ((Power(x, 2) + Power(y, 2)) *
                       (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
                   4 * Power(r0 - r1, 2) * Power(z, 2)))) /
                (Power(r0 - r1, 2) *
                 Sqrt((Power(x, 2) + Power(y, 2)) *
                          (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
                      4 * Power(r0 - r1, 2) * Power(z, 2)) *
                 ((Power(x, 2) + Power(y, 2)) *
                      (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
                  4 * Power(r0 - r1, 2) * Power(z, 2)));
  dJ(2)(1, 1) =
      (2 +
       (2 * (2 * r0 * (r0 - r1) + Power(x, 2) + 3 * Power(y, 2))) /
           Sqrt((Power(x, 2) + Power(y, 2)) *
                    (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
                4 * Power(r0 - r1, 2) * Power(z, 2)) -
       (4 * Power(y, 2) *
        Power(2 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2), 2)) /
           (Sqrt((Power(x, 2) + Power(y, 2)) *
                     (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
                 4 * Power(r0 - r1, 2) * Power(z, 2)) *
            ((Power(x, 2) + Power(y, 2)) *
                 (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
             4 * Power(r0 - r1, 2) * Power(z, 2)))) /
      Power(r0 - r1, 2);
  dJ(2)(1, 2) =
      (-8 * y * (2 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) * z) /
      (Sqrt((Power(x, 2) + Power(y, 2)) *
                (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
            4 * Power(r0 - r1, 2) * Power(z, 2)) *
       ((Power(x, 2) + Power(y, 2)) *
            (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
        4 * Power(r0 - r1, 2) * Power(z, 2)));
  dJ(2)(2, 0) =
      (-8 * x * (2 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) * z) /
      (Sqrt((Power(x, 2) + Power(y, 2)) *
                (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
            4 * Power(r0 - r1, 2) * Power(z, 2)) *
       ((Power(x, 2) + Power(y, 2)) *
            (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
        4 * Power(r0 - r1, 2) * Power(z, 2)));
  dJ(2)(2, 1) =
      (-8 * y * (2 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) * z) /
      (Sqrt((Power(x, 2) + Power(y, 2)) *
                (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
            4 * Power(r0 - r1, 2) * Power(z, 2)) *
       ((Power(x, 2) + Power(y, 2)) *
            (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
        4 * Power(r0 - r1, 2) * Power(z, 2)));
  dJ(2)(2, 2) = (4 * (-4 * Power(r0 - r1, 2) * Power(z, 2) +
                      ((Power(x, 2) + Power(y, 2)) *
                           (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
                       4 * Power(r0 - r1, 2) * Power(z, 2)))) /
                (Sqrt((Power(x, 2) + Power(y, 2)) *
                          (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
                      4 * Power(r0 - r1, 2) * Power(z, 2)) *
                 ((Power(x, 2) + Power(y, 2)) *
                      (4 * r0 * (r0 - r1) + Power(x, 2) + Power(y, 2)) +
                  4 * Power(r0 - r1, 2) * Power(z, 2)));

  return std::make_tuple(J, dJ);
}

} // namespace Cake
} // namespace MultiPatch

#endif // MULTIPATCH_CAKE_MINUS_Z_JAC_HXX
