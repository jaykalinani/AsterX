#ifndef MULTIPATCH_CAKE_JACOBIANS_HXX
#define MULTIPATCH_CAKE_JACOBIANS_HXX

#define CAKE_MINUS_X_JACOBIAN                                                  \
  J(0)(0) = -d_core_da;                                                        \
  J(0)(1) = -d_core_db;                                                        \
  J(0)(2) = -d_core_dc;                                                        \
  J(1)(0) = -(b * d_core_da);                                                  \
  J(1)(1) = -core - b * d_core_db;                                             \
  J(1)(2) = -(b * d_core_dc);                                                  \
  J(2)(0) = c * d_core_da;                                                     \
  J(2)(1) = c * d_core_db;                                                     \
  J(2)(2) = core + c * d_core_dc;

#define CAKE_MINUS_X_DJACOBIAN                                                 \
  dJ(0)(0, 0) = -d2_core_da_da;                                                \
  dJ(0)(0, 1) = -d2_core_da_db;                                                \
  dJ(0)(0, 2) = -d2_core_da_dc;                                                \
  dJ(0)(1, 0) = -d2_core_da_db;                                                \
  dJ(0)(1, 1) = -d2_core_db_db;                                                \
  dJ(0)(1, 2) = -d2_core_db_dc;                                                \
  dJ(0)(2, 0) = -d2_core_da_dc;                                                \
  dJ(0)(2, 1) = -d2_core_db_dc;                                                \
  dJ(0)(2, 2) = -d2_core_dc_dc;                                                \
  dJ(1)(0, 0) = -(b * d2_core_da_da);                                          \
  dJ(1)(0, 1) = -d_core_da - b * d2_core_da_db;                                \
  dJ(1)(0, 2) = -(b * d2_core_da_dc);                                          \
  dJ(1)(1, 0) = -d_core_da - b * d2_core_da_db;                                \
  dJ(1)(1, 1) = -2 * d_core_db - b * d2_core_db_db;                            \
  dJ(1)(1, 2) = -d_core_dc - b * d2_core_db_dc;                                \
  dJ(1)(2, 0) = -(b * d2_core_da_dc);                                          \
  dJ(1)(2, 1) = -d_core_dc - b * d2_core_db_dc;                                \
  dJ(1)(2, 2) = -(b * d2_core_dc_dc);                                          \
  dJ(2)(0, 0) = c * d2_core_da_da;                                             \
  dJ(2)(0, 1) = c * d2_core_da_db;                                             \
  dJ(2)(0, 2) = d_core_da + c * d2_core_da_dc;                                 \
  dJ(2)(1, 0) = c * d2_core_da_db;                                             \
  dJ(2)(1, 1) = c * d2_core_db_db;                                             \
  dJ(2)(1, 2) = d_core_db + c * d2_core_db_dc;                                 \
  dJ(2)(2, 0) = d_core_da + c * d2_core_da_dc;                                 \
  dJ(2)(2, 1) = d_core_db + c * d2_core_db_dc;                                 \
  dJ(2)(2, 2) = 2 * d_core_dc + c * d2_core_dc_dc;

#define CAKE_PLUS_X_JACOBIAN                                                   \
  J(0)(0) = d_core_da;                                                         \
  J(0)(1) = d_core_db;                                                         \
  J(0)(2) = d_core_dc;                                                         \
  J(1)(0) = b * d_core_da;                                                     \
  J(1)(1) = core + b * d_core_db;                                              \
  J(1)(2) = b * d_core_dc;                                                     \
  J(2)(0) = c * d_core_da;                                                     \
  J(2)(1) = c * d_core_db;                                                     \
  J(2)(2) = core + c * d_core_dc;

#define CAKE_PLUS_X_DJACOBIAN                                                  \
  dJ(0)(0, 0) = d2_core_da_da;                                                 \
  dJ(0)(0, 1) = d2_core_da_db;                                                 \
  dJ(0)(0, 2) = d2_core_da_dc;                                                 \
  dJ(0)(1, 0) = d2_core_da_db;                                                 \
  dJ(0)(1, 1) = d2_core_db_db;                                                 \
  dJ(0)(1, 2) = d2_core_db_dc;                                                 \
  dJ(0)(2, 0) = d2_core_da_dc;                                                 \
  dJ(0)(2, 1) = d2_core_db_dc;                                                 \
  dJ(0)(2, 2) = d2_core_dc_dc;                                                 \
  dJ(1)(0, 0) = b * d2_core_da_da;                                             \
  dJ(1)(0, 1) = d_core_da + b * d2_core_da_db;                                 \
  dJ(1)(0, 2) = b * d2_core_da_dc;                                             \
  dJ(1)(1, 0) = d_core_da + b * d2_core_da_db;                                 \
  dJ(1)(1, 1) = 2 * d_core_db + b * d2_core_db_db;                             \
  dJ(1)(1, 2) = d_core_dc + b * d2_core_db_dc;                                 \
  dJ(1)(2, 0) = b * d2_core_da_dc;                                             \
  dJ(1)(2, 1) = d_core_dc + b * d2_core_db_dc;                                 \
  dJ(1)(2, 2) = b * d2_core_dc_dc;                                             \
  dJ(2)(0, 0) = c * d2_core_da_da;                                             \
  dJ(2)(0, 1) = c * d2_core_da_db;                                             \
  dJ(2)(0, 2) = d_core_da + c * d2_core_da_dc;                                 \
  dJ(2)(1, 0) = c * d2_core_da_db;                                             \
  dJ(2)(1, 1) = c * d2_core_db_db;                                             \
  dJ(2)(1, 2) = d_core_db + c * d2_core_db_dc;                                 \
  dJ(2)(2, 0) = d_core_da + c * d2_core_da_dc;                                 \
  dJ(2)(2, 1) = d_core_db + c * d2_core_db_dc;                                 \
  dJ(2)(2, 2) = 2 * d_core_dc + c * d2_core_dc_dc;

#define CAKE_MINUS_Y_JACOBIAN                                                  \
  J(0)(0) = b * d_core_da;                                                     \
  J(0)(1) = core + b * d_core_db;                                              \
  J(0)(2) = b * d_core_dc;                                                     \
  J(1)(0) = -d_core_da;                                                        \
  J(1)(1) = -d_core_db;                                                        \
  J(1)(2) = -d_core_dc;                                                        \
  J(2)(0) = c * d_core_da;                                                     \
  J(2)(1) = c * d_core_db;                                                     \
  J(2)(2) = core + c * d_core_dc;

#define CAKE_MINUS_Y_DJACOBIAN                                                 \
  dJ(0)(0, 0) = b * d2_core_da_da;                                             \
  dJ(0)(0, 1) = d_core_da + b * d2_core_da_db;                                 \
  dJ(0)(0, 2) = b * d2_core_da_dc;                                             \
  dJ(0)(1, 0) = d_core_da + b * d2_core_da_db;                                 \
  dJ(0)(1, 1) = 2 * d_core_db + b * d2_core_db_db;                             \
  dJ(0)(1, 2) = d_core_dc + b * d2_core_db_dc;                                 \
  dJ(0)(2, 0) = b * d2_core_da_dc;                                             \
  dJ(0)(2, 1) = d_core_dc + b * d2_core_db_dc;                                 \
  dJ(0)(2, 2) = b * d2_core_dc_dc;                                             \
  dJ(1)(0, 0) = -d2_core_da_da;                                                \
  dJ(1)(0, 1) = -d2_core_da_db;                                                \
  dJ(1)(0, 2) = -d2_core_da_dc;                                                \
  dJ(1)(1, 0) = -d2_core_da_db;                                                \
  dJ(1)(1, 1) = -d2_core_db_db;                                                \
  dJ(1)(1, 2) = -d2_core_db_dc;                                                \
  dJ(1)(2, 0) = -d2_core_da_dc;                                                \
  dJ(1)(2, 1) = -d2_core_db_dc;                                                \
  dJ(1)(2, 2) = -d2_core_dc_dc;                                                \
  dJ(2)(0, 0) = c * d2_core_da_da;                                             \
  dJ(2)(0, 1) = c * d2_core_da_db;                                             \
  dJ(2)(0, 2) = d_core_da + c * d2_core_da_dc;                                 \
  dJ(2)(1, 0) = c * d2_core_da_db;                                             \
  dJ(2)(1, 1) = c * d2_core_db_db;                                             \
  dJ(2)(1, 2) = d_core_db + c * d2_core_db_dc;                                 \
  dJ(2)(2, 0) = d_core_da + c * d2_core_da_dc;                                 \
  dJ(2)(2, 1) = d_core_db + c * d2_core_db_dc;                                 \
  dJ(2)(2, 2) = 2 * d_core_dc + c * d2_core_dc_dc;

#define CAKE_PLUS_Y_JACOBIAN                                                   \
  J(0)(0) = -(b * d_core_da);                                                  \
  J(0)(1) = -core - b * d_core_db;                                             \
  J(0)(2) = -(b * d_core_dc);                                                  \
  J(1)(0) = d_core_da;                                                         \
  J(1)(1) = d_core_db;                                                         \
  J(1)(2) = d_core_dc;                                                         \
  J(2)(0) = c * d_core_da;                                                     \
  J(2)(1) = c * d_core_db;                                                     \
  J(2)(2) = core + c * d_core_dc;

#define CAKE_PLUS_Y_DJACOBIAN                                                  \
  dJ(0)(0, 0) = -(b * d2_core_da_da);                                          \
  dJ(0)(0, 1) = -d_core_da - b * d2_core_da_db;                                \
  dJ(0)(0, 2) = -(b * d2_core_da_dc);                                          \
  dJ(0)(1, 0) = -d_core_da - b * d2_core_da_db;                                \
  dJ(0)(1, 1) = -2 * d_core_db - b * d2_core_db_db;                            \
  dJ(0)(1, 2) = -d_core_dc - b * d2_core_db_dc;                                \
  dJ(0)(2, 0) = -(b * d2_core_da_dc);                                          \
  dJ(0)(2, 1) = -d_core_dc - b * d2_core_db_dc;                                \
  dJ(0)(2, 2) = -(b * d2_core_dc_dc);                                          \
  dJ(1)(0, 0) = d2_core_da_da;                                                 \
  dJ(1)(0, 1) = d2_core_da_db;                                                 \
  dJ(1)(0, 2) = d2_core_da_dc;                                                 \
  dJ(1)(1, 0) = d2_core_da_db;                                                 \
  dJ(1)(1, 1) = d2_core_db_db;                                                 \
  dJ(1)(1, 2) = d2_core_db_dc;                                                 \
  dJ(1)(2, 0) = d2_core_da_dc;                                                 \
  dJ(1)(2, 1) = d2_core_db_dc;                                                 \
  dJ(1)(2, 2) = d2_core_dc_dc;                                                 \
  dJ(2)(0, 0) = c * d2_core_da_da;                                             \
  dJ(2)(0, 1) = c * d2_core_da_db;                                             \
  dJ(2)(0, 2) = d_core_da + c * d2_core_da_dc;                                 \
  dJ(2)(1, 0) = c * d2_core_da_db;                                             \
  dJ(2)(1, 1) = c * d2_core_db_db;                                             \
  dJ(2)(1, 2) = d_core_db + c * d2_core_db_dc;                                 \
  dJ(2)(2, 0) = d_core_da + c * d2_core_da_dc;                                 \
  dJ(2)(2, 1) = d_core_db + c * d2_core_db_dc;                                 \
  dJ(2)(2, 2) = 2 * d_core_dc + c * d2_core_dc_dc;

#define CAKE_MINUS_Z_JACOBIAN                                                  \
  J(0)(0) = c * d_core_da;                                                     \
  J(0)(1) = c * d_core_db;                                                     \
  J(0)(2) = core + c * d_core_dc;                                              \
  J(1)(0) = b * d_core_da;                                                     \
  J(1)(1) = core + b * d_core_db;                                              \
  J(1)(2) = b * d_core_dc;                                                     \
  J(2)(0) = -d_core_da;                                                        \
  J(2)(1) = -d_core_db;                                                        \
  J(2)(2) = -d_core_dc;

#define CAKE_MINUS_Z_DJACOBIAN                                                 \
  dJ(0)(0, 0) = c * d2_core_da_da;                                             \
  dJ(0)(0, 1) = c * d2_core_da_db;                                             \
  dJ(0)(0, 2) = d_core_da + c * d2_core_da_dc;                                 \
  dJ(0)(1, 0) = c * d2_core_da_db;                                             \
  dJ(0)(1, 1) = c * d2_core_db_db;                                             \
  dJ(0)(1, 2) = d_core_db + c * d2_core_db_dc;                                 \
  dJ(0)(2, 0) = d_core_da + c * d2_core_da_dc;                                 \
  dJ(0)(2, 1) = d_core_db + c * d2_core_db_dc;                                 \
  dJ(0)(2, 2) = 2 * d_core_dc + c * d2_core_dc_dc;                             \
  dJ(1)(0, 0) = b * d2_core_da_da;                                             \
  dJ(1)(0, 1) = d_core_da + b * d2_core_da_db;                                 \
  dJ(1)(0, 2) = b * d2_core_da_dc;                                             \
  dJ(1)(1, 0) = d_core_da + b * d2_core_da_db;                                 \
  dJ(1)(1, 1) = 2 * d_core_db + b * d2_core_db_db;                             \
  dJ(1)(1, 2) = d_core_dc + b * d2_core_db_dc;                                 \
  dJ(1)(2, 0) = b * d2_core_da_dc;                                             \
  dJ(1)(2, 1) = d_core_dc + b * d2_core_db_dc;                                 \
  dJ(1)(2, 2) = b * d2_core_dc_dc;                                             \
  dJ(2)(0, 0) = -d2_core_da_da;                                                \
  dJ(2)(0, 1) = -d2_core_da_db;                                                \
  dJ(2)(0, 2) = -d2_core_da_dc;                                                \
  dJ(2)(1, 0) = -d2_core_da_db;                                                \
  dJ(2)(1, 1) = -d2_core_db_db;                                                \
  dJ(2)(1, 2) = -d2_core_db_dc;                                                \
  dJ(2)(2, 0) = -d2_core_da_dc;                                                \
  dJ(2)(2, 1) = -d2_core_db_dc;                                                \
  dJ(2)(2, 2) = -d2_core_dc_dc;

#define CAKE_PLUS_Z_JACOBIAN                                                   \
  J(0)(0) = -(c * d_core_da);                                                  \
  J(0)(1) = -(c * d_core_db);                                                  \
  J(0)(2) = -core - c * d_core_dc;                                             \
  J(1)(0) = b * d_core_da;                                                     \
  J(1)(1) = core + b * d_core_db;                                              \
  J(1)(2) = b * d_core_dc;                                                     \
  J(2)(0) = d_core_da;                                                         \
  J(2)(1) = d_core_db;                                                         \
  J(2)(2) = d_core_dc;

#define CAKE_PLUS_Z_DJACOBIAN                                                  \
  dJ(0)(0, 0) = -(c * d2_core_da_da);                                          \
  dJ(0)(0, 1) = -(c * d2_core_da_db);                                          \
  dJ(0)(0, 2) = -d_core_da - c * d2_core_da_dc;                                \
  dJ(0)(1, 0) = -(c * d2_core_da_db);                                          \
  dJ(0)(1, 1) = -(c * d2_core_db_db);                                          \
  dJ(0)(1, 2) = -d_core_db - c * d2_core_db_dc;                                \
  dJ(0)(2, 0) = -d_core_da - c * d2_core_da_dc;                                \
  dJ(0)(2, 1) = -d_core_db - c * d2_core_db_dc;                                \
  dJ(0)(2, 2) = -2 * d_core_dc - c * d2_core_dc_dc;                            \
  dJ(1)(0, 0) = b * d2_core_da_da;                                             \
  dJ(1)(0, 1) = d_core_da + b * d2_core_da_db;                                 \
  dJ(1)(0, 2) = b * d2_core_da_dc;                                             \
  dJ(1)(1, 0) = d_core_da + b * d2_core_da_db;                                 \
  dJ(1)(1, 1) = 2 * d_core_db + b * d2_core_db_db;                             \
  dJ(1)(1, 2) = d_core_dc + b * d2_core_db_dc;                                 \
  dJ(1)(2, 0) = b * d2_core_da_dc;                                             \
  dJ(1)(2, 1) = d_core_dc + b * d2_core_db_dc;                                 \
  dJ(1)(2, 2) = b * d2_core_dc_dc;                                             \
  dJ(2)(0, 0) = d2_core_da_da;                                                 \
  dJ(2)(0, 1) = d2_core_da_db;                                                 \
  dJ(2)(0, 2) = d2_core_da_dc;                                                 \
  dJ(2)(1, 0) = d2_core_da_db;                                                 \
  dJ(2)(1, 1) = d2_core_db_db;                                                 \
  dJ(2)(1, 2) = d2_core_db_dc;                                                 \
  dJ(2)(2, 0) = d2_core_da_dc;                                                 \
  dJ(2)(2, 1) = d2_core_db_dc;                                                 \
  dJ(2)(2, 2) = d2_core_dc_dc;

#endif // MULTIPATCH_CAKE_JACOBIANS_HXX