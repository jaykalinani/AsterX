#ifndef LOOPCONTROL_H
#define LOOPCONTROL_H

#ifdef __cplusplus
extern "C" {
#endif

// this must minic the GridDescBase in loop.hxx

#define LC_DIM 3

typedef struct {
  int gsh[LC_DIM];
  int lbnd[LC_DIM], ubnd[LC_DIM];
  int lsh[LC_DIM];
  int ash[LC_DIM];
  int bbox[2 * LC_DIM];
  int nghostzones[LC_DIM];
  int tmin[LC_DIM], tmax[LC_DIM];
  CCTK_REAL x0[LC_DIM];
  CCTK_REAL dx[LC_DIM];
} GridDescBase_t;

GridDescBase_t LC_CreateGridDesc(const cGH *cctkGH);

// replicate looping constructs in loop.hxx as good as we can, given that there
// is only one CCTK_LOOP macro thta cannot distinguish between cell and vertex
// variables
//
// this is not quite the way LoopControl does it which would involve changing
// LC_LOOP3STROFF_NORMAL which no longer knows what it needs to loop over and
// is only given the imin/imax ranges that Cactus thinks are required but has
// no idea about tiles

#define LC_LOOP3_BOX(name, i,j,k, grid, imin,imax, inormal) \
for (int d = 0; d < dim; ++d)                               \
  if (imin[d] >= imax[d])                                   \
    return;                                                 \
                                                            \
for (int k = imin[2]; k < imax[2]; ++k) {                   \
  for (int j = imin[1]; j < imax[1]; ++j) {                 \
    for (int i = imin[0]; i < imax[0]; ++i) {

#define LC_ENDLOOP3_BOX(name)                               \
    }                                                       \
  }                                                         \
}

/* Replace CCTK_LOOP macros */
#if !defined CCTK_LOOP3STROFF_NORMAL || !defined CCTK_ENDLOOP3STROFF_NORMAL
#error "internal error"
#endif
#undef CCTK_LOOP3STROFF_NORMAL
#undef CCTK_ENDLOOP3STROFF_NORMAL

#undef CCTK_LOOP3_ALL
#undef CCTK_ENDLOOP3_ALL

#define CCTK_LOOP3_ALL(name, cctki3_cctkGH_, i,j,k)                             \
do {                                                                            \
  GridDescBase_t grid = LC_CreateGridDesc(cctki3_cctkGH_);                      \
  int imin[LC_DIM], imax[LC_DIM], inormal[LC_DIM];                              \
  for (int d = 0; d < dim; ++d) {                                               \
    imin[d] = grid.tmin[d];                                                     \
    imax[d] = grid.tmax[d] + (grid.tmax[d] >= grid.lsh[d] ? grid.offset[d] : 0);\
    if(imax[d] > grid.lsh[d] + grid.offset[d]) {                                \
      imax[d] = grid.lsh[d] + grid.offset[d];                                   \
    }                                                                           \
  }                                                                             \
  int inormal[3] = {0, 0, 0};                                                   \
  LC_LOOP3_BOX(name, i,j,k, grid, imin,imax, inormal)

#define CCTK_ENDLOOP3_ALL(name)                                                 \
  LC_ENDLOOP3_BOX(name)                                                         \
}

#ifdef __cplusplus
}
#endif

#endif // LOOPCONTROL_H
