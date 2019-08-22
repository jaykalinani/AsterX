#include "prolongate_3d_rf2.hxx"

#include "cctk.h"

#include "typeprops.hh"

#include "vectors.h"

#include <AMReX_Vector.H>

#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>

namespace CarpetLib {
using namespace std;

static inline size_t index3(size_t const i, size_t const j, size_t const k,
                            size_t const padexti, size_t const padextj,
                            size_t const padextk, size_t const exti,
                            size_t const extj, size_t const extk) {
#ifdef CARPET_DEBUG
  assert(i < exti);
  assert(j < extj);
  assert(k < extk);
#endif

  return i + padexti * (j + padextj * k);
}

template <typename T>
inline T ipow_helper(T x, unsigned int y) CCTK_ATTRIBUTE_CONST;
template <typename T> inline T ipow_helper(T x, unsigned int y) {                                         T z = y & 1 ? x : T(1);
  while (y >>= 1) {
    x *= x;
    if (y & 1)                                                                                                z *= x;
  }
  return z;
}                                                                                                       
template <class T> T ipow(T const x, int const y) {
  if (y < 0)
    return T(1) / ipow_helper(x, -y);
  else
    return ipow_helper(x, y);
}

#define SRCIND3(i, j, k)                                                       \
  index3(i, j, k, srcipadext, srcjpadext, srckpadext, srciext, srcjext, srckext)
#define DSTIND3(i, j, k)                                                       \
  index3(i, j, k, dstipadext, dstjpadext, dstkpadext, dstiext, dstjext, dstkext)

namespace coeffs_3d_rf2 {

// Names for types

#ifdef HAVE_CCTK_INT1
inline const char *typestring(const CCTK_INT1 &) { return "CCTK_INT1"; }
#endif

#ifdef HAVE_CCTK_INT2
inline const char *typestring(const CCTK_INT2 &) { return "CCTK_INT2"; }
#endif

#ifdef HAVE_CCTK_INT4
inline const char *typestring(const CCTK_INT4 &) { return "CCTK_INT4"; }
#endif

#ifdef HAVE_CCTK_INT8
inline const char *typestring(const CCTK_INT8 &) { return "CCTK_INT8"; }
#endif

#ifdef HAVE_CCTK_INT16
inline const char *typestring(const CCTK_INT16 &) { return "CCTK_INT16"; }
#endif

#ifdef HAVE_CCTK_REAL4
inline const char *typestring(const CCTK_REAL4 &) { return "CCTK_REAL4"; }
#endif

#ifdef HAVE_CCTK_REAL8
inline const char *typestring(const CCTK_REAL8 &) { return "CCTK_REAL8"; }
#endif

#ifdef HAVE_CCTK_REAL16
inline const char *typestring(const CCTK_REAL16 &) { return "CCTK_REAL16"; }
#endif

#ifdef HAVE_CCTK_REAL4
inline const char *typestring(const CCTK_COMPLEX8 &) { return "CCTK_COMPLEX8"; }
#endif

#ifdef HAVE_CCTK_REAL8
inline const char *typestring(const CCTK_COMPLEX16 &) {
  return "CCTK_COMPLEX16";
}
#endif

#ifdef HAVE_CCTK_REAL16
inline const char *typestring(const CCTK_COMPLEX32 &) {
  return "CCTK_COMPLEX32";
}
#endif

// 1D interpolation coefficients

template <typename RT, int ORDER> struct coeffs1d {
  static const RT
      coeffs[];

  static ptrdiff_t const ncoeffs = ORDER + 1;
  static ptrdiff_t const imin = -ncoeffs / 2 + 1;
  static ptrdiff_t const imax = imin + ncoeffs;

  static inline RT const &get(ptrdiff_t const i) {
    static_assert(ncoeffs == sizeof coeffs / sizeof *coeffs,
                  "coefficient array has wrong size");
#ifdef CARPET_DEBUG
    assert(i >= imin and i < imax);
#endif
    ptrdiff_t const j = i - imin;
#ifdef CARPET_DEBUG
    assert(j >= 0 and j < ncoeffs);
#endif
    return coeffs[j];
  }

  // Test coefficients
  static void test() {
    static bool tested = false;
    if (tested)
      return;
    tested = true;

    static_assert(ncoeffs == sizeof coeffs / sizeof *coeffs,
                  "coefficient array has wrong size");

    // Do not test integer operators (they should be disabled anyway)
    if (std::fabs(RT(0.5) - 0.5) > 1.0e-5)
      return;

    // Test all orders
    bool error = false;
    for (int n = 0; n <= ORDER; ++n) {
      RT res = RT(0.0);
      for (ptrdiff_t i = imin; i < imax; ++i) {
        RT const x = RT(i) - RT(0.5);
        RT const y = ipow(x, n);
        res += get(i) * y;
      }
      RT const x0 = RT(0.0);
      RT const y0 = ipow(x0, n);
      // Allow losing 3 digits:
      CCTK_REAL const eps = RT(1.0e+3) * numeric_limits<RT>::epsilon();
      if (not(std::fabs(res - y0) < eps)) {
        RT rt;
        ostringstream buf;
        buf << "Error in prolongate_3d_rf2::coeffs_3d_rf2\n"
            << "   RT=" << typestring(rt) << "\n"
            << "   ORDER=" << ORDER << "\n"
            << "   n=" << n << "\n"
            << "   y0=" << y0 << ", res=" << res;
        CCTK_WARN(CCTK_WARN_ALERT, buf.str().c_str());
        error = true;
      }
    } // for n
    if (error)
      CCTK_ERROR("Aborting.");
  }
};

#define TYPECASE(N, RT)                                                        \
                                                                               \
  template <>                                                                  \
  const RT coeffs1d<RT, 1>::coeffs[] = {+1 / RT(2.0), +1 / RT(2.0)};           \
                                                                               \
  template <>                                                                  \
  const RT coeffs1d<RT, 3>::coeffs[] = {-1 / RT(16.0), +9 / RT(16.0),          \
                                        +9 / RT(16.0), -1 / RT(16.0)};         \
                                                                               \
  template <>                                                                  \
  const RT coeffs1d<RT, 5>::coeffs[] = {+3 / RT(256.0),   -25 / RT(256.0),     \
                                        +150 / RT(256.0), +150 / RT(256.0),    \
                                        -25 / RT(256.0),  +3 / RT(256.0)};     \
                                                                               \
  template <>                                                                  \
  const RT coeffs1d<RT, 7>::coeffs[] = {                                       \
      -5 / RT(2048.0),    +49 / RT(2048.0),   -245 / RT(2048.0),               \
      +1225 / RT(2048.0), +1225 / RT(2048.0), -245 / RT(2048.0),               \
      +49 / RT(2048.0),   -5 / RT(2048.0)};                                    \
                                                                               \
  template <>                                                                  \
  const RT coeffs1d<RT, 9>::coeffs[] = {                                       \
      +35 / RT(65536.0),   -405 / RT(65536.0),   +2268 / RT(65536.0),          \
      -8820 / RT(65536.0), +39690 / RT(65536.0), +39690 / RT(65536.0),         \
      -8820 / RT(65536.0), +2268 / RT(65536.0),  -405 / RT(65536.0),           \
      +35 / RT(65536.0)};                                                      \
                                                                               \
  template <>                                                                  \
  const RT coeffs1d<RT, 11>::coeffs[] = {                                      \
      -63 / RT(524288.0),     +847 / RT(524288.0),   -5445 / RT(524288.0),     \
      +22869 / RT(524288.0),  -76230 / RT(524288.0), +320166 / RT(524288.0),   \
      +320166 / RT(524288.0), -76230 / RT(524288.0), +22869 / RT(524288.0),    \
      -5445 / RT(524288.0),   +847 / RT(524288.0),   -63 / RT(524288.0)};

TYPECASE(CCTK_VARIABLE_REAL8, CCTK_REAL8);
#undef TYPECASE

} // namespace coeffs_3d_rf2

using namespace coeffs_3d_rf2;

// 0D "interpolation"
template <typename T, int ORDER>
static inline T const &interp0(T const *restrict const p) {
  return *p;
}

// 1D interpolation
template <typename T, int ORDER, int di>
static inline T interp1(T const *restrict const p, size_t const d1) {
  typedef typeprops<T> typ;
  typedef typename typ::real RT;
  typedef coeffs1d<RT, ORDER> coeffs;
  if (di == 0) {
    return interp0<T, ORDER>(p);
  } else {
    if (d1 == 1) {
#if 0
        T res = typ::fromreal (0);
        for (ptrdiff_t i=coeffs::imin; i<coeffs::imax; ++i) {
          res += coeffs::get(i) * interp0<T,ORDER> (p + i);
        }
        return res;
#endif
      typedef vecprops<T> VP;
      typedef typename VP::vector_t VT;
      ptrdiff_t const vsize = VP::size();
      ptrdiff_t i = coeffs::imin;
      T res = typ::fromreal(0);
      if (coeffs::ncoeffs >= vsize) {
        VT vres = VP::mul(VP::load(typ::fromreal(coeffs::get(i))),
                          VP::loadu(interp0<T, ORDER>(p + i)));
        i += vsize;
#if defined(__INTEL_COMPILER)
        // Unroll the loop manually to help the Intel compiler
        // (This manual unrolling hurts with other compilers, e.g. PGI)
        assert(coeffs::ncoeffs / vsize <= 12);
        switch (coeffs::ncoeffs / vsize) {
        // Note that all case statements fall through
        case 12:
          vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                          VP::loadu(interp0<T, ORDER>(p + i)), vres);
          i += vsize;
        case 11:
          vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                          VP::loadu(interp0<T, ORDER>(p + i)), vres);
          i += vsize;
        case 10:
          vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                          VP::loadu(interp0<T, ORDER>(p + i)), vres);
          i += vsize;
        case 9:
          vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                          VP::loadu(interp0<T, ORDER>(p + i)), vres);
          i += vsize;
        case 8:
          vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                          VP::loadu(interp0<T, ORDER>(p + i)), vres);
          i += vsize;
        case 7:
          vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                          VP::loadu(interp0<T, ORDER>(p + i)), vres);
          i += vsize;
        case 6:
          vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                          VP::loadu(interp0<T, ORDER>(p + i)), vres);
          i += vsize;
        case 5:
          vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                          VP::loadu(interp0<T, ORDER>(p + i)), vres);
          i += vsize;
        case 4:
          vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                          VP::loadu(interp0<T, ORDER>(p + i)), vres);
          i += vsize;
        case 3:
          vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                          VP::loadu(interp0<T, ORDER>(p + i)), vres);
          i += vsize;
        case 2:
          vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                          VP::loadu(interp0<T, ORDER>(p + i)), vres);
          i += vsize;
        }
#else
        for (; i + vsize <= ptrdiff_t(coeffs::imax); i += vsize) {
          vres = VP::madd(VP::load(typ::fromreal(coeffs::get(i))),
                          VP::loadu(interp0<T, ORDER>(p + i)), vres);
        }
#endif
        for (int d = 0; d < vsize; ++d) {
          res += VP::elt(vres, d);
        }
      }
      assert(i == ptrdiff_t(coeffs::imax) - ptrdiff_t(coeffs::ncoeffs % vsize));
      for (i = coeffs::imax - ptrdiff_t(coeffs::ncoeffs % vsize);
           i < coeffs::imax; ++i) {
        res += coeffs::get(i) * interp0<T, ORDER>(p + i * d1);
      }
      return res;
    } else {
      assert(0); // why would d1 have a non-unit stride?
      T res = typ::fromreal(0);
      for (ptrdiff_t i = coeffs::imin; i < coeffs::imax; ++i) {
        res += coeffs::get(i) * interp0<T, ORDER>(p + i * d1);
      }
      return res;
    }
  }
}

// 2D interpolation
template <typename T, int ORDER, int di, int dj>
static inline T interp2(T const *restrict const p, size_t const d1,
                        size_t const d2) {
  typedef typeprops<T> typ;
  typedef typename typ::real RT;
  typedef coeffs1d<RT, ORDER> coeffs;
  if (dj == 0) {
    return interp1<T, ORDER, di>(p, d1);
  } else {
    T res = typeprops<T>::fromreal(0);
    for (ptrdiff_t i = coeffs::imin; i < coeffs::imax; ++i) {
      res += coeffs::get(i) * interp1<T, ORDER, di>(p + i * d2, d1);
    }
    return res;
  }
}

// 3D interpolation
template <typename T, int ORDER, int di, int dj, int dk>
static inline T interp3(T const *restrict const p, size_t const d1,
                        size_t const d2, size_t const d3) {
  typedef typeprops<T> typ;
  typedef typename typ::real RT;
  typedef coeffs1d<RT, ORDER> coeffs;
  if (dk == 0) {
    return interp2<T, ORDER, di, dj>(p, d1, d2);
  } else {
    T res = typeprops<T>::fromreal(0);
    for (ptrdiff_t i = coeffs::imin; i < coeffs::imax; ++i) {
      res += coeffs::get(i) * interp2<T, ORDER, di, dj>(p + i * d3, d1, d2);
    }
    return res;
  }
}

// Check interpolation index ranges
template <typename T, int ORDER> static inline void check_indices0() {}

template <typename T, int ORDER, int di>
static inline void check_indices1(ptrdiff_t const is, ptrdiff_t const srciext) {
  static_assert(di == 0 or di == 1, "di must be 0 or 1");
#ifdef CARPET_DEBUG
  typedef typename typeprops<T>::real RT;
  typedef coeffs1d<RT, ORDER> coeffs;
  if (di == 0) {
    assert(is >= 0);
    assert(is < srciext);
  } else {
    assert(is + coeffs::imin >= 0);
    assert(is + coeffs::imax <= srciext);
  }
  check_indices0<T, ORDER>();
#endif
}

template <typename T, int ORDER, int di, int dj>
static inline void check_indices2(ptrdiff_t const is, ptrdiff_t const js,
                                  ptrdiff_t const srciext,
                                  ptrdiff_t const srcjext) {
  static_assert(dj == 0 or dj == 1, "dj must be 0 or 1");
#ifdef CARPET_DEBUG
  typedef typename typeprops<T>::real RT;
  typedef coeffs1d<RT, ORDER> coeffs;
  if (dj == 0) {
    assert(js >= 0);
    assert(js < srcjext);
  } else {
    assert(js + coeffs::imin >= 0);
    assert(js + coeffs::imax <= srcjext);
  }
  check_indices1<T, ORDER, di>(is, srciext);
#endif
}

template <typename T, int ORDER, int di, int dj, int dk>
static inline void check_indices3(ptrdiff_t const is, ptrdiff_t const js,
                                  ptrdiff_t const ks, ptrdiff_t const srciext,
                                  ptrdiff_t const srcjext,
                                  ptrdiff_t const srckext) {
  static_assert(dk == 0 or dk == 1, "dk must be 0 or 1");
#ifdef CARPET_DEBUG
  typedef typename typeprops<T>::real RT;
  typedef coeffs1d<RT, ORDER> coeffs;
  if (dk == 0) {
    assert(ks >= 0);
    assert(ks < srckext);
  } else {
    assert(ks + coeffs::imin >= 0);
    assert(ks + coeffs::imax <= srckext);
  }
  check_indices2<T, ORDER, di, dj>(is, js, srciext, srcjext);
#endif
}

template <typename T, int ORDER>
void prolongate_3d_rf2(T *restrict const dst,
                       int const *restrict flo, int const *restrict fhi, //dstpadext,
                       int const *restrict fblo, int const *restrict fbhi, //dstext,
                       int const ncomp, int const ratioX, int const ratioY, int const ratioZ,
                       T const *restrict const src,
                       int const *restrict clo, int const *restrict chi//, //srcpadext,
                       //int const *restrict cblo, int const *restrict cbhi, //srcext,
                       /*ibbox3 const &restrict srcbbox,
                       ibbox3 const &restrict dstbbox, ibbox3 const &restrict,
                       ibbox3 const &restrict regbbox, void *extraargs*/) {

  static_assert(ORDER >= 0 and ORDER % 2 == 1,
                "ORDER must be non-negative and odd");

  typedef typename typeprops<T>::real RT;
  coeffs1d<RT, ORDER>::test();

  //if (any(srcbbox.stride() <= regbbox.stride() or
  //        dstbbox.stride() != regbbox.stride()))
  //  CCTK_ERROR("Internal error: strides disagree");

  //if (any(srcbbox.stride() != reffact2 * dstbbox.stride()))
  //  CCTK_ERROR(
  //      "Internal error: source strides are not twice the destination strides");

  //if (any(srcbbox.lower() % srcbbox.stride() != 0))
  //  CCTK_ERROR("Internal error: source bbox is not aligned with vertices");
  //if (any(dstbbox.lower() % dstbbox.stride() != 0))
  //  CCTK_ERROR("Internal error: destination bbox is not aligned with vertices");
  //if (any(regbbox.lower() % regbbox.stride() != 0))
  //  CCTK_ERROR("Internal error: prolongation region bbox is not aligned with "
  //             "vertices");

  // This could be handled, but is likely to point to an error elsewhere
  //if (regbbox.empty())
  //  CCTK_ERROR("Internal error: region extent is empty");

  int const regext[3] = {fbhi[0] - fblo[0] + 1, fbhi[1] - fblo[1] + 1, fbhi[2] - fblo[2] + 1};
  //assert(all((regbbox.lower() - srcbbox.lower()) % regbbox.stride() == 0));
  int const srcoff[3] = {0,0,0};
  //assert(all((regbbox.lower() - dstbbox.lower()) % regbbox.stride() == 0));
  int const dstoff[3] = {0,0,0};

  //bvect3 const needoffsetlo = srcoff % reffact2 != 0;
  //bvect3 const needoffsethi = (srcoff + regext - 1) % reffact2 != 0;
  //int const offsetlo[3] = {0,0,0};
  //int const offsethi[3] = {0,0,0};

  //if (not regbbox.expand(offsetlo, offsethi).is_contained_in(srcbbox) or
  //    not regbbox.is_contained_in(dstbbox)) {
  //  cerr << "ORDER=" << ORDER << "\n"
  //       << "offsetlo=" << offsetlo << "\n"
  //       << "offsethi=" << offsethi << "\n"
  //       << "regbbox=" << regbbox << "\n"
  //       << "dstbbox=" << dstbbox << "\n"
  //       << "regbbox.expand=" << regbbox.expand(offsetlo, offsethi) << "\n"
  //       << "srcbbox=" << srcbbox << "\n";
  //  CCTK_ERROR(
  //      "Internal error: region extent is not contained in array extent");
  //}

  size_t const srcipadext = chi[0] - clo[0] + 1; //srcpadext[0];
  size_t const srcjpadext = chi[1] - clo[1] + 1; //srcpadext[1];
  size_t const srckpadext = chi[2] - clo[2] + 1; //srcpadext[2];

  size_t const dstipadext = fhi[0] - flo[0] + 1; // dstpadext[0];
  size_t const dstjpadext = fhi[1] - flo[1] + 1; // dstpadext[1];
  size_t const dstkpadext = fhi[2] - flo[2] + 1; // dstpadext[2];

  size_t const srciext = chi[0] - clo[0] + 1; // rcext[0];
  size_t const srcjext = chi[1] - clo[1] + 1; // rcext[1];
  size_t const srckext = chi[2] - clo[2] + 1; // rcext[2];

  size_t const dstiext = fbhi[0] - fblo[0] + 1; // stext[0];
  size_t const dstjext = fbhi[1] - fblo[1] + 1; // stext[1];
  size_t const dstkext = fbhi[2] - fblo[2] + 1; // stext[2];

  size_t const regiext = regext[0];
  size_t const regjext = regext[1];
  size_t const regkext = regext[2];

  size_t const srcioff = srcoff[0];
  size_t const srcjoff = srcoff[1];
  size_t const srckoff = srcoff[2];

  size_t const dstioff = dstoff[0];
  size_t const dstjoff = dstoff[1];
  size_t const dstkoff = dstoff[2];

  size_t const fi = srcioff % 2;
  size_t const fj = srcjoff % 2;
  size_t const fk = srckoff % 2;

  size_t const i0 = srcioff / 2;
  size_t const j0 = srcjoff / 2;
  size_t const k0 = srckoff / 2;

  // size_t const srcdi = SRCIND3(1,0,0) - SRCIND3(0,0,0);
  size_t const srcdi = 1;
  assert(srcdi == (srciext > 1 ? SRCIND3(1, 0, 0) - SRCIND3(0, 0, 0) : 1));
  size_t const srcdj = srcjext > 1 ? SRCIND3(0, 1, 0) - SRCIND3(0, 0, 0) : 0;
  size_t const srcdk = srckext > 1 ? SRCIND3(0, 0, 1) - SRCIND3(0, 0, 0) : 0;

  // Loop over fine region
  // Label scheme: l 8 fk fj fi

  size_t is, js, ks;
  size_t id, jd, kd;
  size_t i, j, k;

  // begin k loop
  k = 0;
  ks = k0;
  kd = dstkoff;
  if (fk == 0)
    goto l80;
  goto l81;

// begin j loop
l80:
  j = 0;
  js = j0;
  jd = dstjoff;
  if (fj == 0)
    goto l800;
  goto l801;

// begin i loop
l800:
  i = 0;
  is = i0;
  id = dstioff;
  if (fi == 0)
    goto l8000;
  goto l8001;

// kernel
l8000:
  check_indices3<T, ORDER, 0, 0, 0>(is, js, ks, srciext, srcjext, srckext);
  dst[DSTIND3(id, jd, kd)] = interp3<T, ORDER, 0, 0, 0>(
      &src[SRCIND3(is, js, ks)], srcdi, srcdj, srcdk);
  i = i + 1;
  id = id + 1;
  if (i < regiext)
    goto l8001;
  goto l900;

// kernel
l8001:
  check_indices3<T, ORDER, 1, 0, 0>(is, js, ks, srciext, srcjext, srckext);
  dst[DSTIND3(id, jd, kd)] = interp3<T, ORDER, 1, 0, 0>(
      &src[SRCIND3(is, js, ks)], srcdi, srcdj, srcdk);
  i = i + 1;
  id = id + 1;
  is = is + 1;
  if (i < regiext)
    goto l8000;
  goto l900;

// end i loop
l900:
  j = j + 1;
  jd = jd + 1;
  if (j < regjext)
    goto l801;
  goto l90;

// begin i loop
l801:
  i = 0;
  is = i0;
  id = dstioff;
  if (fi == 0)
    goto l8010;
  goto l8011;

// kernel
l8010:
  check_indices3<T, ORDER, 0, 1, 0>(is, js, ks, srciext, srcjext, srckext);
  dst[DSTIND3(id, jd, kd)] = interp3<T, ORDER, 0, 1, 0>(
      &src[SRCIND3(is, js, ks)], srcdi, srcdj, srcdk);
  i = i + 1;
  id = id + 1;
  if (i < regiext)
    goto l8011;
  goto l901;

// kernel
l8011:
  check_indices3<T, ORDER, 1, 1, 0>(is, js, ks, srciext, srcjext, srckext);
  dst[DSTIND3(id, jd, kd)] = interp3<T, ORDER, 1, 1, 0>(
      &src[SRCIND3(is, js, ks)], srcdi, srcdj, srcdk);
  i = i + 1;
  id = id + 1;
  is = is + 1;
  if (i < regiext)
    goto l8010;
  goto l901;

// end i loop
l901:
  j = j + 1;
  jd = jd + 1;
  js = js + 1;
  if (j < regjext)
    goto l800;
  goto l90;

// end j loop
l90:
  k = k + 1;
  kd = kd + 1;
  if (k < regkext)
    goto l81;
  goto l9;

// begin j loop
l81:
  j = 0;
  js = j0;
  jd = dstjoff;
  if (fj == 0)
    goto l810;
  goto l811;

// begin i loop
l810:
  i = 0;
  is = i0;
  id = dstioff;
  if (fi == 0)
    goto l8100;
  goto l8101;

// kernel
l8100:
  check_indices3<T, ORDER, 0, 0, 1>(is, js, ks, srciext, srcjext, srckext);
  dst[DSTIND3(id, jd, kd)] = interp3<T, ORDER, 0, 0, 1>(
      &src[SRCIND3(is, js, ks)], srcdi, srcdj, srcdk);
  i = i + 1;
  id = id + 1;
  if (i < regiext)
    goto l8101;
  goto l910;

// kernel
l8101:
  check_indices3<T, ORDER, 1, 0, 1>(is, js, ks, srciext, srcjext, srckext);
  dst[DSTIND3(id, jd, kd)] = interp3<T, ORDER, 1, 0, 1>(
      &src[SRCIND3(is, js, ks)], srcdi, srcdj, srcdk);
  i = i + 1;
  id = id + 1;
  is = is + 1;
  if (i < regiext)
    goto l8100;
  goto l910;

// end i loop
l910:
  j = j + 1;
  jd = jd + 1;
  if (j < regjext)
    goto l811;
  goto l91;

// begin i loop
l811:
  i = 0;
  is = i0;
  id = dstioff;
  if (fi == 0)
    goto l8110;
  goto l8111;

// kernel
l8110:
  check_indices3<T, ORDER, 0, 1, 1>(is, js, ks, srciext, srcjext, srckext);
  dst[DSTIND3(id, jd, kd)] = interp3<T, ORDER, 0, 1, 1>(
      &src[SRCIND3(is, js, ks)], srcdi, srcdj, srcdk);
  i = i + 1;
  id = id + 1;
  if (i < regiext)
    goto l8111;
  goto l911;

// kernel
l8111:
  check_indices3<T, ORDER, 1, 1, 1>(is, js, ks, srciext, srcjext, srckext);
  dst[DSTIND3(id, jd, kd)] = interp3<T, ORDER, 1, 1, 1>(
      &src[SRCIND3(is, js, ks)], srcdi, srcdj, srcdk);
  i = i + 1;
  id = id + 1;
  is = is + 1;
  if (i < regiext)
    goto l8110;
  goto l911;

// end i loop
l911:
  j = j + 1;
  jd = jd + 1;
  js = js + 1;
  if (j < regjext)
    goto l810;
  goto l91;

// end j loop
l91:
  k = k + 1;
  kd = kd + 1;
  ks = ks + 1;
  if (k < regkext)
    goto l80;
  goto l9;

// end k loop
l9:;
}

} // namespace CarpetLib

namespace amrex {
template <int ORDER>
Prolongate3D_VC_RF2<ORDER>::~Prolongate3D_VC_RF2 () { }

template <int ORDER>
Box
Prolongate3D_VC_RF2<ORDER>::CoarseBox (const Box& fine, int ratio)
{
  assert(ORDER % 2);

  Box b = amrex::coarsen(fine, ratio);
  b.grow(ORDER/2);

  return b;
}

template <int ORDER>
Box
Prolongate3D_VC_RF2<ORDER>::CoarseBox (const Box& fine, const IntVect& ratio)
{
  assert(ORDER % 2);

  Box b = amrex::coarsen(fine, ratio);
  b.grow(ORDER/2);

  return b;
}

template <int ORDER>
void
Prolongate3D_VC_RF2<ORDER>::interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec> const& bcr,
                         int              actual_comp,
                         int              actual_state,
                         RunOn            gpu_or_cpu)
{
    assert(bcr.size() >= ncomp);
    //
    // Make box which is intersection of fine_region and domain of fine.
    //
    Box target_fine_region = fine_region & fine.box();

    //Box crse_bx(amrex::coarsen(target_fine_region,ratio));
    
    //
    // Get coarse and fine edge-centered volume coordinates.
    //
    // TODO: do not use AMREX macro
    /*
    Vector<CCTK_REAL> fvc[AMREX_SPACEDIM];
    Vector<CCTK_REAL> cvc[AMREX_SPACEDIM];
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        fine_geom.GetEdgeVolCoord(fvc[dir],target_fine_region,dir);
        crse_geom.GetEdgeVolCoord(cvc[dir],crse_bx,dir);
    }
    */
    CCTK_REAL* fdat        = fine.dataPtr(fine_comp);
    const CCTK_REAL* cdat  = crse.dataPtr(crse_comp);
    const int* flo    = fine.loVect();
    const int* fhi    = fine.hiVect();
    const int* fblo   = target_fine_region.loVect();
    const int* fbhi   = target_fine_region.hiVect();
    const int* clo    = crse.loVect();
    const int* chi    = crse.hiVect();
    //const int* cblo   = crse_bx.loVect();
    //const int* cbhi   = crse_bx.hiVect();
    Vector<int> bc     = GetBCArray(bcr);
    const int* ratioV = ratio.getVect();

    CarpetLib::prolongate_3d_rf2<CCTK_REAL, ORDER>(fdat,flo,fhi,
                   fblo, fbhi,
                   ncomp,ratioV[0],ratioV[1],ratioV[2],
                   cdat,clo,chi//,
                   //cblo, cbhi,
                   /*fslo,fshi,
                   bc.dataPtr(), 
                   fvc[0].dataPtr(),fvc[1].dataPtr(),fvc[2].dataPtr(),
                   cvc[0].dataPtr(),cvc[1].dataPtr(),cvc[2].dataPtr(),
                   actual_comp,actual_state*/);
}

Prolongate3D_VC_RF2<5> prolongate3d_vc_rf2_o5;
} // namespace amrex

