#ifndef LOOP_HXX
#define LOOP_HXX

#include <cctk.h>

#include <array>
#include <iostream>
#include <string>

namespace Loop {
using namespace std;

constexpr int dim = 3;

struct GridDescBase {
  array<int, dim> gsh;
  array<int, dim> lbnd, ubnd;
  array<int, dim> lsh;
  array<int, dim> ash;
  array<int, 2 * dim> bbox;
  array<int, dim> nghostzones;
  array<int, dim> tmin, tmax;

  template <typename T, size_t N>
  static void output(ostream &os, const string &str, const array<T, N> &arr) {
    os << str << ":[";
    for (size_t n = 0; n < N; ++n) {
      if (n > 0)
        os << ",";
      os << arr[n];
    }
    os << "]";
  }
  friend ostream &operator<<(ostream &os, const GridDescBase &grid) {
    os << "GridDescBase{";
    output(os, "gsh", grid.gsh);
    output(os, ",lbnd", grid.lbnd);
    output(os, ",ubnd", grid.ubnd);
    output(os, ",lsh", grid.lsh);
    output(os, ",bbox", grid.bbox);
    output(os, ",nghostzones", grid.nghostzones);
    output(os, ",tmin", grid.tmin);
    output(os, ",tmax", grid.tmax);
    os << "}";
    return os;
  }

protected:
  GridDescBase();

public:
  GridDescBase(const cGH *cctkGH);

  // Loop over a given box
  template <typename F>
  void loop_box(const F &f, const array<int, dim> &restrict imin,
                const array<int, dim> &restrict imax) const {
    // cout << *this;
    // output(cout, ",imin", imin);
    // output(cout, ",imax", imax);
    // cout << "\n";

    for (int d = 0; d < dim; ++d)
      if (imin[d] >= imax[d])
        return;

    constexpr int di = 1;
    const int dj = di * ash[0];
    const int dk = dj * ash[1];

    for (int k = imin[2]; k < imax[2]; ++k) {
      for (int j = imin[1]; j < imax[1]; ++j) {
#pragma omp simd
        for (int i = imin[0]; i < imax[0]; ++i) {
          int idx = i * di + j * dj + k * dk;
          f(i, j, k, idx);
        }
      }
    }
  }

  // Loop over all points
  template <typename F> void loop_all(const F &f) const {
    array<int, dim> imin, imax;
    for (int d = 0; d < dim; ++d) {
      imin[d] = max(tmin[d], 0);
      imax[d] = min(tmax[d], lsh[d]);
    }

    loop_box(f, imin, imax);
  }

  // Loop over all interior points
  template <typename F> void loop_int(const F &f) const {
    array<int, dim> imin, imax;
    for (int d = 0; d < dim; ++d) {
      imin[d] = max(tmin[d], nghostzones[d]);
      imax[d] = min(tmax[d], lsh[d] - nghostzones[d]);
    }

    loop_box(f, imin, imax);
  }

  // Loop over all outer boundary points. This excludes ghost faces, but
  // includes ghost edges/corners on non-ghost faces.
  template <typename F> void loop_bnd(const F &f) const {
    for (int dir = 0; dir < dim; ++dir) {
      for (int face = 0; face < 2; ++face) {
        if (bbox[2 * dir + face]) {

          array<int, dim> imin, imax;
          for (int d = 0; d < dim; ++d) {
            // by default, include interior and outer boundaries and ghosts
            imin[d] = 0;
            imax[d] = lsh[d];

            // avoid covering edges and corners multiple times
            if (d < dir) {
              if (bbox[2 * d])
                imin[d] = nghostzones[d]; // only interior
              if (bbox[2 * d + 1])
                imax[d] = lsh[d] - nghostzones[d]; // only interior
            }
          }
          // only one face on outer boundary
          if (face == 0)
            imax[dir] = nghostzones[dir];
          else
            imin[dir] = lsh[dir] - nghostzones[dir];

          for (int d = 0; d < dim; ++d) {
            imin[d] = max(tmin[d], imin[d]);
            imax[d] = min(tmax[d], imax[d]);
          }

          loop_box(f, imin, imax);
        }
      }
    }
  }
};

template <typename F> void loop_all(const cGH *cctkGH, const F &f) {
  GridDescBase(cctkGH).loop_all(f);
}

template <typename F> void loop_int(const cGH *cctkGH, const F &f) {
  GridDescBase(cctkGH).loop_int(f);
}

template <typename F> void loop_bnd(const cGH *cctkGH, const F &f) {
  GridDescBase(cctkGH).loop_bnd(f);
}

} // namespace Loop

#endif // #ifndef LOOP_HXX
