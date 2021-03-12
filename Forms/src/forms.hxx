#include <cxxabi.h>

#include <array>
#include <cassert>
#include <cstddef>
#include <initializer_list>
#include <memory>
#include <ostream>
#include <string>
#include <type_traits>

namespace Forms {
using namespace std;

constexpr int bitsign(bool s) { return s ? -1 : +1; }
constexpr int bitsign(int s) { return bitsign(bool(s & 1)); }

constexpr size_t count_bits(size_t x) {
  size_t n = 0;
  while (x > 0) {
    n += x & 1;
    x >>= 1;
  }
  return n;
}

// Taken from
// <https://stackoverflow.com/questions/9640256/define-multiple-methods-with-parameters-from-variadic-templates>
inline string demangle(const char *mangled) {
  int status;
  unique_ptr<char[], void (*)(void *)> result(
      abi::__cxa_demangle(mangled, 0, 0, &status), std::free);
  return result.get() ? string(result.get()) : "<unknown type>";
}

template <typename T, size_t D, size_t R> struct form;

template <typename T, size_t D, size_t R> struct form {
  static_assert(D >= 0, "");
  static_assert(R >= 0 && R <= D, "");

  constexpr static size_t ncomponents() {
    size_t ncomps = 0;
    for (size_t idx = 0; idx < (1 << D); ++idx)
      ncomps += count_bits(idx) == R;
    return ncomps;
  }

  constexpr static size_t idx2comp(const size_t idx0) {
    assert(count_bits(idx0) == R);
    size_t comp = 0;
    for (size_t idx = 0; idx < (1 << D); ++idx) {
      if (idx == idx0)
        break;
      comp += count_bits(idx) == R;
    }
    return comp;
  }

  constexpr static size_t comp2idx(const size_t comp0) {
    size_t comp = 0;
    size_t idx;
    for (idx = 0; idx < (1 << D); ++idx) {
      if (count_bits(idx) == R) {
        if (comp == comp0)
          break;
        ++comp;
      }
    }
    assert(count_bits(idx) == R);
    return idx;
  }

  static constexpr size_t N = ncomponents();
  array<T, N> comps;

  form() : comps() {}
  form(const array<T, N> &comps) : comps(comps) {}

  static constexpr form pure(const T &rc) {
    form r;
    for (size_t n = 0; n < N; ++n)
      r.comps[n] = rc;
    return r;
  }

  template <typename T1, size_t D1, size_t R1>
  friend constexpr form<T1, D1, R1> make_pure(const T1 &rc);

  template <typename T1, size_t D1, size_t R1, typename F, typename U>
  friend constexpr form<U, D1, R1> fmap(const F &f, const form<T1, D1, R1> &x);
  template <typename T1, size_t D1, size_t R1, typename F, typename U>
  friend constexpr form<U, D1, R1> fmap(const F &f, const form<T1, D1, R1> &x,
                                        const form<T1, D1, R1> &y);

  template <typename T1, size_t D1, size_t R1, typename F>
  friend constexpr T1 foldl(const F &f, const T1 &z, const form<T1, D1, R1> &x);

  template <typename T1, size_t D1, size_t R1, typename F, typename G,
            typename U>
  friend constexpr U foldmap(const F &f, const G &g, const U &z,
                             const form<T1, D1, R1> &x);
  template <typename T1, size_t D1, size_t R1, typename F, typename G,
            typename U>
  friend constexpr U foldmap(const F &f, const G &g, const U &z,
                             const form<T1, D1, R1> &x,
                             const form<T1, D1, R1> &y);

  friend /*constexpr*/ form operator+(const form &x) {
    return fmap([](const T &xc) { return +xc; }, x);
  }
  friend /*constexpr*/ form operator-(const form &x) {
    return fmap([](const T &xc) { return -xc; }, x);
  }

  friend /*constexpr*/ form operator+(const form &x, const form &y) {
    return fmap([](const T &xc, const T &yc) { return xc + yc; }, x, y);
  }
  friend /*constexpr*/ form operator-(const form &x, const form &y) {
    return fmap([](const T &xc, const T &yc) { return xc - yc; }, x, y);
  }

  friend /*constexpr*/ form operator*(const T &a, const form &x) {
    return fmap([&](const T &xc) { return a * xc; }, x);
  }
  friend /*constexpr*/ form operator*(const form &x, const T &a) {
    return fmap([&](const T &xc) { return xc * a; }, x);
  }
  friend /*constexpr*/ form operator/(const form &x, const T &a) {
    return fmap([&](const T &xc) { return xc / a; }, x);
  }

  friend constexpr bool operator==(const form &x, const form &y) {
    return fmap([](const T &xc, const T &yc) { return xc == yc; }, x, y).all();
  }
  friend constexpr bool operator!=(const form &x, const form &y) {
    return !(x == y);
  }

  constexpr bool all() const {
    return foldl([](const bool x, const bool y) { return x && y; }, true,
                 *this);
  }

  friend ostream &operator<<(ostream &os, const form &x) {
    os << "form<" << demangle(typeid(T).name()) << "," << D << "," << R
       << ">{comps:[";
    for (size_t n = 0; n < N; ++n) {
      if (n > 0)
        os << ",";
      const size_t idx = comp2idx(n);
      os << "(" << n << "/" << idx << ")";
      for (size_t d = 0; d < D; ++d)
        os << int((idx & (1 << d)) != 0);
      os << ":" << x.comps[n];
    }
    os << "]}";
    return os;
  }
};

template <typename T, size_t D, size_t R>
constexpr form<T, D, R> make_pure(const T &rc) {
  constexpr size_t N = form<T, D, R>::N;
  form<T, D, R> r;
  for (size_t n = 0; n < N; ++n)
    r.comps[n] = rc;
  return r;
}

template <typename T, size_t D, size_t R, typename F,
          typename U = result_of_t<F(T)> >
constexpr form<U, D, R> fmap(const F &f, const form<T, D, R> &x) {
  constexpr size_t N = form<T, D, R>::N;
  form<U, D, R> r;
  for (size_t n = 0; n < N; ++n)
    r.comps[n] = f(x.comps[n]);
  return r;
}

template <typename T, size_t D, size_t R, typename F,
          typename U = result_of_t<F(T, T)> >
constexpr form<U, D, R> fmap(const F &f, const form<T, D, R> &x,
                             const form<T, D, R> &y) {
  constexpr size_t N = form<T, D, R>::N;
  form<U, D, R> r;
  for (size_t n = 0; n < N; ++n)
    r.comps[n] = f(x.comps[n], y.comps[n]);
  return r;
}

template <typename T, size_t D, size_t R, typename F>
constexpr T foldl(const F &f, const T &z, const form<T, D, R> &x) {
  constexpr size_t N = form<T, D, R>::N;
  T r(z);
  for (size_t n = 0; n < N; ++n)
    r = f(r, x.comps[n]);
  return r;
}

template <typename T, size_t D, size_t R, typename F, typename G,
          typename U = result_of_t<F(T)> >
constexpr U foldmap(const F &f, const G &g, const U &z,
                    const form<T, D, R> &x) {
  constexpr size_t N = form<T, D, R>::N;
  U r(z);
  for (size_t n = 0; n < N; ++n)
    r = g(r, f(x.comps[n]));
  return r;
}

template <typename T, size_t D, size_t R, typename F, typename G,
          typename U = result_of_t<F(T, T)> >
constexpr U foldmap(const F &f, const G &g, const U &z, const form<T, D, R> &x,
                    const form<T, D, R> &y) {
  constexpr size_t N = form<T, D, R>::N;
  U r(z);
  for (size_t n = 0; n < N; ++n)
    r = g(r, f(x.comps[n], y.comps[n]));
  return r;
}

template <typename T, size_t D, size_t Rx, size_t Ry, size_t R = Rx + Ry>
constexpr form<T, D, R> wedge(const form<T, D, Rx> &x,
                              const form<T, D, Ry> &y) {
  constexpr size_t Nx = form<T, D, Rx>::N;
  constexpr size_t Ny = form<T, D, Ry>::N;
  form<T, D, R> r;
  for (size_t ny = 0; ny < Ny; ++ny) {
    for (size_t nx = 0; nx < Nx; ++nx) {
      const auto idxx = form<T, D, Rx>::comp2idx(nx);
      const auto idxy = form<T, D, Ry>::comp2idx(ny);
      const bool skip = (idxx & idxy) != 0;
      if (!skip) {
        const size_t idx = idxx | idxy;
        const size_t n = form<T, D, R>::idx2comp(idx);
        bool parity = false;
        bool count = false;
        for (size_t d = 0; d < D; ++d) {
          parity ^= (idxx & (1 << d)) != 0 && count;
          count ^= (idxy & (1 << d)) != 0;
        }
        r.comps[n] += bitsign(parity) * x.comps[nx] * y.comps[ny];
      }
    }
  }
  return r;
}

} // namespace Forms
