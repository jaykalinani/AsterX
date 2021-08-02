#ifndef SMALLVECTOR_HXX
#define SMALLVECTOR_HXX

#include "defs.hxx"
#include "vect.hxx"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <exception>
#include <initializer_list>
#include <iostream>
#include <limits>

#ifdef CCTK_DEBUG
#define debug_assert(arg) assert(arg)
#else
#define debug_assert(arg) ((void)0)
#endif

namespace Arith {

template <typename T, std::size_t N> class smallvector {
  std::array<T, N> elts;
  std::size_t sz;

public:
  using value_type = typename std::array<T, N>::value_type;
  using size_type = typename std::array<T, N>::size_type;
  using difference_type = typename std::array<T, N>::difference_type;
  using reference = typename std::array<T, N>::reference;
  using const_reference = typename std::array<T, N>::const_reference;
  using pointer = typename std::array<T, N>::pointer;
  using const_pointer = typename std::array<T, N>::const_pointer;
  using iterator = typename std::array<T, N>::iterator;
  using const_iterator = typename std::array<T, N>::const_iterator;
  // using reverse_iterator = typename std::array<T, N>::reverse_iterator;
  // using const_reverse_iterator =
  //     typename std::array<T, N>::const_reverse_iterator;

  smallvector(const smallvector &) = default;
  smallvector(smallvector &&) = default;
  smallvector &operator=(const smallvector &) = default;
  smallvector &operator=(smallvector &&) = default;

  constexpr smallvector() : elts{}, sz(0) {}

  // constexpr smallvector(const std::initializer_list<T> &lst)
  //     : elts(std::array<T, N>(vect<T, N>(lst))), sz(lst.size()) {}
  constexpr smallvector(const std::initializer_list<T> &lst)
      : elts(construct_array<T, N>([&lst](const size_type i) {
          return i < lst.size() ? lst.begin()[i] : T{};
        })),
        sz(lst.size()) {}
  template <size_type N1>
  constexpr smallvector(const std::array<T, N1> &elts_)
      : elts(construct_array<T, N>(
            [&elts_](const size_type i) { return i < N1 ? elts_[i] : T{}; })),
        sz(N1) {}

  constexpr bool empty() const { return sz == 0; }
  constexpr size_type size() const { return sz; }
  constexpr size_type max_size() const { return N; }

  void resize(const size_type newsz) {
    debug_assert(newsz <= N);
    sz = newsz;
  }

  const_reference at(const size_type i) const {
    if (i >= sz)
      throw std::out_of_range();
    return elts[i];
  }
  reference at(const size_type i) {
    if (i >= sz)
      throw std::out_of_range();
    return elts[i];
  }
  constexpr const_reference operator[](const size_type i) const {
    debug_assert(i < sz);
    return elts[i];
  }
  constexpr reference operator[](const size_type i) {
    debug_assert(i < sz);
    return elts[i];
  }

  constexpr const_reference front() const {
    debug_assert(sz > 0);
    return elts[0];
  }
  constexpr reference front() {
    debug_assert(sz > 0);
    return elts[0];
  }
  constexpr const_reference back() const {
    debug_assert(sz > 0);
    return elts[sz - 1];
  }
  constexpr reference back() {
    debug_assert(sz > 0);
    return elts[sz - 1];
  }

  const T *data() const { return elts.data(); }
  T *data() { return elts.data(); }

  const_iterator begin() const { return &elts[0]; }
  iterator begin() { return &elts[0]; }
  const_iterator end() const { return &elts[sz]; }
  iterator end() { return &elts[sz]; }

  // constexpr const T *rbegin() const { return &elts[0]; }
  // constexpr T *rbegin() { return &elts[0]; }
  // constexpr const T *rend() const { return &elts[sz]; }
  // constexpr T *rend() { return &elts[sz]; }

  void clear() { sz = 0; }

  void push_back(const T &x) {
    debug_assert(sz < N);
    elts[sz++] = x;
  }
  void push_back(T &&x) {
    debug_assert(sz < N);
    elts[sz++] = x;
  }
  void pop_back() {
    debug_assert(sz > 0);
    --sz;
  }

  template <typename... Args> reference emplace_back(Args &&...args) {
    debug_assert(sz < N);
    return elts[sz++] = T(std::forward<Args>(args)...);
  }

  friend constexpr smallvector operator+(const smallvector &x) {
    smallvector r;
    r.sz = x.sz;
    for (size_type i = 0; i < r.sz; ++i)
      r.elts[i] = +x.elts[i];
    return r;
  }
  friend constexpr smallvector operator-(const smallvector &x) {
    smallvector r;
    r.sz = x.sz;
    for (size_type i = 0; i < r.sz; ++i)
      r.elts[i] = -x.elts[i];
    return r;
  }

  friend constexpr smallvector operator+(const smallvector &x,
                                         const smallvector &y) {
    smallvector r;
    debug_assert(x.sz == y.sz);
    using std::min;
    r.sz = min(x.sz, y.sz);
    for (size_type i = 0; i < r.sz; ++i)
      r.elts[i] = x.elts[i] + y.elts[i];
    return r;
  }
  friend constexpr smallvector operator-(const smallvector &x,
                                         const smallvector &y) {
    smallvector r;
    debug_assert(x.sz == y.sz);
    using std::min;
    r.sz = min(x.sz, y.sz);
    for (size_type i = 0; i < r.sz; ++i)
      r.elts[i] = x.elts[i] - y.elts[i];
    return r;
  }

  friend constexpr smallvector operator*(const T &a, const smallvector &y) {
    smallvector r;
    r.sz = y.sz;
    for (size_type i = 0; i < r.sz; ++i)
      r.elts[i] = a * y.elts[i];
    return r;
  }
  friend constexpr smallvector operator*(const smallvector &x, const T &b) {
    smallvector r;
    r.sz = x.sz;
    for (size_type i = 0; i < r.sz; ++i)
      r.elts[i] = x.elts[i] * b;
    return r;
  }
  friend constexpr smallvector operator/(const smallvector &x, const T &b) {
    smallvector r;
    r.sz = x.sz;
    for (size_type i = 0; i < r.sz; ++i)
      r.elts[i] = x.elts[i] / b;
    return r;
  }
  template <typename T1 = T,
            typename = std::enable_if_t<std::is_integral_v<T1> > >
  friend constexpr smallvector operator%(const smallvector &x, const T &b) {
    smallvector r;
    r.sz = x.sz;
    for (size_type i = 0; i < r.sz; ++i)
      r.elts[i] = x.elts[i] % b;
    return r;
  }

  friend constexpr smallvector abs(const smallvector &x) {
    smallvector r;
    r.sz = x.sz;
    using std::abs;
    for (size_type i = 0; i < r.sz; ++i)
      r.elts[i] = abs(x.elts[i]);
    return r;
  }

  friend constexpr T maximum(const smallvector &x) {
    T r = -std::numeric_limits<T>::infinity();
    for (size_type i = 0; i < x.sz; ++i)
      r = max1(r, x.elts[i]);
    return r;
  }
  friend constexpr T minimum(const smallvector &x) {
    T r = std::numeric_limits<T>::infinity();
    for (size_type i = 0; i < x.sz; ++i)
      r = min1(r, x.elts[i]);
    return r;
  }
  friend constexpr T prod(const smallvector &x) {
    T r = 1;
    for (size_type i = 0; i < x.sz; ++i)
      r *= x.elts[i];
    return r;
  }
  friend constexpr T sum(const smallvector &x) {
    T r = 0;
    for (size_type i = 0; i < x.sz; ++i)
      r += x.elts[i];
    return r;
  }

  friend constexpr bool operator==(const smallvector &x, const smallvector &y) {
    if (x.sz != y.sz)
      return false;
    for (size_type i = 0; i < x.sz; ++i)
      if (x.elts[i] != y.elts[i])
        return false;
    return true;
  }
  friend constexpr bool operator!=(const smallvector &x, const smallvector &y) {
    return !(x == y);
  }
  friend constexpr bool operator<(const smallvector &x, const smallvector &y) {
    using std::min;
    for (size_type i = 0; i < min(x.sz, y.sz); ++i) {
      if (x.elts[i] < y.elts[i])
        return true;
      if (y.elts[i] < x.elts[i])
        return false;
    }
    return x.sz < y.sz;
  }
  friend constexpr bool operator>(const smallvector &x, const smallvector &y) {
    return y < x;
  }
  friend constexpr bool operator<=(const smallvector &x, const smallvector &y) {
    return !(y < x);
  }
  friend constexpr bool operator>=(const smallvector &x, const smallvector &y) {
    return !(x < y);
  }

  friend std::ostream &operator<<(std::ostream &os, const smallvector &x) {
    os << "[";
    for (const auto &a : x)
      os << a << ",";
    os << "]";
    return os;
  }
};

} // namespace Arith

#endif // #ifndef SMALLVECTOR_HXX
