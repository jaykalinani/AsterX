#ifndef MEMPOOL_HXX
#define MEMPOOL_HXX

#include <cctk.h>

#include <atomic>
#include <cstddef>
#include <vector>

namespace Loop {
using namespace std;

class mempool_t {
  size_t arena_size;
  vector<unsigned char> arena;
  vector<vector<unsigned char> > old_arenas;
  size_t total_size;

public:
  mempool_t();
  mempool_t(const mempool_t &) = delete;
  mempool_t(mempool_t &&mp);
  mempool_t &operator=(const mempool_t &) = delete;
  mempool_t &operator=(mempool_t &&mp);
  friend void swap(mempool_t &mp1, mempool_t &mp2);
  void reset();

private:
  unsigned char *restrict alloc_bytes(size_t count);

public:
  template <typename T> T *restrict alloc(size_t count) {
    return reinterpret_cast<T *>(alloc_bytes(count * sizeof(T)));
  }
};

class mempool_set_t {
  atomic<bool> have_mempools;
  vector<mempool_t> mempools;

public:
  mempool_set_t() : have_mempools{false} {}
  mempool_t &restrict get_mempool();
};

} // namespace Loop

#endif // #ifndef MEMPOOL_HXX
