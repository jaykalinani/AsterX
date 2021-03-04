#ifndef MEMPOOL_HXX
#define MEMPOOL_HXX

#include <cctk.h>

#include <atomic>
#include <cstddef>
#include <vector>

namespace Loop {
using namespace std;

class mempool_t {
  size_t next_size;

  void *arena;
  size_t arena_size;
  size_t arena_used;

  vector<void *> old_arenas;
  size_t total_size; // arena + all old arenas

public:
  mempool_t();
  mempool_t(const mempool_t &) = delete;
  mempool_t(mempool_t &&mp);
  mempool_t &operator=(const mempool_t &) = delete;
  mempool_t &operator=(mempool_t &&mp);
  friend void swap(mempool_t &mp1, mempool_t &mp2);
  void reset();
  ~mempool_t();

private:
  static void *alloc_arena(size_t count);
  static void free_arena(void *);

  void *restrict alloc_bytes(size_t count);

public:
  template <typename T> T *restrict alloc(size_t count) {
    return reinterpret_cast<T *>(alloc_bytes(count * sizeof(T)));
  }
};

class mempool_set_t {
  vector<mempool_t> mempools;

public:
  mempool_set_t() = default;
  mempool_set_t(const mempool_set_t &) = delete;
  mempool_set_t &operator=(const mempool_set_t &) = delete;

  void reset();

  mempool_t &restrict get_mempool(size_t id);
};

extern mempool_set_t mempools;

} // namespace Loop

#endif // #ifndef MEMPOOL_HXX
