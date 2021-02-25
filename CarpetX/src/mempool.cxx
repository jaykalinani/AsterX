#include "mempool.hxx"

#ifdef _OPENMP
#include <omp.h>
#else
extern "C" {
static inline int omp_get_max_threads(void) { return 1; }
static inline int omp_get_num_threads(void) { return 1; }
static inline int omp_get_thread_num(void) { return 0; }
static inline int omp_in_parallel(void) { return 0; }
}
#endif

#include <algorithm>
#include <cassert>

namespace Loop {

const size_t cacheline_size = 64;

mempool_t::mempool_t() : arena_size(1024 * 1024), total_size(0) {}

mempool_t::mempool_t(mempool_t &&mp) : mempool_t() { swap(*this, mp); }

mempool_t &mempool_t::operator=(mempool_t &&mp) {
  {
    mempool_t tmp;
    swap(*this, tmp);
  }
  swap(*this, mp);
  return *this;
}

void swap(mempool_t &mp1, mempool_t &mp2) {
  using std::swap;
  swap(mp1.arena_size, mp2.arena_size);
  swap(mp1.arena, mp2.arena);
  swap(mp1.old_arenas, mp2.old_arenas);
  swap(mp1.total_size, mp2.total_size);
}

void mempool_t::reset() {
  old_arenas.clear();
  arena_size = max(arena_size, total_size);
  total_size = 0;
  arena.clear();
}

unsigned char *restrict mempool_t::alloc_bytes(size_t count) {
  // Handle empty objects
  if (count == 0)
    return nullptr;
  // Align count to cache line size
  count = (count + cacheline_size - 1) / cacheline_size * cacheline_size;
  // Ensure the new object fits into the arena
  if (!(arena.size() + count <= arena.capacity())) {
    if (!arena.empty())
      old_arenas.push_back(move(arena));
    arena_size = max(count, 2 * arena_size);
    arena.reserve(arena_size);
  }
  // Check alignment
  const auto ptr = arena.data() + arena.size();
  if (uintptr_t(ptr) % cacheline_size != 0)
    CCTK_VERROR("A large system memory allocation returned an unaligned "
                "pointer (pointer=0x%p, required alignment=%td)",
                arena.data(), cacheline_size);
  // Allocate
  assert(arena.size() + count <= arena.capacity());
  auto old_data = arena.data();
  auto old_size = arena.size();
  arena.resize(arena.size() + count);
  assert(arena.data() == old_data);
  assert(arena.size() == old_size + count);
  assert(ptr + count <= arena.data() + arena.size());
  total_size += count;
  return ptr;
}

// TODO: Keep a global mempool here, not in the callers
mempool_t &restrict mempool_set_t::get_mempool() {
#pragma omp critical(CarpetX_mempool_get)
  if (!have_mempools) {
    const int max_threads = omp_get_max_threads();
    mempools.resize(max_threads);
    have_mempools = true;
  }
  const int thread_num = omp_get_thread_num();
  mempool_t &restrict mempool = mempools.at(thread_num);
  mempool.reset(); // TODO: reset at end of calling function, not here
  return mempool;
}

} // namespace Loop
