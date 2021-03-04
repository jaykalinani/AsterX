#include "mempool.hxx"

#include <cctk.h>

#include <AMReX.H>
#include <AMReX_Arena.H>

#include <algorithm>
#include <cassert>
#include <stdlib.h>

namespace Loop {

#ifndef AMREX_USE_GPU
const size_t cacheline_size = 64; // CPU
#else
const size_t cacheline_size = 256; // GPU
#endif

mempool_t::mempool_t()
    : next_size(1024 * 1024), arena(nullptr), arena_size(0), arena_used(0),
      total_size(0) {}

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
  swap(mp1.next_size, mp2.next_size);
  swap(mp1.arena, mp2.arena);
  swap(mp1.arena_size, mp2.arena_size);
  swap(mp1.arena_used, mp2.arena_used);
  swap(mp1.old_arenas, mp2.old_arenas);
  swap(mp1.total_size, mp2.total_size);
}

void mempool_t::reset() {
  next_size = max(next_size / 2, total_size);

  if (arena) {
    free_arena(arena);
    arena = nullptr;
  }
  arena_size = 0;
  arena_used = 0;

  for (void *old_arena : old_arenas)
    free_arena(old_arena);
  old_arenas.clear();

  total_size = 0;
}

mempool_t::~mempool_t() { reset(); }

#ifndef AMREX_USE_GPU
// CPU

void *mempool_t::alloc_arena(size_t count) {
  // const auto ptr = new unsigned char[count];
  // // Check alignment
  // if (uintptr_t(ptr) % cacheline_size != 0)
  //   CCTK_VERROR("A large system memory allocation returned an unaligned "
  //               "pointer (pointer=0x%p, required alignment=%td)",
  //               ptr, cacheline_size);
  void *ptr;
  int ierr = posix_memalign(&ptr, cacheline_size, count);
  if (ierr) {
    CCTK_VINFO("posix_memalign(...,%td,%td)=%d", cacheline_size, count, ierr);
    ptr = malloc(count);
    if (!ptr)
      CCTK_VINFO("malloc(%td) failed", count);
  }
  return ptr;
}

void mempool_t::free_arena(void *arena) {
  // delete[] arena;
  free(arena);
}

#else
// GPU

void *mempool_t::alloc_arena(size_t count) {
  return amrex::The_Arena()->alloc(count);
}

void mempool_t::free_arena(void *arena) { amrex::The_Arena()->free(arena); }

#endif

void *restrict mempool_t::alloc_bytes(size_t count) {
  // Handle empty objects
  if (count == 0)
    return nullptr;
  // Align count to cache line size
  count = (count + cacheline_size - 1) / cacheline_size * cacheline_size;
  // Ensure the new object fits into the arena
  if (!(arena_used + count <= arena_size)) {
    if (arena)
      old_arenas.push_back(arena);
    next_size = max(next_size, count);
    arena_size = next_size;
    arena = alloc_arena(arena_size);
    arena_used = 0;
    next_size *= 2;
  }
  // Allocate
  assert(arena_used + count <= arena_size);
  const auto ptr = static_cast<unsigned char *>(arena) + arena_used;
  arena_used += count;
  total_size += count;
  return ptr;
}

// TODO: Keep a global mempool here, not in the callers
mempool_t &restrict mempool_set_t::get_mempool(const size_t id) {
  assert(id >= 0);
  mempool_t *pmempool = nullptr;
#pragma omp critical(CarpetX_mempool_get)
  {
    if (id >= mempools.size()) {
      if (id >= mempools.capacity())
        mempools.reserve(max(id + 1, 2 * mempools.capacity() + 100));
      mempools.resize(id + 1);
    }
    pmempool = &mempools.at(id);
  }
  pmempool->reset();
  return *pmempool;
}

void mempool_set_t::reset() { mempools.clear(); }

mempool_set_t mempools;

} // namespace Loop
