#ifndef SILO_HXX
#define SILO_HXX

#include <silo.h>

#include <cstdlib>
#include <memory>
#include <string>
#include <type_traits>

namespace DB {

namespace detail {
inline void free(DBfile *const obj) { DBClose(obj); }
inline void free(DBgroupelmap *const obj) { DBFreeGroupelmap(obj); }
inline void free(DBmrgtree *const obj) { DBFreeMrgtree(obj); }
inline void free(DBmrgvar *const obj) { DBFreeMrgvar(obj); }
inline void free(DBmultimesh *const obj) { DBFreeMultimesh(obj); }
inline void free(DBmultivar *const obj) { DBFreeMultivar(obj); }
inline void free(DBoptlist *const obj) { DBFreeOptlist(obj); }
inline void free(DBquadmesh *const obj) { DBFreeQuadmesh(obj); }
inline void free(DBquadvar *const obj) { DBFreeQuadvar(obj); }
inline void free(DBtoc *const obj) {}
} // namespace detail

template <typename T> using ptr = std::shared_ptr<T>;

template <typename T> ptr<T> make(T *const obj) {
  return std::shared_ptr<T>(
      obj,
      static_cast<void (*)(typename std::remove_cv<T>::type *)>(detail::free));
}

std::string legalize_name(const std::string &name);

} // namespace DB

#endif // #ifndef SILO_HXX
