#ifndef VALID_HXX
#define VALID_HXX

#include <yaml-cpp/yaml.h>

#include <functional>
#include <iostream>
#include <sstream>
#include <tuple>

namespace CarpetX {
using namespace std;

struct valid_t {
  bool valid_int, valid_outer, valid_ghosts;
  constexpr valid_t() : valid_t(false) {}

  explicit constexpr valid_t(bool b)
      : valid_int(b), valid_outer(b), valid_ghosts(b) {}
  friend constexpr valid_t operator~(const valid_t &x) {
    valid_t r;
    r.valid_int = !x.valid_int;
    r.valid_outer = !x.valid_outer;
    r.valid_ghosts = !x.valid_ghosts;
    return r;
  }
  friend constexpr valid_t operator&(const valid_t &x, const valid_t &y) {
    valid_t r;
    r.valid_int = x.valid_int && y.valid_int;
    r.valid_outer = x.valid_outer && y.valid_outer;
    r.valid_ghosts = x.valid_ghosts && y.valid_ghosts;
    return r;
  }
  friend constexpr valid_t operator|(const valid_t &x, const valid_t &y) {
    valid_t r;
    r.valid_int = x.valid_int || y.valid_int;
    r.valid_outer = x.valid_outer || y.valid_outer;
    r.valid_ghosts = x.valid_ghosts || y.valid_ghosts;
    return r;
  }
  valid_t &operator&=(const valid_t &x) { return *this = *this & x; }
  valid_t &operator|=(const valid_t &x) { return *this = *this | x; }

  constexpr bool valid_all() const {
    return valid_int && valid_outer && valid_ghosts;
  }
  constexpr bool valid_any() const {
    return valid_int || valid_outer || valid_ghosts;
  }

  friend bool operator==(const valid_t &x, const valid_t &y) {
    return make_tuple(x.valid_int, x.valid_outer, x.valid_ghosts) ==
           make_tuple(y.valid_int, y.valid_outer, y.valid_ghosts);
  }
  friend bool operator<(const valid_t &x, const valid_t &y) {
    return make_tuple(x.valid_int, x.valid_outer, x.valid_ghosts) <
           make_tuple(y.valid_int, y.valid_outer, y.valid_ghosts);
  }

  // TODO: Move I/O functions to valid.cxx
  friend ostream &operator<<(ostream &os, const valid_t v) {
    auto str = [](bool v) { return v ? "VAL" : "INV"; };
    return os << "[int:" << str(v.valid_int) << ",outer:" << str(v.valid_outer)
              << ",ghosts:" << str(v.valid_ghosts) << "]";
  }
  operator string() const {
    ostringstream buf;
    buf << *this;
    return buf.str();
  }

  friend YAML::Emitter &operator<<(YAML::Emitter &yaml, const valid_t v) {
    yaml << YAML::LocalTag("valid-1.0.0");
    yaml << YAML::Flow << YAML::BeginMap;
    yaml << YAML::Key << "int" << YAML::Value << v.valid_int;
    yaml << YAML::Key << "outer" << YAML::Value << v.valid_outer;
    yaml << YAML::Key << "ghosts" << YAML::Value << v.valid_ghosts;
    yaml << YAML::EndMap;
    return yaml;
  }
};

inline constexpr valid_t make_valid_int() {
  valid_t valid;
  valid.valid_int = true;
  return valid;
}
inline constexpr valid_t make_valid_outer() {
  valid_t valid;
  valid.valid_outer = true;
  return valid;
}
inline constexpr valid_t make_valid_ghosts() {
  valid_t valid;
  valid.valid_ghosts = true;
  return valid;
}
inline constexpr valid_t make_valid_all() { return ~valid_t(); }

} // namespace CarpetX
namespace std {
using namespace CarpetX;
template <> struct equal_to<valid_t> {
  constexpr bool operator()(const valid_t &x, const valid_t &y) const {
    return make_tuple(x.valid_int, x.valid_outer, x.valid_ghosts) ==
           make_tuple(y.valid_int, y.valid_outer, y.valid_ghosts);
  }
};
template <> struct less<valid_t> {
  constexpr bool operator()(const valid_t &x, const valid_t &y) const {
    return make_tuple(x.valid_int, x.valid_outer, x.valid_ghosts) <
           make_tuple(y.valid_int, y.valid_outer, y.valid_ghosts);
  }
};
} // namespace std
namespace CarpetX {

class why_valid_t {
  valid_t valid;
  function<string()> why_int, why_outer, why_ghosts;

public:
  // The constructor that doesn't give a reason should never be called
  why_valid_t() = delete;
  why_valid_t(const function<string()> &why) : why_valid_t(false, why) {}
  why_valid_t(bool b, const function<string()> &why)
      : why_valid_t(valid_t(b), why) {}
  why_valid_t(const valid_t &val, const function<string()> &why)
      : valid(val), why_int(why), why_outer(why), why_ghosts(why) {}

  const valid_t &get() const { return valid; }

  void set(const valid_t &val, const function<string()> &why) {
    auto old_valid = valid;
    valid = val;
    auto new_valid = valid;
    if (new_valid.valid_int != old_valid.valid_int)
      why_int = why;
    if (new_valid.valid_outer != old_valid.valid_outer)
      why_outer = why;
    if (new_valid.valid_ghosts != old_valid.valid_ghosts)
      why_ghosts = why;
  }
  void set_int(bool b, const function<string()> &why) {
    auto val = valid;
    val.valid_int = b;
    set(val, why);
  }
  void set_outer(bool b, const function<string()> &why) {
    auto val = valid;
    val.valid_outer = b;
    set(val, why);
  }
  void set_ghosts(bool b, const function<string()> &why) {
    auto val = valid;
    val.valid_ghosts = b;
    set(val, why);
  }
  void set_and(const valid_t &val, const function<string()> &why) {
    set(valid & val, why);
  }
  void set_or(const valid_t &val, const function<string()> &why) {
    set(valid | val, why);
  }

  friend ostream &operator<<(ostream &os, const why_valid_t &why) {
    return os << why.valid << ","
              << "why{int:" << why.why_int() << ","
              << "outer:" << why.why_outer() << ","
              << "ghosts:" << why.why_ghosts() << "}";
  }
  operator string() const {
    ostringstream buf;
    buf << *this;
    return buf.str();
  }

  friend YAML::Emitter &operator<<(YAML::Emitter &yaml,
                                   const why_valid_t &why) {
    yaml << YAML::LocalTag("why_valid-1.0.0");
    yaml << YAML::BeginMap;
    yaml << YAML::Key << "valid" << YAML::Value << why.valid;
    yaml << YAML::Key << "why" << YAML::Value << YAML::BeginMap;
    yaml << YAML::Key << "int" << YAML::Value << why.why_int();
    yaml << YAML::Key << "outer" << YAML::Value << why.why_outer();
    yaml << YAML::Key << "ghosts" << YAML::Value << why.why_ghosts();
    yaml << YAML::EndMap;
    yaml << YAML::EndMap;
    return yaml;
  }
};

} // namespace CarpetX

#endif // #ifndef VALID_HXX
