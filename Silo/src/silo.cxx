#include "silo.hxx"

#include <cctype>
#include <sstream>

namespace DB {
using namespace std;

string legalize_name(const string &name) {
  ostringstream buf;
  for (const char c : name) {
    if (isalnum(c) || c == '_')
      buf << c;
    else
      buf << '_';
  }
  return buf.str();
}

} // namespace DB
