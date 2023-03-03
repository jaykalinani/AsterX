#include "forms.hxx"

namespace Forms {

void q() {
  form<double, 3, 0> a0({0});
  form<double, 3, 1> a1({1, 2, 3});
  form<double, 3, 2> a2({4, 5, 6});
  form<double, 3, 3> a3({7});

  +a0;
  +a1;
  +a2;
  +a3;

  -a0;
  -a1;
  -a2;
  -a3;

  a0 + a0;
  a1 + a1;
  a2 + a2;
  a3 + a3;

  a0 - a0;
  a1 - a1;
  a2 - a2;
  a3 - a3;

  2 * a0;
  2 * a1;
  2 * a2;
  2 * a3;

  a0 * 2;
  a1 * 2;
  a2 * 2;
  a3 * 2;

  a0 / 2;
  a1 / 2;
  a2 / 2;
  a3 / 2;

  wedge(a0, a0);
  wedge(a0, a1);
  wedge(a0, a2);
  wedge(a0, a3);
  wedge(a1, a0);
  wedge(a1, a1);
  wedge(a1, a2);
  wedge(a2, a0);
  wedge(a2, a1);
  wedge(a3, a0);
}

} // namespace Forms
