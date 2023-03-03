#include "logo.hxx"

#include <sstream>
#include <string>
#include <vector>

namespace CarpetX {
using namespace std;
string logo() {

  // buf << "                                                               \n";
  // buf << "    ██████╗ █████╗ ██████╗ ██████╗ ███████╗████████╗██╗  ██╗   \n";
  // buf << "   ██╔════╝██╔══██╗██╔══██╗██╔══██╗██╔════╝╚══██╔══╝╚██╗██╔╝   \n";
  // buf << "   ██║     ███████║██████╔╝██████╔╝█████╗     ██║    ╚███╔╝    \n";
  // buf << "   ██║     ██╔══██║██╔══██╗██╔═══╝ ██╔══╝     ██║    ██╔██╗    \n";
  // buf << "   ╚██████╗██║  ██║██║  ██║██║     ███████╗   ██║   ██╔╝ ██╗   \n";
  // buf << "    ╚═════╝╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚══════╝   ╚═╝   ╚═╝  ╚═╝   \n";
  // buf << "                                                               \n";

  // colours
  const string blk{"\e[30m"};
  const string red{"\e[31m"};
  const string grn{"\e[32m"};
  const string yel{"\e[33m"};
  const string blu{"\e[34m"};
  const string mag{"\e[35m"};
  const string cyn{"\e[36m"};
  const string whi{"\e[37m"};
  const string reset{"\e[39m"};

  // lines and boxes
  const string space{"  "};
  const string vline{"▕▏"};
  const string hline{"──"};
  const string hline2{"━━"};
  const string block{"██"};

  const vector<string> blu_box{
      blu + space + vline + space,
      blu + hline + block + hline,
      blu + space + vline + space,
  };
  const vector<string> grn_box{
      grn + block + block + block,
      grn + block + block + block,
      grn + block + block + block,
  };
  const vector<string> red_box{
      red + block + block + block,
      red + block + block + block,
      red + block + block + block,
  };

  const vector<string> outer_box{
      blu_box[0] + blu_box[0] + blu_box[0],
      blu_box[1] + blu_box[1] + blu_box[1],
      blu_box[2] + blu_box[2] + blu_box[2],
      blu_box[0] + grn_box[0] + blu_box[0],
      blu_box[1] + grn_box[1] + blu_box[1],
      blu_box[2] + grn_box[2] + blu_box[2],
      blu_box[0] + blu_box[0] + blu_box[0],
      blu_box[1] + blu_box[1] + blu_box[1],
      blu_box[2] + blu_box[2] + blu_box[2],
  };
  const vector<string> inner_box{
      red_box[0] + red_box[0] + red_box[0],
      red_box[1] + red_box[1] + red_box[1],
      red_box[2] + red_box[2] + red_box[2],
      red_box[0] + red_box[0] + red_box[0],
      red_box[1] + red_box[1] + red_box[1],
      red_box[2] + red_box[2] + red_box[2],
      red_box[0] + red_box[0] + red_box[0],
      red_box[1] + red_box[1] + red_box[1],
      red_box[2] + red_box[2] + red_box[2],
  };

  vector<string> logo;
  for (int line = 0; line < 9; ++line)
    logo.push_back(outer_box[line] + outer_box[line] + outer_box[line]);
  for (int line = 0; line < 9; ++line)
    logo.push_back(outer_box[line] + inner_box[line] + outer_box[line]);
  for (int line = 0; line < 9; ++line)
    logo.push_back(outer_box[line] + outer_box[line] + outer_box[line]);

  ostringstream buf;
  buf << "\n";
  for (const auto &str : logo)
    buf << "  " << str << reset << "  \n";
  buf << "\n";
  return buf.str();
}
} // namespace CarpetX
