#include "ArgParser.h"

#include "costh_fit.C"

#include <string>

#if !(defined(__CINT__) or defined(__CLING__))
int main(int argc, char* argv[])
{
  ArgParser parser(argc, argv);
  const auto outfile = parser.getOptionVal<std::string>("--outfile", "results/costh_scan.root");

  costh_fit(outfile);
  return 0;
}
#endif
