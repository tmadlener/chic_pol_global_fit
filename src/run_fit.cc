#include "global_fit.C"
#include "ArgParser.h"


#if !(defined(__CINT__) or defined(__CLING__))
int main(int argc, char* argv[])
{
  const auto parser = ArgParser(argc, argv);
  const auto nScanPoints = parser.getOptionVal<unsigned>("--nscan", 5000000);
  const auto outfile = parser.getOptionVal<std::string>("--outfile", "results/scan_file.root");


  global_fit(nScanPoints, outfile);
  return 0;
}
#endif
