#include "global_fit.C"
#include "ArgParser.h"

#include <iostream>

#if !(defined(__CINT__) or defined(__CLING__))
int main(int argc, char* argv[])
{
  const auto parser = ArgParser(argc, argv);
  auto nSamplePoints = parser.getOptionVal<unsigned>("--npoints", 1000000);
  const auto outfile = parser.getOptionVal<std::string>("--outfile", "results/scan_file.root");
  const auto graphfile = parser.getOptionVal<std::string>("--graphfile", "results/fit_result_graphs.root");
  const auto useCosthRatios = parser.getOptionVal<bool>("--useCosthRatios", true);
  const auto noscan = parser.getOptionVal<bool>("--noscan", false);
  const auto nographs = parser.getOptionVal<bool>("--nographs", false);

  if (noscan) {
    nSamplePoints = 0;
  }

  global_fit(outfile, graphfile, nSamplePoints, useCosthRatios, !nographs);
  return 0;
}
#endif
