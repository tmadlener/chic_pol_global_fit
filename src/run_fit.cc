#include "global_fit.C"
#include "ArgParser.h"

#include <iostream>

#if !(defined(__CINT__) or defined(__CLING__))
int main(int argc, char* argv[])
{
  const auto parser = ArgParser(argc, argv);
  auto nScanPoints1 = parser.getOptionVal<unsigned>("--nscan1", 51);
  auto nScanPoints2 = parser.getOptionVal<unsigned>("--nscan2", 51);
  const auto outfile = parser.getOptionVal<std::string>("--outfile", "results/scan_file.root");
  const auto graphfile = parser.getOptionVal<std::string>("--graphfile", "results/fit_result_graphs.root");
  const auto flow1 = parser.getOptionVal<double>("--flow1", -1);
  const auto fhigh1 = parser.getOptionVal<double>("--fhigh1", 1);
  const auto flow2 = parser.getOptionVal<double>("--flow2", -1);
  const auto fhigh2 = parser.getOptionVal<double>("--fhigh2", 1);
  const auto useCosthRatios = parser.getOptionVal<bool>("--useCosthRatios", true);
  const auto noscan = parser.getOptionVal<bool>("--noscan", false);
  const auto nographs = parser.getOptionVal<bool>("--nographs", false);

  if (noscan) {
    nScanPoints1 = 0;
    nScanPoints2 = 0;
  }

  global_fit(outfile, graphfile,
             nScanPoints1, nScanPoints2,
             flow1, fhigh1, flow2, fhigh2, useCosthRatios, !nographs);
  return 0;
}
#endif
