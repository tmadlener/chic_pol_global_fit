#include "global_fit.C"
#include "ArgParser.h"


#if !(defined(__CINT__) or defined(__CLING__))
int main(int argc, char* argv[])
{
  const auto parser = ArgParser(argc, argv);
  const auto nScanPoints1 = parser.getOptionVal<unsigned>("--nscan1", 51);
  const auto nScanPoints2 = parser.getOptionVal<unsigned>("--nscan2", 51);
  const auto outfile = parser.getOptionVal<std::string>("--outfile", "results/scan_file.root");
  const auto flow1 = parser.getOptionVal<double>("--flow1", -1);
  const auto fhigh1 = parser.getOptionVal<double>("--fhigh1", 1);
  const auto flow2 = parser.getOptionVal<double>("--flow2", -1);
  const auto fhigh2 = parser.getOptionVal<double>("--fhigh2", 1);


  global_fit(outfile, nScanPoints1, nScanPoints2, flow1, fhigh1, flow2, fhigh2);
  return 0;
}
#endif
