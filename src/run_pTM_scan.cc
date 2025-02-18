#include "pTM_scan.C"
#include "ArgParser.h"

#if !(defined(__CINT__) or defined(__CLING__))
int main(int argc, char* argv[])
{
  const auto parser = ArgParser(argc, argv);
  const auto fitresult = parser.getOptionVal<std::string>("--fitresult");
  const auto outfile = parser.getOptionVal<std::string>("--outfile", "ptm_dep_scan_results.root");
  const auto ptmmin = parser.getOptionVal<double>("--ptmmin", 0);
  const auto ptmmax = parser.getOptionVal<double>("--ptmmax", 10.0);
  const auto npoints = parser.getOptionVal<size_t>("--npoints", 50);
  const auto nscans = parser.getOptionVal<size_t>("--nscans", 20000);
  const auto noparams = parser.getOptionVal<bool>("--noparams", false);


  pTM_scan(fitresult, outfile, ptmmin, ptmmax, npoints, nscans, !noparams);

  return 0;
}
#endif
