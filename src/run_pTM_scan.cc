#ifndef NRQCD_FIT
#define NRQCD_FIT 0
#endif

#if NRQCD_FIT
#ifndef DO_COMBINATION_FIT
#define DO_COMBINATION_FIT 1
#endif
#include "read_sdcs.h"
#endif

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
  const auto sdcDir = parser.getOptionVal<std::string>("--sdcdir");

#if NRQCD_FIT
#if DO_COMBINATION_FIT
  const auto lpfraction = parser.getOptionVal<double>("--lpfraction", 0);

  pTM_scan(fitresult, outfile, ptmmin, ptmmax, npoints, nscans, sdcDir, !noparams, lpfraction);
#else
  const auto sdcOrder = parser.getOptionVal<std::string>("--order", "NLO");

  const auto order = asOrderEnum(sdcOrder);
  pTM_scan(fitresult, outfile, ptmmin, ptmmax, npoints, nscans, sdcDir, !noparams, order);
#endif // DO_COMBINATION_FIT
#else
  pTM_scan(fitresult, outfile, ptmmin, ptmmax, npoints, nscans, !noparams);
#endif

  return 0;
}
#endif
