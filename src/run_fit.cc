#include "global_fit.C"
#include "ArgParser.h"

#include <iostream>

#if !(defined(__CINT__) or defined(__CLING__))
int main(int argc, char* argv[])
{
  const auto parser = ArgParser(argc, argv);
  const auto nSamplePoints = parser.getOptionVal<unsigned>("--npoints", 1000000);
  const auto outfile = parser.getOptionVal<std::string>("--outfile", "results/scan_file.root");
  const auto graphfile = parser.getOptionVal<std::string>("--graphfile", "results/fit_result_graphs.root");
  const auto datadir = parser.getOptionVal<std::string>("--datadir", "./data/");
  const auto useCosthRatios = parser.getOptionVal<bool>("--useCosthRatios", true);
  const auto noscan = parser.getOptionVal<bool>("--noscan", false);
  const auto nographs = parser.getOptionVal<bool>("--nographs", false);
  const auto nSamplePoints1 = parser.getOptionVal<unsigned>("--nscan1", 51);
  const auto nSamplePoints2 = parser.getOptionVal<unsigned>("--nscan2", 51);
  const auto flow1 = parser.getOptionVal<double>("--flow1", -1);
  const auto fhigh1 = parser.getOptionVal<double>("--fhigh1", 1);
  const auto flow2 = parser.getOptionVal<double>("--flow2", -1);
  const auto fhigh2 = parser.getOptionVal<double>("--fhigh2", 1);
  const auto randomScan = parser.getOptionVal<bool>("--randomScan", false);
  const auto usePsiPolarizations = parser.getOptionVal<bool>("--usePsiPolarizations", true);
  const auto usePhysicalCorrections = parser.getOptionVal<bool>("--usePhysicalCorrections", false);
  const auto fixLambdasChicBestFit = parser.getOptionVal<bool>("--fixLambdasChicBestFit", false);
  const auto minPtM = parser.getOptionVal<double>("--minPtM", 2.0);


  ScanArguments scanArgs;
  scanArgs.randomScan = randomScan;
  scanArgs.nSamplePoints = nSamplePoints;
  scanArgs.nScan1 = nSamplePoints1;
  scanArgs.nScan2 = nSamplePoints2;
  scanArgs.flow1 = flow1;
  scanArgs.fhigh1 = fhigh1;
  scanArgs.flow2 = flow2;
  scanArgs.fhigh2 = fhigh2;

  if (noscan) {
    scanArgs.nSamplePoints = 0;
    scanArgs.nScan1 = 0;
    scanArgs.nScan2 = 0;
  }

  if (scanArgs.randomScan) {
    scanArgs.nScan1 = 0;
    scanArgs.nScan2 = 0;
  }

  global_fit(outfile, graphfile, datadir, scanArgs, useCosthRatios, !nographs, usePsiPolarizations,
             usePhysicalCorrections, fixLambdasChicBestFit, minPtM);
  return 0;
}
#endif
