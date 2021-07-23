#ifndef DO_COMBINATION_FIT
// Combination fit: Instead of using the inputs with a fixed order, use the NLO
// SDCs as baseline and add a (configurable) fraction of the LP correction on
// top of it
#define DO_COMBINATION_FIT 1
#endif

#include "ArgParser.h"
#include "global_fit_nrqcd.C"

#include <iostream>
#include <string>
#include <filesystem>
#include <iostream>
#include <cstdlib>


int main(int argc, char* argv[]) {
  const auto parser = ArgParser(argc, argv);
  const auto outdir = parser.getOptionVal<std::string>("--outdir", "results/");
  const auto paramsSetFile = parser.getOptionVal<std::string>("--paramsSettings", "");
#if DO_COMBINATION_FIT
  const auto lpFraction = parser.getOptionVal<double>("--lpfraction", 0);
#else
  const auto sdcOrder = parser.getOptionVal<std::string>("--order", "NLO");
#endif
  const auto sdcDir = parser.getOptionVal<std::string>("--sdcdir", "../SDCs_Chung/");
  const auto datadir = parser.getOptionVal<std::string>("--datadir", "./data/");

  std::filesystem::create_directory(outdir);
  if (not paramsSetFile.empty()) {
    // Normalize the resulting file name here for easier scripting
    const auto paramsFileCopy = outdir + "/params_settings.txt";
    if (std::filesystem::exists(paramsFileCopy)) {
      std::filesystem::remove(paramsFileCopy);
    }
    std::filesystem::copy_file(paramsSetFile, paramsFileCopy);
  }

#if DO_COMBINATION_FIT
  global_fit_nrqcd(outdir, datadir, sdcDir, lpFraction, paramsSetFile);
#else
  const auto order = asOrderEnum(sdcOrder);
  global_fit_nrqcd(outdir, datadir, sdcDir, order, paramsSetFile);
#endif

  return 0;
}
