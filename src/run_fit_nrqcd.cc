#include "misc_utils.h"

#include "ArgParser.h"
#include "global_fit_nrqcd.C"

#include <iostream>
#include <string>
#include <filesystem>
#include <iostream>
#include <cstdlib>

sdc::SDCType asOrderEnum(const std::string& sdcOrder) {
  const auto it = std::find(sdc::SDCTypeNames.begin(), sdc::SDCTypeNames.end(), sdcOrder);
  if (it == sdc::SDCTypeNames.end()) {
    std::cerr << "Argument to --order is not valid. Must be one of \'LO\', \'NLO\' or \'LP+NLO\'" << std::endl;
    std::exit(1);
  }

  return util::to_enum<sdc::SDCType>(std::distance(sdc::SDCTypeNames.begin(), it));
}

int main(int argc, char* argv[]) {
  const auto parser = ArgParser(argc, argv);
  const auto outdir = parser.getOptionVal<std::string>("--outdir", "results/");
  const auto paramsSetFile = parser.getOptionVal<std::string>("--paramsSettings", "");
  const auto sdcOrder = parser.getOptionVal<std::string>("--order", "NLO");

  const std::string outFile = outdir + "/fit_results_nrqcd_global_fit.root";
  const std::string graphFile = outdir + "/fit_graphs_and_models_nrqcd_global_fit.root";

  std::filesystem::create_directory(outdir);
  if (not paramsSetFile.empty()) {
    const auto paramsFileCopy = outdir + "/" + paramsSetFile;
    if (std::filesystem::exists(paramsFileCopy)) {
      std::filesystem::remove(paramsFileCopy);
    }
    std::filesystem::copy_file(paramsSetFile, paramsFileCopy);
  }

  const auto order = asOrderEnum(sdcOrder);

  global_fit_nrqcd(outFile, graphFile, order, paramsSetFile);

  return 0;
}
