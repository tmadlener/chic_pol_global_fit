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
  const auto sdcOrder = parser.getOptionVal<std::string>("--order", "NLO");

  std::filesystem::create_directory(outdir);
  if (not paramsSetFile.empty()) {
    const auto paramsFileCopy = outdir + "/" + paramsSetFile;
    if (std::filesystem::exists(paramsFileCopy)) {
      std::filesystem::remove(paramsFileCopy);
    }
    std::filesystem::copy_file(paramsSetFile, paramsFileCopy);
  }

  const auto order = asOrderEnum(sdcOrder);

  global_fit_nrqcd(outdir, order, paramsSetFile);

  return 0;
}
