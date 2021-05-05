#include "ArgParser.h"
#include "global_fit_nrqcd.C"

#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
  const auto parser = ArgParser(argc, argv);
  const auto outdir = parser.getOptionVal<std::string>("--outdir", "results/");
  const auto paramsSetFile = parser.getOptionVal<std::string>("--paramsSettings", "");

  const std::string outFile = outdir + "/fit_results_nrqcd_global_fit.root";
  const std::string graphFile = outdir + "/fit_graphs_and_models_nrqcd_global_fit.root";

  global_fit_nrqcd(outFile, graphFile, paramsSetFile);

  return 0;
}
