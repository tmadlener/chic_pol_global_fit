#include "data_reading.h"
#include "read_sdcs.h"

#include "likelihood_nrqcd.h"
#include "fitter.h"

GlobalLikelihoodNRQCD get_likelihood(std::string dataDir, std::string sdcDir) {
  const auto data = read_all_data(dataDir, 2); // pT/M >= 2
  // TODO: read sdcs and construct likelihood proper with it

  return GlobalLikelihoodNRQCD(std::move(data));
}

void global_fit_nrqcd() {
  auto likelihood = get_likelihood("./data", "");



  LikelihoodFitter fitter;
  fitter.Fit(likelihood);
}
