#include "data_reading.h"
#include "read_sdcs.h"

#include "fitter.h"
#include "likelihood_nrqcd.h"

GlobalLikelihoodNRQCD get_likelihood(std::string dataDir, std::string sdcDir) {
  const auto data = read_all_data(dataDir, 2); // pT/M >= 2
  // TODO: read sdcs and construct likelihood proper with it
  const auto psi_SDC = readPsiSDCs(sdcDir);
  const auto chic1_SDC = readChic1SDCs(sdcDir);
  const auto chic2_SDC = readChic2SDCs(sdcDir);

  return GlobalLikelihoodNRQCD(std::move(data), std::move(psi_SDC), std::move(chic1_SDC), std::move(chic2_SDC));
}

void global_fit_nrqcd() {
  auto likelihood = get_likelihood("./data", "../SDCs_Chung");

  LikelihoodFitter fitter;
  fitter.Fit(likelihood);
}
