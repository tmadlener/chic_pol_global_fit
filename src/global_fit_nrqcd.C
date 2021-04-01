#include "data_reading.h"

#include "likelihood_nrqcd.h"
#include "fitter.h"

void global_fit_nrqcd() {
  auto likelihood = get_likelihood<GlobalLikelihoodNRQCD>(true, true, false, "./data");

  LikelihoodFitter fitter;
  fitter.Fit(likelihood);
}
