#ifndef H_FITTER__
#define H_FITTER__

#include "data_structures.h"
#include "misc_utils.h"
#include "progress.h"

#include "Math/Minimizer.h"
#include "Fit/Fitter.h"
#include "Math/MinimizerOptions.h"
#include "TTree.h"
#include "TRandom3.h"

#include <iostream>

class LikelihoodFitter {
public:
  LikelihoodFitter() {
    ROOT::Math::MinimizerOptions minOpt;
    minOpt.SetMinimizerType("Minuit2");
    minOpt.SetPrintLevel(2);
    minOpt.SetDefaultErrorDef(0.5); // likelihood fit
    fitter.Config().SetMinimizerOptions(minOpt);
    // fitter.Config().SetMinosErrors(true);
  }


  /**
   * Run the fit
   */
  template<typename LLH>
  bool Fit(const LLH& llh)
  {
    const auto ret = FitFromParams(llh, llh.getStartParams());
    PrintResults();
    return ret;
  }


  template<typename LLH>
  void Scan(const LLH& llh, const ScanSettings& scanSettings, TTree* tree, const size_t nSamples=1000);

  void PrintResults() { fitter.GetMinimizer()->PrintResults(); }

private:

  template<typename LLH>
  bool FitFromParams(const LLH& llh, const ParamsSettings& params)
  {
    fitter.Config().SetParamsSettings(params);

    const int nPar = llh.nPars();     // default arguments are the ones
                                      // we want in this case
    return fitter.FitFCN(nPar, llh);
  }

  ROOT::Fit::Fitter fitter{};
};


template<typename LLH>
void LikelihoodFitter::Scan(const LLH& llh, const ScanSettings& scanSettings, TTree* tree, const size_t nSamples)
{
  // avoid too much output
  const int oldPrintLevel = fitter.Config().MinimizerOptions().PrintLevel();
  fitter.Config().MinimizerOptions().SetPrintLevel(0);

  auto& params = fitter.Config().ParamsSettings();

  std::vector<double> values(params.size() + 1);

  for (size_t i=0; i < params.size(); ++i) {
    tree->Branch(params[i].Name().c_str(), &values[i]);
  }
  tree->Branch("llh", &values[params.size()]);

  // first do a fit with the usual settings just to make sure that the scan is done around the minimum
  if (!FitFromParams(llh, params)) {
    std::cerr << "Could not find a valid minimum to scan around\n";
    return;
  }
  values.assign(fitter.Result().Parameters().cbegin(), fitter.Result().Parameters().cend());
  const double minimum = fitter.Result().MinFcnValue();
  values[params.size()] = minimum;
  tree->Fill();

  // get the internal indices of the parameters that are varied randomly now to
  // avoid having to get them over and over again
  std::vector<int> freeParams;
  for (size_t i=0; i < llh.nPars(); ++i) {freeParams.push_back(i);}
  std::pair<int, int> parIdcs;

  parIdcs.first = llh.getParIdx(scanSettings.first.name);
  freeParams.erase(std::find(freeParams.cbegin(), freeParams.cend(), parIdcs.first));
  parIdcs.second = llh.getParIdx(scanSettings.second.name);
  freeParams.erase(std::find(freeParams.cbegin(), freeParams.cend(), parIdcs.second));


  delete gRandom;
  gRandom = new TRandom3(0);

  // mean and sigma of normal distribution for sampling for all the parameters
  const auto& par = fitter.Result().Parameters();
  const auto& errs = fitter.Result().Errors();

  auto startTime = ProgressClock::now();
  size_t count{};
  const size_t total = scanSettings.first.n * scanSettings.second.n * nSamples;

  for (const double p1 : linspace(scanSettings.first.min, scanSettings.first.max, scanSettings.first.n)) {
    values[parIdcs.first] = p1;
    for (const double p2 : linspace(scanSettings.second.min, scanSettings.second.max, scanSettings.second.n)) {
      values[parIdcs.second] = p2;


      for (size_t iSamp=0; iSamp < nSamples; ++iSamp) {
        for (const int iPar : freeParams) {
          values[iPar] = std::abs(gRandom->Gaus(par[iPar], errs[iPar]));
        }

        const double val = llh(values.data());
        // Only storing "reasonable" candidate events. Using 2 * nPars should ensure
        // that even the 99 % CLs can be drawn for fits with up to a few hundred
        // parameters
        if (val < minimum + params.size() * 2) {
          values[params.size()] = val;
          tree->Fill();
        }
        count++;
        printProgress(count, total - 1, startTime, total / nSamples / 5);
      }

    }
  }

  fitter.Config().MinimizerOptions().SetPrintLevel(oldPrintLevel);
}

#endif
