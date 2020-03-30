#ifndef H_FITTER__
#define H_FITTER__

#include "data_structures.h"
#include "progress.h"

#include "Math/Minimizer.h"
#include "Fit/Fitter.h"
#include "Math/MinimizerOptions.h"
#include "TTree.h"
#include "TRandom3.h"

#include <iostream>
#include <algorithm>
#include <cmath>

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
  void Scan(const LLH& llh, const ScanSettings& scanSettings, TTree* tree);

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
void LikelihoodFitter::Scan(const LLH& llh, const ScanSettings& scanSettings, TTree* tree)
{
  // avoid too much output from the fitter
  const int oldPrintLevel = fitter.Config().MinimizerOptions().PrintLevel();
  fitter.Config().MinimizerOptions().SetPrintLevel(0);

  auto& params = fitter.Config().ParamsSettings();

  // first do a fit with the usual settings just to make sure that the scan is done around the minimum
  if (!FitFromParams(llh, params)) {
    std::cerr << "Could not find a valid minimum to scan around\n";
    return;
  }

  std::vector<double> values(params.size());

  for (size_t i=0; i < params.size(); ++i) {
    tree->Branch(params[i].Name().c_str(), &values[i]);
  }
  double llh_val;
  tree->Branch("llh", &llh_val);


  // store all the obtained values, but make it possible to filter out
  // non-successful fits
  int goodFit = 1;
  tree->Branch("goodFit", &goodFit);

  values.assign(fitter.Result().Parameters().cbegin(), fitter.Result().Parameters().cend());
  llh_val = fitter.Result().MinFcnValue();
  tree->Fill();

  // get the internal indices of the parameters that are not fixed fro scanning
  std::vector<int> freeParams;
  for (size_t i=0; i < llh.nPars(); ++i) {freeParams.push_back(i);}
  std::pair<int, int> parIdcs;

  parIdcs.first = llh.getParIdx(scanSettings.first.name);
  freeParams.erase(std::find(freeParams.cbegin(), freeParams.cend(), parIdcs.first));
  parIdcs.second = llh.getParIdx(scanSettings.second.name);
  freeParams.erase(std::find(freeParams.cbegin(), freeParams.cend(), parIdcs.second));


  delete gRandom;
  gRandom = new TRandom3(0);

  auto startTime = ProgressClock::now();
  size_t count{};
  const size_t total = scanSettings.first.values.size() * scanSettings.second.values.size();

  // For each parameter set, find the minimal chi2 value with fixed values. If
  // this value is above the chi2 level corresponding to a given CL, then no
  // configuration of the other parameters can exist, for which a smaller chi2
  // value can be achieved (given the fixed values of the scan parameters).
  // Hence, this point is definitely outside of the projected contour. All other
  // points must be within, since such a configuration must exist, but since only
  // the 2d projection is of interest the exact configuration is of no concern.
  for (const double p1 : scanSettings.first.values) {
    params[parIdcs.first].Fix();
    params[parIdcs.first].SetValue(p1);

    for (const double p2 : scanSettings.second.values) {
      params[parIdcs.second].Fix();
      params[parIdcs.second].SetValue(p2);

      goodFit = 0;

      const auto& parResults = fitter.Result().Parameters();
      if (FitFromParams(llh, params) &&
          std::none_of(parResults.cbegin(), parResults.cend(), [] (double d) { return std::isnan(d); })) {
        goodFit = 1;
      }

      values.assign(parResults.cbegin(), parResults.cend());
      llh_val = fitter.Result().MinFcnValue();
      tree->Fill();

      count++;
      printProgress<PrintStyle::ProgressText>(count, total - 1, startTime, 100);
    }
  }

  fitter.Config().MinimizerOptions().SetPrintLevel(oldPrintLevel);
}

#endif
