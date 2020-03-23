#ifndef H_FITTER__
#define H_FITTER__

#include "data_structures.h"

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
  void Scan(const LLH& llh, const ScanSettings& scanSettings, TTree* tree, const size_t nPoints=1000);

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
void LikelihoodFitter::Scan(const LLH& llh, const ScanSettings& scanSettings, TTree* tree, const size_t nPoints)
{
  // avoid too much output
  const int oldPrintLevel = fitter.Config().MinimizerOptions().PrintLevel();
  fitter.Config().MinimizerOptions().SetPrintLevel(0);

  auto& params = fitter.Config().ParamsSettings();

  std::vector<double> values;
  values.reserve(params.size() + 1);

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

  // get the internal indices now to avoid having to get them over and over again
  std::vector<int> freeParams;
  for (size_t i=0; i < llh.nPars(); ++i) {freeParams.push_back(i);}

  std::vector<int> parIdcs;
  for (const auto & s : scanSettings) {
    const int idx = llh.getParIdx(s.name);
    parIdcs.push_back(idx);
    freeParams.erase(std::find(freeParams.cbegin(), freeParams.cend(), idx));
  }


  delete gRandom;
  gRandom = new TRandom3(0);

  const auto& par = fitter.Result().Parameters();
  const auto& errs = fitter.Result().Errors();

  for (size_t iScan=0; iScan < nPoints; ++iScan) {

    for (size_t iPar=0; iPar < scanSettings.size(); ++iPar) {
      values[parIdcs[iPar]] = gRandom->Uniform(scanSettings[iPar].min, scanSettings[iPar].max);
    }

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
  }

  fitter.Config().MinimizerOptions().SetPrintLevel(oldPrintLevel);
}

#endif
