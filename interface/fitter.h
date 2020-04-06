#ifndef H_FITTER__
#define H_FITTER__

#include "data_structures.h"
#include "multivariate_normal_distribution.h"
#include "progress.h"
#include "matrix_helper.h"

#include "Math/Minimizer.h"
#include "Fit/Fitter.h"
#include "Math/MinimizerOptions.h"
#include "TTree.h"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <numeric>

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

  template<typename LLH>
  void RandomScan(const LLH& llh, TTree* tree, const size_t nSamples=1000000,
                  const std::vector<std::pair<int, double>>& covRedFactors={},
                  const double maxDeltaLLH=std::numeric_limits<double>::max());

  void PrintResults() { fitter.GetMinimizer()->PrintResults(); }

  const ROOT::Fit::FitResult& Result() const { return fitter.Result(); }

  std::vector<double> GetCovMatrix() const {
    if (fitter.Result().CovMatrixStatus() != 3) {
      std::cerr << "WARNING: cov matrix status != 3\n";
      return {};
    }

    // fill the covariance matrix fully, instead of only elements below the
    // diagonal. This is done because than the transformation to an
    // Eigen::Matrix can be done very easily.
    const size_t nPars = fitter.Result().Parameters().size();
    std::vector<double> covMatrix(nPars * nPars);
    for (size_t i = 0; i < nPars; ++i) {
      for (size_t j = 0; j < nPars; ++j) {
        covMatrix[i * nPars + j] = fitter.Result().CovMatrix(i, j);
      }
    }

    return covMatrix;
  }

  /**
   * store the parameter values, the covariance matrix (as 1D vector) and the
   * likelihood value at the minimum
   */
  void storeFitResult(TTree* tree)
  {
    auto* params = new std::vector<double>(fitter.Result().Parameters());
    tree->Branch("parameters", &params);

    auto* covMatrix = new std::vector<double>(GetCovMatrix());
    tree->Branch("cov_matrix", &covMatrix);

    double minVal = fitter.Result().MinFcnValue();
    tree->Branch("min_llh", &minVal);

    tree->Fill();

    delete params;
    delete covMatrix;
  }

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

  // get a COPY here, in order to always start the fit at some "reasonable"
  // starting point
  auto params = fitter.Config().ParamsSettings();

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

  // get the internal indices of the parameters that are fixed
  const std::pair<int, int> parIdcs = {
    llh.getParIdx(scanSettings.first.name),
    llh.getParIdx(scanSettings.second.name)
  };

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
      printProgress<PrintStyle::ProgressText>(count, total, startTime, 100);
    }
  }

  fitter.Config().MinimizerOptions().SetPrintLevel(oldPrintLevel);
}


template<typename LLH>
void LikelihoodFitter::RandomScan(const LLH& llh, TTree* tree, const size_t nSamples,
                                  const std::vector<std::pair<int, double>>& covRedFactors,
                                  const double maxDeltaLLH)
{
  // avoid too much output from the fitter
  const int oldPrintLevel = fitter.Config().MinimizerOptions().PrintLevel();
  fitter.Config().MinimizerOptions().SetPrintLevel(0);

  const auto& params = fitter.Config().ParamsSettings();

  // first do a fit with the usual settings just to make sure that the scan is done around the minimum
  if (!FitFromParams(llh, params)) {
    std::cerr << "Could not find a valid minimum to scan around\n";
    return;
  }
  fitter.Config().MinimizerOptions().SetPrintLevel(oldPrintLevel);

  std::vector<double> values(params.size());
  for (size_t i=0; i < params.size(); ++i) {
    tree->Branch(params[i].Name().c_str(), &values[i]);
  }
  double llh_val;
  tree->Branch("llh", &llh_val);

  // fill the minimum in any case
  values.assign(fitter.Result().Parameters().cbegin(), fitter.Result().Parameters().cend());
  llh_val = fitter.Result().MinFcnValue();
  tree->Fill();

  const double min_llh = llh_val;

  std::vector<double> covMatrix = GetCovMatrix();
  if (covMatrix.empty()) {
    return;
  }

  covMatrix = change_variance(covMatrix, covRedFactors);

  auto means = fitter.Result().Parameters();
  const MultivariateNormalDistribution multiVarNorm(means, covMatrix);
  size_t stored = 0;

  const auto startTime = ProgressClock::now();
  for (size_t i = 0; i < nSamples; ++i) {
    const auto pars = multiVarNorm();
    llh_val = llh(pars.data());

    if (!std::isnan(llh_val) && (llh_val - min_llh) < maxDeltaLLH) {
      values.assign(pars.cbegin(), pars.cend());
      tree->Fill();
      stored++;
    }

    printProgress(i, nSamples - 1, startTime, 100);
  }
  std::cout << "Generated " << nSamples << " random sampling points and stored "
            << stored << " evaluations of the likelihood\n";
}

#endif
