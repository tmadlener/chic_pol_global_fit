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

bool isGoodFit(int goodFit, const std::vector<double>& parValues) {
  return (goodFit &&
          std::none_of(parValues.cbegin(), parValues.cend(), [] (double d) { return std::isnan(d); }));
}

struct MaxDeltaLLH {
  bool operator()(double llh, double min_llh) const { return (llh - min_llh) < minDelta; }

  double minDelta{std::numeric_limits<double>::max()};
};

struct DeltaChi2 {
  DeltaChi2(double dChi2=2.2977) : high(dChi2 + 0.1), low(dChi2 - 0.1) {}

  bool operator()(double llh, double min_llh) const {
    const double delta = 2 * (llh - min_llh);
    return (delta > low && delta < high);
  }

  double high;
  double low;
};

class LikelihoodFitter {
public:
  LikelihoodFitter(bool withMinos=false) {
    ROOT::Math::MinimizerOptions minOpt;
    minOpt.SetMinimizerType("Minuit2");
    minOpt.SetPrintLevel(2);
    minOpt.SetErrorDef(0.5); // likelihood fit
    fitter.Config().SetMinimizerOptions(minOpt);
    fitter.Config().SetMinosErrors(withMinos);
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

  template<typename LLH, typename StoreCondF=decltype(MaxDeltaLLH{})>
  void RandomScan(const LLH& llh, TTree* tree, const size_t nSamples=1000000,
                  const double varRedFactor=1.0,
                  StoreCondF storeCondF=MaxDeltaLLH{});

  void PrintResults() {
    fitter.GetMinimizer()->PrintResults();
    if (fitter.Config().MinosErrors()) {
      std::cout << "MINOS Errors:\n";
      const auto& fitRes = fitter.Result();
      const auto& parValues = fitRes.Parameters();
      for (size_t i = 0; i < fitRes.NTotalParameters(); ++i) {
        const auto name = fitRes.GetParameterName(i);
        std::cout << name << "\t = " << parValues[i] << "\t";
        if (fitRes.HasMinosError(i)) {
          std::cout << "+" << fitRes.UpperError(i) << "\t" << fitRes.LowerError(i) << "\n";
        }
      }
    }
  }

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

    std::vector<double>* lowErrors = nullptr;
    std::vector<double>* upErrors = nullptr;

    if (fitter.Config().MinosErrors()) {
      const auto& fitRes = fitter.Result();

      lowErrors = new std::vector<double>(fitRes.Parameters().size());
      upErrors = new std::vector<double>(fitRes.Parameters().size());

      for (size_t i = 0; i < fitRes.NTotalParameters(); ++i) {
        if (fitRes.HasMinosError(i)) {
          lowErrors->at(i) = fitRes.LowerError(i);
          upErrors->at(i) = fitRes.UpperError(i);
        }
      }

      tree->Branch("minos_error_low", &lowErrors);
      tree->Branch("minos_error_up", &upErrors);
    }

    tree->Fill();

    delete params;
    delete covMatrix;

    if (lowErrors) delete lowErrors;
    if (upErrors) delete upErrors;
  }

  /**
   * Set up the tree for scanning
   */
  void setupScanTree(TTree* tree, bool goodFitBranch=true)
  {
    const auto& params = fitter.Config().ParamsSettings();
    m_values.reserve(params.size());

    for (size_t i=0; i < params.size(); ++i) {
      tree->Branch(params[i].Name().c_str(), &m_values[i]);
    }
    tree->Branch("llh", &m_llh_val);

    if (goodFitBranch) {
      tree->Branch("goodFit", &m_goodFit);
    }

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

  // for scanning. These need to be class variables in order to be able to use
  // them as TTree branches
  double m_llh_val{0};
  int m_goodFit{0};
  std::vector<double> m_values;
};


template<typename LLH>
void LikelihoodFitter::Scan(const LLH& llh, const ScanSettings& scanSettings, TTree* tree)
{
  // avoid too much output from the fitter
  const int oldPrintLevel = fitter.Config().MinimizerOptions().PrintLevel();
  fitter.Config().MinimizerOptions().SetPrintLevel(0);
  // disable minos as we are really only interested in the minimum here, not the
  // uncertainties
  const bool oldMinos = fitter.Config().MinosErrors();
  fitter.Config().SetMinosErrors(false);

  // first do a fit with the usual settings just to make sure that the scan is done around the minimum
  if (!FitFromParams(llh, llh.getStartParams())) {
    std::cerr << "Could not find a valid minimum to scan around\n";
    return;
  }

  m_goodFit = 1;
  m_values.assign(fitter.Result().Parameters().cbegin(), fitter.Result().Parameters().cend());
  m_llh_val = fitter.Result().MinFcnValue();
  tree->Fill();

  // get the internal indices of the parameters that are fixed
  const std::pair<int, int> parIdcs = {
    llh.getParIdx(scanSettings.first.name),
    llh.getParIdx(scanSettings.second.name)
  };

  auto startTime = ProgressClock::now();
  size_t count{};
  const size_t total = scanSettings.first.values.size() * scanSettings.second.values.size();

  // get a COPY here, in order to always start the fit at some "reasonable"
  // starting point
  auto params = fitter.Config().ParamsSettings();

  // keep a list of parameter settings that have been used successfully in a
  // previous fit in case the fit from the "standard" starting point does not
  // work converge
  auto lastGoodParams = params;

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


      m_goodFit = FitFromParams(llh, params);
      const auto& parResults = fitter.Result().Parameters();

      if (isGoodFit(m_goodFit, parResults)) {
        m_values.assign(parResults.cbegin(), parResults.cend());
        m_llh_val = fitter.Result().MinFcnValue();

        lastGoodParams = fitter.Config().ParamsSettings();
      } else {
        // two more options at our disposal to get a fit
        // 1) retry from the last values that lead to a converging fit
        lastGoodParams[parIdcs.first].Fix();
        lastGoodParams[parIdcs.first].SetValue(p1);
        lastGoodParams[parIdcs.second].Fix();
        lastGoodParams[parIdcs.second].SetValue(p2);

        std::cout << "Retrying fitting from last good parameter values\n";

        m_goodFit = 2 * FitFromParams(llh, lastGoodParams); // indicate a refit in the variable
        const auto& parResults2 = fitter.Result().Parameters();

        if (isGoodFit(m_goodFit, parResults2)) {
          m_values.assign(parResults2.cbegin(), parResults2.cend());
          m_llh_val = fitter.Result().MinFcnValue();

          lastGoodParams = fitter.Config().ParamsSettings();
        } else {
          // 2) start from the very starting parameters again
          auto startParams = llh.getStartParams();
          startParams[parIdcs.first].Fix();
          startParams[parIdcs.first].SetValue(p1);
          startParams[parIdcs.second].Fix();
          startParams[parIdcs.second].SetValue(p2);

          std::cout << "Retrying fitting from the starting parameter values\n";

          m_goodFit = 3 * FitFromParams(llh, startParams);
          const auto& parResults3 = fitter.Result().Parameters();

          if (isGoodFit(m_goodFit, parResults3)) {
            m_values.assign(parResults3.cbegin(), parResults3.cend());
            m_llh_val = fitter.Result().MinFcnValue();
          } else {
            std::cout << "None of the three trials of fitting resulted in a converging fit for "
                      << scanSettings.first.name << " = " << p1 << ", "
                      << scanSettings.second.name << " = " << p2 << "\n";
            m_goodFit = 0;
            m_llh_val = llh(parResults3.data());
          }
        }
      }

      tree->Fill();

      count++;
      printProgress<PrintStyle::ProgressText>(count, total, startTime, 100);
    }
  }

  fitter.Config().MinimizerOptions().SetPrintLevel(oldPrintLevel);
  fitter.Config().SetMinosErrors(oldMinos);
}


template<typename LLH, typename StoreCondF>
void LikelihoodFitter::RandomScan(const LLH& llh, TTree* tree, const size_t nSamples,
                                  const double varRedFactor,
                                  StoreCondF storeCond)
{
  // avoid too much output from the fitter
  const int oldPrintLevel = fitter.Config().MinimizerOptions().PrintLevel();
  fitter.Config().MinimizerOptions().SetPrintLevel(0);

  // first do a fit with the usual settings just to make sure that the scan is done around the minimum
  if (!FitFromParams(llh, fitter.Config().ParamsSettings())) {
    std::cerr << "Could not find a valid minimum to scan around\n";
    return;
  }
  fitter.Config().MinimizerOptions().SetPrintLevel(oldPrintLevel);

  // fill the minimum in any case
  m_values.assign(fitter.Result().Parameters().cbegin(), fitter.Result().Parameters().cend());
  m_llh_val = fitter.Result().MinFcnValue();
  tree->Fill();

  const double min_llh = m_llh_val;

  std::vector<double> covMatrix = GetCovMatrix();
  if (covMatrix.empty()) {
    return;
  }

  covMatrix = reduce_variance(covMatrix, varRedFactor);

  auto means = fitter.Result().Parameters();
  const MultivariateNormalDistribution<> multiVarNorm(means, covMatrix);
  size_t stored = 0;

  const auto startTime = ProgressClock::now();
  for (size_t i = 0; i < nSamples; ++i) {
    const auto pars = multiVarNorm();
    m_llh_val = llh(pars.data());

    if (!std::isnan(m_llh_val) && storeCond(m_llh_val, min_llh)) {
      m_values.assign(pars.cbegin(), pars.cend());
      tree->Fill();
      stored++;
    }

    printProgress(i, nSamples - 1, startTime, 100);
  }
  std::cout << "Generated " << nSamples << " random sampling points and stored "
            << stored << " evaluations of the likelihood\n";
}

#endif
