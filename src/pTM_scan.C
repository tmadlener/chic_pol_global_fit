#if NRQCD_FIT
#include "likelihood_nrqcd.h"
#include "read_sdcs.h"
#else
#include "likelihood.h"
#endif

#include "multivariate_normal_distribution.h"
#include "store_variables.h"


#include "progress.h"
#include "misc_utils.h"

#include "TFile.h"
#include "TTree.h"

#include <vector>
#include <utility>
#include <string>
#include <iostream>


std::pair<std::vector<double>, std::vector<double>> readFitResults(const std::string& filename)
{
  TFile* file = TFile::Open(filename.c_str());
  TTree* tree = static_cast<TTree*>(file->Get("fit_result"));
  if (!tree) {
    std::cerr << "Cannot find \'fit_result\' tree in file: " << filename << '\n';
  }

  std::vector<double>* par = new std::vector<double>();
  std::vector<double>* cov = new std::vector<double>();

  tree->SetBranchAddress("parameters", &par);
  tree->SetBranchAddress("cov_matrix", &cov);

  tree->GetEntry(0);

  file->Close();

  return {*par, *cov};
}

template<typename LLH>
size_t scan_ptM_point(const double ptM, StoreVariables<LLH>& storeVars,
                      const MultivariateNormalDistribution<>& mvn, TTree* tree,
                      size_t nscans=1000)
{
  size_t nValid = 0;
  for (size_t iScan = 0; iScan < nscans; ++iScan) {
    if (storeVars(ptM, mvn())) {
      tree->Fill();
      nValid++;
    }
  }
  return nValid;
}



void pTM_scan(const std::string& fitResultsFile,
              const std::string& scanFileName,
              const double ptMmin, const double ptMmax, const size_t nPoints, const size_t nScans,
#if NRQCD_FIT
              const std::string& sdcDir,
#if DO_COMBINATION_FIT
              const bool storeParams=true, const double lpFraction=0.)
#else
              const bool storeParams=true, sdc::SDCType order=sdc::SDCType::NLO)
#endif // DO_COMBINATION_FIT
#else
              const bool storeParams=true)
#endif

{
  std::vector<double> params;
  std::vector<double> cov;

  std::tie(params, cov) = readFitResults(fitResultsFile);
  const MultivariateNormalDistribution<> multiVarNorm(params, cov);

#if NRQCD_FIT
#if DO_COMBINATION_FIT
  GlobalLikelihoodNRQCD likelihood(readPsiSDCs(sdcDir, lpFraction),
                                   readChic1SDCs(sdcDir, lpFraction),
                                   readChic2SDCs(sdcDir, lpFraction));
#else
  GlobalLikelihoodNRQCD likelihood(readPsiSDCs(sdcDir, order),
                                   readChic1SDCs(sdcDir, order),
                                   readChic2SDCs(sdcDir, order));
#endif // DO_COMBINATION_FIT
  StoreVariables<GlobalLikelihoodNRQCD> storeVars(likelihood);
#else
  GlobalLikelihood likelihood;
  StoreVariables<GlobalLikelihood> storeVars(likelihood);
#endif
  TFile* scanFile = new TFile(scanFileName.c_str(), "recreate");
  TTree* scanTree = storeVars.create(storeParams);

  const size_t nEvaluations = nPoints * nScans;
  size_t valid = 0;
  size_t nProcessed = 0;

  const auto ptMScanPoints = nPoints > 1 ? linspace(ptMmin, ptMmax, nPoints) : std::vector<double>{ptMmin};
  std::cout << "Starting scan of " << nPoints << " ptM points between " << ptMmin << " and " << ptMmax
            << " doing " << nScans << " random evaluations at each point\n";
  const auto startTime = ProgressClock::now();

  for (const double ptM : ptMScanPoints) {
    valid += scan_ptM_point(ptM, storeVars, multiVarNorm, scanTree, nScans);
    nProcessed += nScans;
    printProgress(nProcessed, nEvaluations, startTime, nEvaluations);
  }

  std::cout << "Stored " << valid << " / " << nEvaluations << " evaluation points ("
            << 100.0 * valid / nEvaluations << " % valid)"<< "\n";

  scanTree->Write("", TObject::kWriteDelete);
  scanFile->Close();
}
