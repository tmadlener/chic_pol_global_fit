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

size_t scan_ptM_point(const double ptM, StoreVariables& storeVars,
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
              const bool storeParams=true)
{
  std::vector<double> params;
  std::vector<double> cov;

  std::tie(params, cov) = readFitResults(fitResultsFile);
  const MultivariateNormalDistribution<> multiVarNorm(params, cov);

  StoreVariables storeVars;
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
