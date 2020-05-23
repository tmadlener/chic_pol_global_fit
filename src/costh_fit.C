#include "data_structures.h"
#include "fitter.h"
#include "costh_ratio_llh.h"

#include "misc_utils.h"

#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <string>


void costh_fit(const std::string& outfile="results/costh_scan.root")
{
  const auto costh_ratio_8_12 = readData<CosthRatioData>("data/CMS_chic_costh_ratio_pt_8_12.dat");
  const auto costh_ratio_12_18 = readData<CosthRatioData>("data/CMS_chic_costh_ratio_pt_12_18.dat");
  const auto costh_ratio_18_30 = readData<CosthRatioData>("data/CMS_chic_costh_ratio_pt_18_30.dat");

  CosthRatioLikelihood likelihood({costh_ratio_8_12, costh_ratio_12_18, costh_ratio_18_30});

  LikelihoodFitter fitter;
  fitter.Fit(likelihood);

  TFile* scanFile = new TFile(outfile.c_str(), "recreate");
  TTree* scanTree = new TTree("log_like_scan", "log likelihood scan values");

  fitter.RandomScan(likelihood, scanTree);

  TTree* resultTree = new TTree("fit_result", "fit result information");
  fitter.storeFitResult(resultTree);
  resultTree->Write();

  fitter.setupScanTree(scanTree, false);
  scanTree->Write("", TObject::kWriteDelete);
  scanFile->Close();
}
