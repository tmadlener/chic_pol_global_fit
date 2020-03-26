#include "data_structures.h"
#include "fitter.h"
#include "costh_ratio_llh.h"

#include "misc_utils.h"

#include "TFile.h"
#include "TTree.h"

#include <iostream>

void costh_fit()
{
  const auto costh_ratio_8_12 = readData<CosthRatioData>("data/CMS_chic_costh_ratio_pt_8_12.dat");
  const auto costh_ratio_12_18 = readData<CosthRatioData>("data/CMS_chic_costh_ratio_pt_12_18.dat");
  const auto costh_ratio_18_30 = readData<CosthRatioData>("data/CMS_chic_costh_ratio_pt_18_30.dat");

  CosthRatioLikelihood likelihood({costh_ratio_8_12, costh_ratio_12_18, costh_ratio_18_30});

  LikelihoodFitter fitter;
  fitter.Fit(likelihood);

  TFile* scanFile = new TFile("results/costh_scan.root", "recreate");
  TTree* scanTree = new TTree("log_like_scan", "log likelihood scan values");

  const ScanSettings scanParameters = {{linspace(-2., 2., 100), "lambda_1"},
                                       {linspace(-2., 2., 100), "lambda_2"}};

  fitter.Scan(likelihood,
              scanParameters,
              scanTree);

  scanTree->Write("", TObject::kWriteDelete);
  scanFile->Close();
}
