#include "data_structures.h"
#include "likelihood.h"
#include "fitter.h"

#include "TFile.h"
#include "TTree.h"

#include <string>

void global_fit(unsigned nScanPoints=5000000,
                const std::string scanFileName="results/scan_file.root")
{
  // cross section data
  const auto psi2S_ATLAS_cs = readData<CrossSectionData>("data/ATLAS_psi2S_cross_section.dat");
  const auto psi2S_CMS_cs = readData<CrossSectionData>("data/CMS_psi2S_cross_section.dat");
  const auto jpsi_CMS_cs = readData<CrossSectionData>("data/CMS_jpsi_cross_section.dat");
  const auto chic2_ATLAS_cs = readData<CrossSectionData>("data/ATLAS_chic2_cross_section.dat");
  const auto chic1_ATLAS_cs = readData<CrossSectionData>("data/ATLAS_chic1_cross_section.dat");
  const auto chic_ratio_CMS_cs = readData<CrossSectionData>("data/CMS_chic_ratio_cross_section.dat");

  // // polarization data
  const auto psi2S_CMS_pol = readData<PolarizationData>("data/CMS_psi2S_polarization.dat");
  const auto jpsi_CMS_pol = readData<PolarizationData>("data/CMS_jpsi_polarization.dat");



  GlobalLikelihood likelihood(psi2S_ATLAS_cs, psi2S_CMS_cs, chic2_ATLAS_cs, chic1_ATLAS_cs,
                              jpsi_CMS_cs, chic_ratio_CMS_cs, psi2S_CMS_pol, jpsi_CMS_pol);

  // fix all nuissances
  for (const auto* par : {"L_CMS", "L_ATLAS", "br_psip_dp", "br_psip_mm",
        "br_psip_c2", "br_psip_c1", "br_psip_jpsi", "br_c2_jpsi", "br_c1_jpsi", "br_jpsi_mm"}) {
    likelihood.fixParameter(par, 1);
  }


  LikelihoodFitter fitter;

  fitter.Fit(likelihood);

  TFile* scanFile = new TFile(scanFileName.c_str(), "recreate");
  TTree* scanTree = new TTree("log_like_scan", "log likelihood scan values");

  const ScanSettings scanParameters = {{0, 1, 101, "f_long_c2"},
                                       {0, 1, 101, "f_long_c1"}};

  fitter.Scan(likelihood, scanParameters, scanTree, nScanPoints);

  scanTree->Write("", TObject::kWriteDelete);
  scanFile->Close();
}
