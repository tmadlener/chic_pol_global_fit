#include "data_structures.h"
#include "likelihood.h"
#include "fitter.h"
#include "misc_utils.h"

#include "TFile.h"
#include "TTree.h"

#include <string>

void global_fit(const std::string scanFileName="results/scan_file.root",
                unsigned nScanPoints1=26,
                unsigned nScanPoints2=26,
                const double fLow1=0, const double fHigh1=1,
                const double fLow2=0, const double fHigh2=1)
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

  const ScanSettings scanParameters = {{linspace(fLow2, fHigh2, nScanPoints2), "f_long_c2"},
                                       {linspace(fLow1, fHigh1, nScanPoints1), "f_long_c1"}};

  fitter.Scan(likelihood, scanParameters, scanTree);

  scanTree->Write("", TObject::kWriteDelete);
  scanFile->Close();
}
