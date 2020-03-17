#include "data_structures.h"
#include "likelihood.h"
#include "fitter.h"

void global_fit()
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

  LikelihoodFitter fitter;

  fitter.Fit(likelihood);
}
