#ifndef DATAREADING_H
#define DATAREADING_H

#include "data_structures.h"

#include <string>
#include <vector>

std::vector<PtCosthRatioMeasurement> read_costh_ratios(const std::string& configfile,
                                                       const std::string& datadir = "data") {
  std::ifstream cfile(datadir + "/" + configfile);
  std::string line;

  if (!cfile) { std::cerr << "Could not open file: " << configfile << std::endl; }

  std::vector<PtCosthRatioMeasurement> data;

  while (std::getline(cfile, line)) {
    std::stringstream sline(line);
    double ptM, low, high;
    std::string ratiofile;
    sline >> ptM >> low >> high >> ratiofile;

    ptM /= M_JPSI; // convert pT from configfile to pt/M
    low /= M_JPSI;
    high /= M_JPSI;

    BinInfo point{ptM, low, high};

    ratiofile = datadir + "/" + ratiofile;
    data.emplace_back(point, readData<CosthRatioData>(ratiofile));
  }

  return data;
}

template <typename LLH>
LLH get_likelihood(bool useCosthRatios, bool usePsiPol = true, bool clipCorrs = false,
                   const std::string& datadir = "./data/", const double minPTM = 2) {
  // cross section data
  const auto psi2S_ATLAS_cs = readData<CrossSectionData>(datadir + "/ATLAS_psi2S_cross_section.dat", minPTM);
  const auto psi2S_CMS_cs = readData<CrossSectionData>(datadir + "/CMS_psi2S_cross_section.dat", minPTM);
  const auto jpsi_CMS_cs = readData<CrossSectionData>(datadir + "/CMS_jpsi_cross_section.dat", minPTM);
  const auto chic2_ATLAS_cs = readData<CrossSectionData>(datadir + "/ATLAS_chic2_cross_section.dat", minPTM);
  const auto chic1_ATLAS_cs = readData<CrossSectionData>(datadir + "/ATLAS_chic1_cross_section.dat", minPTM);
  const auto chic_ratio_CMS_cs = readData<CrossSectionData>(datadir + "/CMS_chic_ratio_cross_section.dat", minPTM);

  // // polarization data
  const auto psi2S_CMS_pol = readData<PolarizationData>(datadir + "/CMS_psi2S_polarization.dat", minPTM);
  const auto jpsi_CMS_pol = readData<PolarizationData>(datadir + "/CMS_jpsi_polarization.dat", minPTM);

  if (!useCosthRatios) {
    return LLH(psi2S_ATLAS_cs, psi2S_CMS_cs, chic2_ATLAS_cs, chic1_ATLAS_cs, jpsi_CMS_cs, chic_ratio_CMS_cs,
               psi2S_CMS_pol, jpsi_CMS_pol, clipCorrs);
  }

  const auto chic_ratios_CMS_pol = read_costh_ratios("chic_ratios.conf", datadir);

  return LLH(psi2S_ATLAS_cs, psi2S_CMS_cs, chic2_ATLAS_cs, chic1_ATLAS_cs, jpsi_CMS_cs, chic_ratio_CMS_cs,
             psi2S_CMS_pol, jpsi_CMS_pol, chic_ratios_CMS_pol, usePsiPol, clipCorrs);
}

#endif
