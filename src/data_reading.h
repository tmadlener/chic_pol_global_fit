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

GlobalFitData read_all_data(std::string datadir = "./data/", const double minPTM = 2) {
  GlobalFitData data;
  // cross section data
  data.psi2S_ATLAS_cs = readData<CrossSectionData>(datadir + "/ATLAS_psi2S_cross_section.dat", minPTM);
  data.psi2S_CMS_cs = readData<CrossSectionData>(datadir + "/CMS_psi2S_cross_section.dat", minPTM);
  data.jpsi_CMS_cs = readData<CrossSectionData>(datadir + "/CMS_jpsi_cross_section.dat", minPTM);
  data.chic2_ATLAS_cs = readData<CrossSectionData>(datadir + "/ATLAS_chic2_cross_section.dat", minPTM);
  data.chic1_ATLAS_cs = readData<CrossSectionData>(datadir + "/ATLAS_chic1_cross_section.dat", minPTM);
  data.chic_ratio_CMS_cs = readData<CrossSectionData>(datadir + "/CMS_chic_ratio_cross_section.dat", minPTM);

  // polarization data
  data.psi2S_CMS_pol = readData<PolarizationData>(datadir + "/CMS_psi2S_polarization.dat", minPTM);
  data.jpsi_CMS_pol = readData<PolarizationData>(datadir + "/CMS_jpsi_polarization.dat", minPTM);

  // polarization data
  data.psi2S_CMS_pol = readData<PolarizationData>(datadir + "/CMS_psi2S_polarization.dat", minPTM);
  data.jpsi_CMS_pol = readData<PolarizationData>(datadir + "/CMS_jpsi_polarization.dat", minPTM);

  data.chic_costh_ratios_CMS = read_costh_ratios("chic_ratios.conf", datadir);

  return data;
}

#endif
