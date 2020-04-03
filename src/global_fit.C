#include "data_structures.h"
#include "likelihood.h"
#include "fitter.h"
#include "misc_utils.h"

#include "TFile.h"
#include "TTree.h"

#include <fstream>
#include <string>
#include <iostream>


std::vector<double> to_frac(const std::vector<double>& lambdas)
{
  std::vector<double> fracs;
  fracs.reserve(lambdas.size());
  for (const double l: lambdas) {
    fracs.push_back((1 - l) / (3 + l));
  }
  return fracs;
}

std::vector<PtCosthRatioMeasurement> read_costh_ratios(const std::string& configfile)
{
  std::ifstream cfile(configfile);
  std::string line;

  if (!cfile) {
    std::cerr << "Could not open file: " << configfile << std::endl;
  }

  std::vector<PtCosthRatioMeasurement> data;

  while(std::getline(cfile, line)) {
    std::stringstream sline(line);
    double ptM, low, high;
    std::string ratiofile;
    sline >> ptM >> low >> high >> ratiofile;

    ptM /= M_JPSI; // convert pT from configfile to pt/M
    low /= M_JPSI;
    high /= M_JPSI;

    BinInfo point{ptM, low, high};

    data.emplace_back(point, readData<CosthRatioData>(ratiofile));
  }

  return data;
}

GlobalLikelihood get_likelihood(bool useCosthRatios) {
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

  if (!useCosthRatios) {
    return GlobalLikelihood(psi2S_ATLAS_cs, psi2S_CMS_cs, chic2_ATLAS_cs, chic1_ATLAS_cs,
                            jpsi_CMS_cs, chic_ratio_CMS_cs, psi2S_CMS_pol, jpsi_CMS_pol);
  }

  const auto chic_ratios_CMS_pol = read_costh_ratios("data/chic_ratios.conf");

  return GlobalLikelihood(psi2S_ATLAS_cs, psi2S_CMS_cs, chic2_ATLAS_cs, chic1_ATLAS_cs,
                          jpsi_CMS_cs, chic_ratio_CMS_cs, psi2S_CMS_pol, jpsi_CMS_pol,
                          chic_ratios_CMS_pol);


}

void global_fit(const std::string& scanFileName="results/scan_file.root",
                const std::string& graphFileName="results/fit_result_graphs.root",
                unsigned nSamplePoints=1000000, const bool useCosthRatios=true,
                const bool storeGraphs=true)
{
  GlobalLikelihood likelihood = get_likelihood(useCosthRatios);

  LikelihoodFitter fitter;

  fitter.Fit(likelihood);

  TFile* scanFile = new TFile(scanFileName.c_str(), "recreate");
  TTree* scanTree = new TTree("log_like_scan", "log likelihood scan values");


  // Using maxDeltaLLH = 50, allows to go to a deltaChi2 value of 100, which
  // should be enough for the most important purposes
  fitter.RandomScan(likelihood, scanTree, nSamplePoints, 50);

  scanTree->Write("", TObject::kWriteDelete);
  scanFile->Close();

  if (storeGraphs) {
    TFile* graphFile = new TFile(graphFileName.c_str(), "recreate");
    for (const auto& graph : likelihood.getDataGraphs(fitter.Result())) {
      graph.Write();
    }

    graphFile->Close();
  }

}
