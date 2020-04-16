#include "data_structures.h"
#include "likelihood.h"
#include "fitter.h"
#include "misc_utils.h"

#include "TFile.h"
#include "TTree.h"

#include <fstream>
#include <string>
#include <iostream>
#include <array>

struct ScanArguments {
  bool randomScan{false};
  unsigned nSamplePoints{0};
  double flow1{-1};
  double fhigh1{1};
  double flow2{-1};
  double fhigh2{1};
  unsigned nScan1{51};
  unsigned nScan2{51};
};


std::vector<double> to_frac(const std::vector<double>& lambdas)
{
  std::vector<double> fracs;
  fracs.reserve(lambdas.size());
  for (const double l: lambdas) {
    fracs.push_back((1 - l) / (3 + l));
  }
  return fracs;
}

std::vector<PtCosthRatioMeasurement> read_costh_ratios(const std::string& configfile,
                                                       const std::string& datadir="data")
{
  std::ifstream cfile(datadir + "/" + configfile);
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

    ratiofile = datadir + "/" + ratiofile;
    data.emplace_back(point, readData<CosthRatioData>(ratiofile));
  }

  return data;
}

GlobalLikelihood get_likelihood(bool useCosthRatios,
                                const std::string& datadir="./data/") {
  // cross section data
  const auto psi2S_ATLAS_cs = readData<CrossSectionData>(datadir + "ATLAS_psi2S_cross_section.dat");
  const auto psi2S_CMS_cs = readData<CrossSectionData>(datadir + "CMS_psi2S_cross_section.dat");
  const auto jpsi_CMS_cs = readData<CrossSectionData>(datadir + "CMS_jpsi_cross_section.dat");
  const auto chic2_ATLAS_cs = readData<CrossSectionData>(datadir + "ATLAS_chic2_cross_section.dat");
  const auto chic1_ATLAS_cs = readData<CrossSectionData>(datadir + "ATLAS_chic1_cross_section.dat");
  const auto chic_ratio_CMS_cs = readData<CrossSectionData>(datadir + "CMS_chic_ratio_cross_section.dat");

  // // polarization data
  const auto psi2S_CMS_pol = readData<PolarizationData>(datadir + "CMS_psi2S_polarization.dat");
  const auto jpsi_CMS_pol = readData<PolarizationData>(datadir + "CMS_jpsi_polarization.dat");

  if (!useCosthRatios) {
    return GlobalLikelihood(psi2S_ATLAS_cs, psi2S_CMS_cs, chic2_ATLAS_cs, chic1_ATLAS_cs,
                            jpsi_CMS_cs, chic_ratio_CMS_cs, psi2S_CMS_pol, jpsi_CMS_pol);
  }

  const auto chic_ratios_CMS_pol = read_costh_ratios("chic_ratios.conf", datadir);

  return GlobalLikelihood(psi2S_ATLAS_cs, psi2S_CMS_cs, chic2_ATLAS_cs, chic1_ATLAS_cs,
                          jpsi_CMS_cs, chic_ratio_CMS_cs, psi2S_CMS_pol, jpsi_CMS_pol,
                          chic_ratios_CMS_pol);
}


template<typename LLH, size_t N>
std::vector<std::pair<int, double>> getCovRedFactors(const LLH& llh, const std::array<const char*, N> params,
                                                     const double redFact=5.0)
{
  std::vector<std::pair<int, double>> redFactors;
  for (const auto* par : params) {
    redFactors.emplace_back(llh.getParIdx(par), 1 / (redFact * redFact));
  }

  return redFactors;
}

constexpr std::array<const char*, 25> lambdaContNuissPars = {
  "sigma_psip", "sigma_chic2", "sigma_chic1", "sigma_jpsi",
  "f_long_psi", "gamma", "beta_long_psi", "beta_trans_psi", "beta_long_c1", "beta_trans_c1",
  "beta_long_c2", "beta_trans_c2",
  "br_psip_dp", "br_psip_mm", "br_psip_c2", "br_psip_c1", "br_psip_jpsi", "br_c2_jpsi",
  "br_c1_jpsi", "br_jpsi_mm", "L_CMS", "L_ATLAS",
  "norm_costh_1", "norm_costh_2", "norm_costh_3"};


void global_fit(const std::string& scanFileName="results/scan_file.root",
                const std::string& graphFileName="results/fit_result_graphs.root",
                const std::string& dataDir="./data/",
                const ScanArguments& scanArgs=ScanArguments{},
                const bool useCosthRatios=true, const bool storeGraphs=true)
{
  GlobalLikelihood likelihood = get_likelihood(useCosthRatios, dataDir);

  LikelihoodFitter fitter;

  fitter.Fit(likelihood);

  TFile* scanFile = new TFile(scanFileName.c_str(), "recreate");

  if (scanArgs.randomScan && scanArgs.nSamplePoints > 0) {
      TTree* scanTree = new TTree("log_like_scan", "log likelihood scan values");

      const auto covRedFactors = getCovRedFactors(likelihood, lambdaContNuissPars);
      // Using maxDeltaLLH = 2u, allows to go to a deltaChi2 value of 25 which
      // should be enough for almost everything
      fitter.RandomScan(likelihood, scanTree, scanArgs.nSamplePoints, covRedFactors, 25);
      scanTree->Write("", TObject::kWriteDelete);
  }

  if (scanArgs.nScan1 > 0 && scanArgs.nScan2 > 0) {
    TTree* scanTree = new TTree("log_like_scan", "log likelihood scan values");

    const auto lambda1 = linspace(scanArgs.flow1, scanArgs.fhigh1, scanArgs.nScan1);
    const auto lambda2 = linspace(scanArgs.flow2, scanArgs.fhigh2, scanArgs.nScan2);

    const ScanSettings scanParameters = {{to_frac(lambda1), "f_long_c1"},
                                         {to_frac(lambda2), "f_long_c2"}};

    fitter.Scan(likelihood, scanParameters, scanTree);
    scanTree->Write("", TObject::kWriteDelete);
  }

  // if (nSamplePoints > 0) {
  // }

  TTree* resultTree = new TTree("fit_result", "fit result information");
  fitter.storeFitResult(resultTree);
  resultTree->Write();

  scanFile->Close();

  if (storeGraphs) {
    TFile* graphFile = new TFile(graphFileName.c_str(), "recreate");
    for (const auto& graph : likelihood.getDataGraphs(fitter.Result())) {
      graph.Write();
    }

    graphFile->Close();
  }

}
