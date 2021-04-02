#include "data_reading.h"

#include "likelihood.h"
#include "fitter.h"
#include "misc_utils.h"

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

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
  "lth_psi", "gamma", "beta_long_psi", "beta_trans_psi", "beta_long_c1", "beta_trans_c1",
  "beta_long_c2", "beta_trans_c2",
  "br_psip_dp", "br_psip_mm", "br_psip_c2", "br_psip_c1", "br_psip_jpsi", "br_c2_jpsi",
  "br_c1_jpsi", "br_jpsi_mm", "L_CMS", "L_ATLAS",
  "norm_costh_1", "norm_costh_2", "norm_costh_3"};


void printResult(const LikelihoodFitter& fitter, const GlobalLikelihood& llh)
{
  const double minChi2 = 2 * fitter.Result().MinFcnValue();
  const int nData = llh.nDataPoints();
  const int nPars = llh.nPars();
  const int nNuiss = llh.nNuissParams();
  // since in the likelihood the nuissance parameters are counted towards the
  // total number of parameters, they have to be added to the data points to
  // arrive at a 0 contribution for each nuissance parameter
  const int ndf = nData + nNuiss - nPars;

  std::ios_base::fmtflags fmtflags(std::cout.flags());

  std::cout << "========================= Fit Results =========================\n"
            << "Number of fitted data points: " << nData << "\n"
            << "Total number of parameters: " << nPars << "\n"
            << "Number of nuissance parameters: " << nNuiss << "\n"
            << "Minimum chi2: " << minChi2 << "\n"
            << " -> chi2 / ndf " << std::setprecision(4) << minChi2 << " / " << ndf
            << " (p = " << TMath::Prob(minChi2, ndf) << ")\n"
            << "===============================================================\n";

  std::cout.flags(fmtflags);
}

GlobalLikelihood get_likelihood(bool useCosthRatios, bool usePsiPol=true, bool clipCorrs=false,
                                const std::string& datadir="./data", const double minPTM=2) {
  const auto data = read_all_data(datadir, minPTM);

  if (!useCosthRatios) {
    return GlobalLikelihood(data.psi2S_ATLAS_cs, data.psi2S_CMS_cs, data.chic2_ATLAS_cs, data.chic1_ATLAS_cs, data.jpsi_CMS_cs,
               data.chic_ratio_CMS_cs, data.psi2S_CMS_pol, data.jpsi_CMS_pol, clipCorrs);
  }

  return GlobalLikelihood(data.psi2S_ATLAS_cs, data.psi2S_CMS_cs, data.chic2_ATLAS_cs, data.chic1_ATLAS_cs, data.jpsi_CMS_cs,
             data.chic_ratio_CMS_cs, data.psi2S_CMS_pol, data.jpsi_CMS_pol, data.chic_costh_ratios_CMS, usePsiPol,
             clipCorrs);
}


void global_fit(const std::string& scanFileName="results/scan_file.root",
                const std::string& graphFileName="results/fit_result_graphs.root",
                const std::string& dataDir="./data/",
                const ScanArguments& scanArgs=ScanArguments{},
                const bool useCosthRatios=true, const bool storeGraphs=true, bool usePsiPol=true,
                const bool clipCorrections=false, const bool fixLambasBestFit=false,
                const double minPtM=2)
{
  auto likelihood = get_likelihood(useCosthRatios, usePsiPol, clipCorrections, dataDir, minPtM);

  if (fixLambasBestFit) {
    LikelihoodFitter fitter;
    std::cout << "Fitting likelihood to obtain best fit values for the longitudinal fractions of chic1 and chic2\n";
    fitter.Fit(likelihood);

    const auto& bfPars = fitter.Result().Parameters();
    likelihood.fixParameter("lth_c2", bfPars[likelihood.getParIdx("lth_c2")]);
    likelihood.fixParameter("lth_c1", bfPars[likelihood.getParIdx("lth_c1")]);
  }

  LikelihoodFitter fitter(true);

  fitter.Fit(likelihood);
  printResult(fitter, likelihood);

  TFile* scanFile = new TFile(scanFileName.c_str(), "recreate");

  if (scanArgs.randomScan && scanArgs.nSamplePoints > 0) {
      TTree* scanTree = new TTree("log_like_scan", "log likelihood scan values");
      fitter.setupScanTree(scanTree, false);

      // do the random sampling scan using different "widths" of the
      // multivariate normal distribution that is used to generate random
      // parameter sampling points. At full width the spread of the parameter
      // values is simply too large to do an efficient scanning. By using
      // narrower multivariate normal distributions (with the same correlations)
      // this efficiency can be increased without generating too many points
      // that are uselessly close to the minimum
      const std::vector<std::pair<double, size_t>> subSampleSettings = {
        {4, scanArgs.nSamplePoints / 8},
        {6.25, scanArgs.nSamplePoints / 8},
        {7.5625, scanArgs.nSamplePoints / 8},
        {9., scanArgs.nSamplePoints / 8},
        {10.5625, scanArgs.nSamplePoints / 8},
        {12.25, scanArgs.nSamplePoints / 8},
        {14.0625, scanArgs.nSamplePoints / 8},
        {16., scanArgs.nSamplePoints / 8}
      };
      for (const auto& subSamp : subSampleSettings) {
        std::cout << "Scanning using multivariate normal distribution with variances reduced by factor " << subSamp.first << "\n";
        fitter.RandomScan(likelihood, scanTree, subSamp.second, 1 / subSamp.first, DeltaChi2(2.2977));
      }
      scanTree->Write("", TObject::kWriteDelete);
  }

  if (scanArgs.nScan1 > 0 && scanArgs.nScan2 > 0) {
    TTree* scanTree = new TTree("log_like_scan", "log likelihood scan values");
    fitter.setupScanTree(scanTree, true);

    const auto lambda1 = linspace(scanArgs.flow1, scanArgs.fhigh1, scanArgs.nScan1);
    const auto lambda2 = linspace(scanArgs.flow2, scanArgs.fhigh2, scanArgs.nScan2);

    const ScanSettings scanParameters = {{lambda1, "lth_c1"},
                                         {lambda2, "lth_c2"}};

    fitter.Scan(likelihood, scanParameters, scanTree);
    scanTree->Write("", TObject::kWriteDelete);
  }

  TTree* resultTree = new TTree("fit_result", "fit result information");
  fitter.storeFitResult(resultTree);
  resultTree->Write();

  auto* parIdxTree = likelihood.storeParameterIndices();
  parIdxTree->Write();

  scanFile->Close();

  if (storeGraphs) {
    TFile* graphFile = new TFile(graphFileName.c_str(), "recreate");
    for (const auto& graph : likelihood.getDataGraphs(fitter.Result())) {
      graph.Write();
    }

    for (const auto& model : likelihood.getBestFitModels(fitter.Result())) {
      model.Write();
    }

    graphFile->Close();
  }




}
