#include "data_reading.h"
#include "read_sdcs.h"

#include "fitter.h"
#include "likelihood_nrqcd.h"
#include "misc_util.h"
#include "nrqcd_helpers.h"
#include "sdc.h"

#include "TFile.h"

#include <array>
#include <string>
#include <fstream>

void storeSDCasTGraph(TFile* file, const sdc::SDC& sdc, const std::string& name) {
  auto graph = sdc.asTGraph();
  graph.SetName(name.c_str());
  file->cd();
  graph.Write();
}

template <typename Contribs, size_t N>
void storeAllContributions(TFile* file, const sdc::StateSDCs& sdcs, const std::array<Contribs, N>& contributions,
                           std::string baseName, const std::array<const char*, N>& names) {
  file->cd();
  for (const auto contrib : contributions) {
    const auto index = util::to_index(contrib);
    const std::string nameTot = "sdc_total_" + baseName + "_" + names[index];
    const auto sdcTot = sdcs.tot[index];
    storeSDCasTGraph(file, sdcTot, nameTot);

    const std::string nameLong = "sdc_long_" + baseName + "_" + names[index];
    const auto sdcLong = sdcs.lng[index];
    storeSDCasTGraph(file, sdcLong, nameLong);

    // These are not technically SDCs, but here we have all the means to treat
    // them as such and the resulting graphs can be useful
    const auto fLong = sdcLong / sdcTot;
    storeSDCasTGraph(file, fLong, "f_long_" + baseName + "_" + names[index]);

    const auto lth = (1 - 3 * fLong) / (1 + fLong);
    storeSDCasTGraph(file, lth, "lth_" + baseName + "_" + names[index]);
  }
}

GlobalLikelihoodNRQCD get_likelihood(std::string const& dataDir, std::string const& sdcDir,
                                     std::string const& outdir,
                                     double const lpFraction) {
  const auto data = read_all_data(dataDir, 2); // pT/M >= 2
  const auto psi_SDC = readPsiSDCs(sdcDir, lpFraction);
  const auto chic1_SDC = readChic1SDCs(sdcDir, lpFraction);
  const auto chic2_SDC = readChic2SDCs(sdcDir, lpFraction);

  // Store the input SDCs for a bit of checking and potentially for plotting
  // later as well
  const auto sdcFileName = outdir + "/input_sdcs.root";
  TFile* inputSDCFile = new TFile(sdcFileName.c_str(), "recreate");
  storeAllContributions(inputSDCFile, psi_SDC, allPsiSDCs, "psi", PsiSDCsNames);
  storeAllContributions(inputSDCFile, chic1_SDC, allChic1SDCs, "chic1", Chic1SDCsNames);
  storeAllContributions(inputSDCFile, chic2_SDC, allChic2SDCs, "chic2", Chic2SDCsNames);

  inputSDCFile->Write("", TObject::kWriteDelete);
  inputSDCFile->Close();
  delete inputSDCFile;

  return GlobalLikelihoodNRQCD(std::move(data), std::move(psi_SDC), std::move(chic1_SDC), std::move(chic2_SDC));
}

GlobalLikelihoodNRQCD get_likelihood(std::string dataDir, std::string sdcDir,
                                     const std::string& outdir,
                                     sdc::SDCType sdcType = sdc::SDCType::LP_NLO) {
  const auto data = read_all_data(dataDir, 2); // pT/M >= 2
  const auto psi_SDC = readPsiSDCs(sdcDir, sdcType);
  const auto chic1_SDC = readChic1SDCs(sdcDir, sdcType);
  const auto chic2_SDC = readChic2SDCs(sdcDir, sdcType);

  const std::string sdcTypeName = sdc::SDCTypeNames[util::to_index(sdcType)];

  // Store the input SDCs for a bit of checking and potentially for plotting
  // later as well
  const auto sdcFileName = outdir + "/input_sdcs.root";
  TFile* inputSDCFile = new TFile(sdcFileName.c_str(), "recreate");
  storeAllContributions(inputSDCFile, psi_SDC, allPsiSDCs, "psi_" + sdcTypeName, PsiSDCsNames);
  storeAllContributions(inputSDCFile, chic1_SDC, allChic1SDCs, "chic1_" + sdcTypeName, Chic1SDCsNames);
  storeAllContributions(inputSDCFile, chic2_SDC, allChic2SDCs, "chic2_" + sdcTypeName, Chic2SDCsNames);

  inputSDCFile->Write("", TObject::kWriteDelete);
  inputSDCFile->Close();
  delete inputSDCFile;

  return GlobalLikelihoodNRQCD(std::move(data), std::move(psi_SDC), std::move(chic1_SDC), std::move(chic2_SDC));
}

void printResult(const LikelihoodFitter& fitter, const GlobalLikelihoodNRQCD& llh,
                 const std::string& outdir) {
  const double minChi2 = 2 * fitter.Result().MinFcnValue();
  const int nData = llh.nDataPoints();
  // const int nPars = llh.nPars();
  const int nPars = fitter.Result().NFreeParameters();
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

  std::ofstream resultJson(outdir + "/fit_result_info.json");
  if (resultJson) {
    resultJson << "{\n"
               << "  \"n_data\": " << nData << ",\n"
               << "  \"n_pars\": " << nPars << ",\n"
               << "  \"n_nuis\": " << nNuiss << ",\n"
               << "  \"chi2\": " << minChi2 << ",\n"
               << "  \"ndf\": " << ndf << "\n}";
  }
}

/**
 * Set the start parameters or fix certain parameters to the values stored in
 * the settings file. Allowed syntax:
 *
 * - empty lines and lines starting with '#' are ignored
 * - To fix a parameter use
 *
 *   fix <paramName> <paramValue> [# comments starting with '#' are allowed]
 *
 * - To set a start value for a parameter (and optionally also an initial step
 *   size) use
 *
 *   set <paramName> <paramValue> [<paramStepSize>] [# comments starting with '#' are allowed]
 *
 */
template<typename LLH>
void setParameters(LLH& llh, const std::string& paramsSetFN) {
  std::ifstream setfile(paramsSetFN);
  if (!setfile) {
    std::cerr << "Cannot open parameters settings file: \'" << paramsSetFN <<"\'. Using the default parameter settings" << std::endl;
    return;
  }

  std::string line;
  std::string paramName, paramSet;
  double value, error;
  while(std::getline(setfile, line)) {
    if (line[0] == '#' || line.empty()) continue;

    std::stringstream linestr{line};
    linestr >> paramSet >> paramName;
    if (paramSet == "fix") {
      linestr >> value;
      llh.fixParameter(paramName, value);
    } else if (paramSet == "set") {
      linestr >> value >> error;
      llh.setStartParam(paramName, value, error);
    }
  }
}

#if DO_COMBINATION_FIT
void global_fit_nrqcd(const std::string& outdir,
                      const std::string& dataDir,
                      const std::string& sdcDir,
                      const double lpFraction,
                      const std::string& paramsSettingsFile="") {
  auto likelihood = get_likelihood(dataDir, sdcDir, outdir, lpFraction);
#else
void global_fit_nrqcd(const std::string& outdir,
                      const std::string& dataDir,
                      const std::string& sdcDir,
                      sdc::SDCType sdcType=sdc::SDCType::NLO,
                      const std::string& paramsSettingsFile="") {
  auto likelihood = get_likelihood(dataDir, sdcDir, outdir, sdcType);
#endif

  if (!paramsSettingsFile.empty()) {
    setParameters(likelihood, paramsSettingsFile);
  }

  LikelihoodFitter fitter;
  fitter.Fit(likelihood);
  printResult(fitter, likelihood, outdir);

  const auto resultsFileName = outdir + "/fit_results_nrqcd_global_fit.root";
  const auto graphFileName = outdir + "/fit_graphs_and_models_nrqcd_global_fit.root";

  TFile* outFile = new TFile(resultsFileName.c_str(), "recreate");
  TTree* resultTree = new TTree("fit_result", "fit result information");
  fitter.storeFitResult(resultTree);
  resultTree->Write();

  auto* parIdxTree = likelihood.storeParameterIndices();
  parIdxTree->Write();

  outFile->Close();

  TFile* graphFile = new TFile(graphFileName.c_str(), "recreate");
  for (const auto& graph : likelihood.getDataGraphs(fitter.Result())) {
    graph.Write();
  }

  for (const auto& model : likelihood.getBestFitModels(fitter.Result())) {
    model.Write();
  }

  graphFile->Close();
}
