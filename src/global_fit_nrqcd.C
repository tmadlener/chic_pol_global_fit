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

GlobalLikelihoodNRQCD get_likelihood(std::string dataDir, std::string sdcDir,
                                     sdc::SDCType sdcType = sdc::SDCType::LP_NLO) {
  const auto data = read_all_data(dataDir, 2); // pT/M >= 2
  const auto psi_SDC = readPsiSDCs(sdcDir, sdcType);
  const auto chic1_SDC = readChic1SDCs(sdcDir, sdcType);
  const auto chic2_SDC = readChic2SDCs(sdcDir, sdcType);

  const std::string sdcTypeName = sdc::SDCTypeNames[util::to_index(sdcType)];

  // Store the input SDCs for a bit of checking and potentially for plotting
  // later as well
  TFile* inputSDCFile = new TFile("./results/nrqcd_fit/input_sdcs.root", "recreate");
  storeAllContributions(inputSDCFile, psi_SDC, allPsiSDCs, "psi_" + sdcTypeName, PsiSDCsNames);
  storeAllContributions(inputSDCFile, chic1_SDC, allChic1SDCs, "chic1_" + sdcTypeName, Chic1SDCsNames);
  storeAllContributions(inputSDCFile, chic2_SDC, allChic2SDCs, "chic2_" + sdcTypeName, Chic2SDCsNames);

  inputSDCFile->Write("", TObject::kWriteDelete);
  inputSDCFile->Close();
  delete inputSDCFile;

  return GlobalLikelihoodNRQCD(std::move(data), std::move(psi_SDC), std::move(chic1_SDC), std::move(chic2_SDC));
}

void global_fit_nrqcd() {
  auto likelihood = get_likelihood("./data", "../SDCs_Chung");

  LikelihoodFitter fitter;
  fitter.Fit(likelihood);
}
