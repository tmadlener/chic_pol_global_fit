#ifndef LIKELIHOOD_NRQCD_H
#define LIKELIHOOD_NRQCD_H

#include "constants.h"
#include "data_structures.h"
#include "likelihood_helpers.h"

#include "Fit/FitResult.h"
#include "Fit/ParameterSettings.h"
#include "TTree.h"

#include <stdexcept>
#include <vector>

// TODO: proper setup of fit parameters and compile time map
#include "simple_compile_time_map.h"
constexpr std::array<ParameterIndex, 1> PARAMETERS = {{{"foopar", 0}}};
constexpr static int IPAR(const char* name) { return getParIdx(PARAMETERS, name); }
constexpr static int NPARS(int max = 0, int index = 0) {
  return (index == PARAMETERS.size()) ? max + 1 : // +1 due to 0-indexing
      PARAMETERS[index].second > max ? NPARS(PARAMETERS[index].second, index + 1) : NPARS(max, index + 1);
}

template <size_t N = NPARS()>
void parametersIndicesAsTTree(TTree* tree, const std::array<ParameterIndex, N>& parameters = PARAMETERS) {
  std::array<int, N> indices;
  int i = 0;
  for (const auto& par : parameters) {
    indices[i] = par.second;
    tree->Branch(par.first, &indices[i]);
    i++;
  }
  tree->Fill();
}

// ----------------------------------------- end params stuff

class GlobalLikelihoodNRQCD {
public:
  /**
   * Constructor with signature that is compatible with the GlobalLikelihood
   * constructor signature
   * */
  GlobalLikelihoodNRQCD(const CrossSectionMeasurement& psi2S_ATLAS, const CrossSectionMeasurement& psi2S_CMS,
                        const CrossSectionMeasurement& chic2_ATLAS, const CrossSectionMeasurement& chic1_ATLAS,
                        const CrossSectionMeasurement& jpsi_CMS, const CrossSectionMeasurement& chic_ratio_CMS,
                        const PolarizationMeasurement& psi2S_CMS_p, const PolarizationMeasurement& jpsi_CMS_p,
                        const std::vector<PtCosthRatioMeasurement> chic_costh_ratios_CMS, bool, bool) :
      m_psi2S_ATLAS_cs(psi2S_ATLAS),
      m_psi2S_CMS_cs(psi2S_CMS),
      m_chic2_ATLAS_cs(chic2_ATLAS),
      m_chic1_ATLAS_cs(chic1_ATLAS),
      m_jpsi_CMS_cs(jpsi_CMS),
      m_chic_ratio_CMS_cs(chic_ratio_CMS),
      m_psi2S_CMS_pol(psi2S_CMS_p),
      m_jpsi_CMS_pol(jpsi_CMS_p),
      m_chic_ratios_CMS_pol(chic_costh_ratios_CMS) {
    setupFit();
  }

  /**
   * Constructor without costh ratios. Only here for compatibility with global likelihood. Entering this will terminate
   * the program
   *
   * TODO: clean up the function that builds the likelihoods to be able to remove this
   */
  GlobalLikelihoodNRQCD(const CrossSectionMeasurement&, const CrossSectionMeasurement&, const CrossSectionMeasurement&,
                        const CrossSectionMeasurement&, const CrossSectionMeasurement&, const CrossSectionMeasurement&,
                        const PolarizationMeasurement&, const PolarizationMeasurement&, bool) {
    throw std::runtime_error("Fitting without costh ratios is not implemented for the NRQCD fit");
  }

  /**
   * Evaluate the likelihood with the given parameter values
   */
  double operator()(const double* p) const;

  /**
   * Get the number of parameters
   */
  size_t nPars() const { return m_startParams.size(); }

  /**
   * Get the number of data points
   */
  size_t nDataPoints() const;

  /**
   * Get the number of nuissance parameters
   */
  size_t nNuissParams() const { return m_nuissParams.size(); }

  /**
   * Get the start parameter settings (necessary to pass them to the fitter)
   */
  ParamsSettings getStartParams() const { return m_startParams; }

  /**
   * Fix a parameter to the passed value
   */
  void fixParameter(const char* name, const double val) {
    auto& par = m_startParams[IPAR(name)];
    par.SetValue(val);
    par.Fix();
  }

  TTree* storeParameterIndices() const {
    TTree* parIdxTree = new TTree("parameter_indices", "tree with the indices of the parameters");
    parametersIndicesAsTTree(parIdxTree);
    return parIdxTree;
  }

  /**
   * Cross section models of the direct cross section as a function of pT/M
   */
  static CSModel getPsi2SXSecModel(const double* p);
  static CSModel getChi1XSecModel(const double* p);
  static CSModel getChi2XSecModel(const double* p);
  static CSModel getJpsiXSecModel(const double* p);

  /**
   * Polarization models of the directly produced states as a function of pT/M
   */
  static PolModel getPsiPolModel(const double* p);
  static PolModel getChi1PolModel(const double* p);
  static PolModel getChi2PolModel(const double* p);

private:
  /**
   * set up everything that is needed so that is needed for the fit
   */
  void setupFit();

  /**
   * define the parameters including their starting values
   */
  void defineStartParams();

  /**
   * Set a parameter in the m_startParams vector.
   * This is just a thin wrapper for slightly less typing and better readability
   */
  template <class... Args>
  void setParam(const char* name, Args&&... args) {
    m_startParams[IPAR(name)] = ROOT::Fit::ParameterSettings(name, args...);
  }

  /**
   * add all the nuissance parameters
   */
  void addNuissances();

  template <typename T>
  void addNuissanceParameter(const std::string& name, const T& par);

  CrossSectionMeasurement m_psi2S_ATLAS_cs;
  CrossSectionMeasurement m_psi2S_CMS_cs;
  CrossSectionMeasurement m_chic2_ATLAS_cs;
  CrossSectionMeasurement m_chic1_ATLAS_cs;
  CrossSectionMeasurement m_jpsi_CMS_cs;
  CrossSectionMeasurement m_chic_ratio_CMS_cs;

  PolarizationMeasurement m_psi2S_CMS_pol;
  PolarizationMeasurement m_jpsi_CMS_pol;

  std::vector<PtCosthRatioMeasurement> m_chic_ratios_CMS_pol;

  ParamsSettings m_startParams{NPARS()}; // default initialize
  std::vector<std::pair<int, NuissanceParameter>> m_nuissParams;
};

double GlobalLikelihoodNRQCD::operator()(const double* p) const {
  double loglike = 0;
  // TODO: implementation of likelihood evaluation. Need fit param names first

  // nuissance parameters
  for (const auto& nuissPar : m_nuissParams) {
    loglike += nuissPar.second(p[nuissPar.first]);
  }

  return -loglike;
}

size_t GlobalLikelihoodNRQCD::nDataPoints() const {
  size_t nData = 0;

  nData += m_psi2S_ATLAS_cs.size();
  nData += m_psi2S_CMS_cs.size();
  nData += m_chic2_ATLAS_cs.size();
  nData += m_chic1_ATLAS_cs.size();
  nData += m_jpsi_CMS_cs.size();
  nData += m_chic_ratio_CMS_cs.size();

  nData += m_psi2S_CMS_pol.size();
  nData += m_jpsi_CMS_pol.size();

  for (const auto& chiRatio : m_chic_ratios_CMS_pol) {
    nData += chiRatio.second.size();
  }

  return nData;
}

void GlobalLikelihoodNRQCD::setupFit() {
  defineStartParams();
  addNuissances();
}

void GlobalLikelihoodNRQCD::defineStartParams() {
  // TODO: Once fit parameters have been defined
}

void GlobalLikelihoodNRQCD::addNuissances() {
  addNuissanceParameter("br_psip_dp", B_PSIP_PIPI);
  addNuissanceParameter("br_psip_mm", B_PSIP_MM);
  addNuissanceParameter("br_psip_c2", B_PSIP_CHIC2);
  addNuissanceParameter("br_psip_c1", B_PSIP_CHIC1);
  addNuissanceParameter("br_psip_jpsi", B_PSIP_JPSI);
  addNuissanceParameter("br_c2_jpsi", B_CHIC2_JPSI);
  addNuissanceParameter("br_c1_jpsi", B_CHIC1_JPSI);
  addNuissanceParameter("br_jpsi_mm", B_JPSI_MM);

  addNuissanceParameter("L_CMS", 2.2e-2);
  addNuissanceParameter("L_ATLAS", 1.8e-2);
}

template <typename T>
void GlobalLikelihoodNRQCD::addNuissanceParameter(const std::string& name, const T& par) {
  const auto index = IPAR(name.c_str());
  m_nuissParams.push_back({index, {par}});
}

#endif
