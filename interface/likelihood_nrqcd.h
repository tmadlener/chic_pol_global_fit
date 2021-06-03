#ifndef LIKELIHOOD_NRQCD_H
#define LIKELIHOOD_NRQCD_H

#include "constants.h"
#include "data_structures.h"
#include "fit_parameters_nrqcd.h"
#include "likelihood_helpers.h"
#include "misc_util.h"
#include "misc_utils.h"
#include "nrqcd_helpers.h"
#include "sdc.h"

#include "Fit/FitResult.h"
#include "Fit/ParameterSettings.h"

#include <stdexcept>
#include <vector>
#include <iostream>

// ----------------------------------------- end params stuff

class GlobalLikelihoodNRQCD {
public:
  GlobalLikelihoodNRQCD(const sdc::StateSDCs&& psi_SDC, const sdc::StateSDCs&& chic1_SDC,
                        const sdc::StateSDCs&& chic2_SDC) :
    m_SDC_psi(std::move(psi_SDC)),
    m_SDC_chic1(std::move(chic1_SDC)),
    m_SDC_chic2(std::move(chic2_SDC)) {}

  /**
   * Constructor with signature that is compatible with the GlobalLikelihood
   * constructor signature
   * */
  GlobalLikelihoodNRQCD(const GlobalFitData&& data, const sdc::StateSDCs&& psi_SDC, const sdc::StateSDCs&& chic1_SDC,
                        const sdc::StateSDCs&& chic2_SDC) :
      m_psi2S_ATLAS_cs(std::move(data.psi2S_ATLAS_cs)),
      m_psi2S_CMS_cs(std::move(data.psi2S_CMS_cs)),
      m_chic2_ATLAS_cs(std::move(data.chic2_ATLAS_cs)),
      m_chic1_ATLAS_cs(std::move(data.chic1_ATLAS_cs)),
      m_jpsi_CMS_cs(std::move(data.jpsi_CMS_cs)),
      m_chic_ratio_CMS_cs(std::move(data.chic_ratio_CMS_cs)),
      m_psi2S_CMS_pol(std::move(data.psi2S_CMS_pol)),
      m_jpsi_CMS_pol(std::move(data.jpsi_CMS_pol)),
      m_chic_ratios_CMS_pol(std::move(data.chic_costh_ratios_CMS)),
      m_SDC_psi(std::move(psi_SDC)),
      m_SDC_chic1(std::move(chic1_SDC)),
      m_SDC_chic2(std::move(chic2_SDC)) {
    setupFit();
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
    std::cout << "Fixing parameter " << name << " to " << val << std::endl;
    auto& par = m_startParams[IPAR(name)];
    par.SetValue(val);
    par.Fix();
  }
  void fixParameter(const std::string& name, const double val) {
    fixParameter(name.c_str(), val);
  }

  /**
   * Set a start parameter value for a parameter. If no error value is passed,
   * determine it from the passed value
   */
  void setStartParam(const char* name, const double val, const double err=0) {
    std::cout << "Setting start value for parameter " << name << " to " << val << " (err = " << err << ")" << std::endl;
    auto& par = m_startParams[IPAR(name)];
    par.SetValue(val);
    par.SetStepSize(err != 0 ? err : val);
  }
  void setStartParam(const std::string& name, const double val, const double err=0) {
    setStartParam(name.c_str(), val, err);
  }

  TTree* storeParameterIndices() const {
    TTree* parIdxTree = new TTree("parameter_indices", "tree with the indices of the parameters");
    parametersIndicesAsTTree(parIdxTree);
    return parIdxTree;
  }

  /**
   * Get the data graphs after applying corrections (for the cross section
   * measurements)
   */
  std::vector<TGraphAsymmErrors> getDataGraphs(const ROOT::Fit::FitResult& fitResult) const;

  /**
   * Get the best fit models as a function of pT/M, respectively as a function
   * of costh for the chic costh ratios
   */
  std::vector<TF1> getBestFitModels(const ROOT::Fit::FitResult& fitResult) const;

  /**
   * Cross section models of the direct cross section as a function of pT/M
   */
  CSModel getPsi2SXSecModel(const double* p) const;
  CSModel getChi1XSecModel(const double* p) const;
  CSModel getChi2XSecModel(const double* p) const;
  CSModel getJpsiXSecModel(const double* p) const;

  /**
   * Polarization models of the directly produced states as a function of pT/M
   */
  PolModel getPsi2SPolModel(const double* p) const;
  PolModel getJpsiPolModel(const double* p) const;
  PolModel getChi1PolModel(const double* p) const;
  PolModel getChi2PolModel(const double* p) const;

private:
  /**
   * set up everything that is needed so that is needed for the fit
   */
  void setupFit();

  /**
   * define the parameters including their default starting values
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

  sdc::StateSDCs m_SDC_psi;
  // Total and longitudinal chic1 SDCs as observed after their decay into a
  // J/psi
  sdc::StateSDCs m_SDC_chic1;
  sdc::StateSDCs m_SDC_chic2;

  ParamsSettings m_startParams{NPARS()}; // default initialize
  std::vector<std::pair<int, NuissanceParameter>> m_nuissParams;
};

double GlobalLikelihoodNRQCD::operator()(const double* p) const {
  double loglike = 0;

  // pulling the different parameters out of the p pointer here to guarantee,
  // that the call to IPAR is evaluated at compile time
  const double L_CMS = p[IPAR("L_CMS")];
  const double L_ATLAS = p[IPAR("L_ATLAS")];
  const double br_psip_mm = p[IPAR("br_psip_mm")];
  const double br_psip_dp = p[IPAR("br_psip_dp")];
  const double br_jpsi_mm = p[IPAR("br_jpsi_mm")];

  const auto psi2SXSecModel = getPsi2SXSecModel(p);
  const auto psi2SPolModel = getPsi2SPolModel(p);

  const Identity<double> id;

  // psi(2S) CMS
  loglike +=
      loglikeCrossSection(m_psi2S_CMS_cs, {psi2SXSecModel}, {psi2SPolModel}, {id}, {1.0}, L_CMS * br_psip_mm, M_PSI2S);

  // psi(2S) ATLAS
  loglike += loglikeCrossSection(m_psi2S_ATLAS_cs, {psi2SXSecModel}, {psi2SPolModel}, {id}, {1.0},
                                 L_ATLAS * br_psip_dp * br_jpsi_mm, M_PSI2S);

  const double br_c2_jpsi = p[IPAR("br_c2_jpsi")];
  const double br_psip_c2 = p[IPAR("br_psip_c2")];

  const auto chic2XSecModel = getChi2XSecModel(p);
  const auto chi2PolModel = getChi2PolModel(p);

  const std::vector<CSModel> chi2CSModels = {chic2XSecModel, psi2SXSecModel};
  const std::vector<PolModel> chi2PolModels = {chi2PolModel, psi2SPolModel};
  const std::vector<PolFeedDownTrafo> chi2FDTrafos = {id, lambdaPsiToChi2};
  const std::vector<double> chi2FDFracs = {1.0, B_PSIP_CHIC2[0] / br_psip_c2};

  // chic2 ATLAS
  loglike += loglikeCrossSection(m_chic2_ATLAS_cs, chi2CSModels, chi2PolModels, chi2FDTrafos, chi2FDFracs,
                                 L_ATLAS * br_c2_jpsi * br_jpsi_mm, M_CHIC2);

  const double br_c1_jpsi = p[IPAR("br_c1_jpsi")];
  const double br_psip_c1 = p[IPAR("br_psip_c1")];
  const auto chic1XSecModel = getChi1XSecModel(p);
  const auto chi1PolModel = getChi1PolModel(p);

  const std::vector<CSModel> chi1CSModels = {chic1XSecModel, psi2SXSecModel};
  const std::vector<PolModel> chi1PolModels = {chi1PolModel, psi2SPolModel};
  const std::vector<PolFeedDownTrafo> chi1FDTrafos = {id, lambdaPsiToChi1};
  const std::vector<double> chi1FDFracs = {1.0, B_PSIP_CHIC1[0] / br_psip_c1};

  // chic1 ATLAS
  loglike += loglikeCrossSection(m_chic1_ATLAS_cs, chi1CSModels, chi1PolModels, chi1FDTrafos, chi1FDFracs,
                                 L_ATLAS * br_c1_jpsi * br_jpsi_mm, M_CHIC1);

  const double br_psip_jpsi = p[IPAR("br_psip_jpsi")];

  const auto jpsiXSecModel = getJpsiXSecModel(p);
  const auto jpsiPolModel = getJpsiPolModel(p);

  // TODO: make it possible to compose one function that reflects the chi
  // polarization dependening on pt/M. Currently, we it is necessary to be very
  // general here and list all the contributions to the j/psi individually,
  // considering double feed-down accordingly. NOTE: The lambdas are only
  // affected in the feed-down decay from psi(2S) -> chi

  const std::vector<CSModel> jpsiCSModels = {jpsiXSecModel,  chic1XSecModel, chic2XSecModel,
                                             psi2SXSecModel, psi2SXSecModel, psi2SXSecModel};
  const std::vector<PolModel> jpsiPolModels = {jpsiPolModel,  chi1PolModel,  chi2PolModel,
                                               psi2SPolModel, psi2SPolModel, psi2SPolModel};
  const std::vector<PolFeedDownTrafo> jpsiFDTrafos = {id, id, id, id, lambdaPsiToChi1, lambdaPsiToChi2};

  const std::vector<double> jpsiFDFracs = {1.0,
                                           B_CHIC1_JPSI[0] / br_c1_jpsi,
                                           B_CHIC2_JPSI[0] / br_c2_jpsi,
                                           B_PSIP_JPSI[0] / br_psip_jpsi,
                                           B_PSIP_CHIC1[0] * B_CHIC1_JPSI[0] / br_c1_jpsi / br_psip_c1,
                                           B_PSIP_CHIC2[0] * B_CHIC2_JPSI[0] / br_c2_jpsi / br_psip_c2};

  // jpsi CMS
  loglike += loglikeCrossSection(m_jpsi_CMS_cs, jpsiCSModels, jpsiPolModels, jpsiFDTrafos, jpsiFDFracs,
                                 L_CMS * br_jpsi_mm, M_JPSI);

  // chic ratio CMS
  loglike +=
      loglikeCrossSectionRatio(m_chic_ratio_CMS_cs, chi2CSModels, chi1CSModels, chi2PolModels, chi1PolModels,
                               chi2FDTrafos, chi1FDTrafos, chi2FDFracs, chi1FDFracs, br_psip_c2 / br_psip_c1, M_JPSI);

  // ======================= POLARIZATONS ================================

  // psi(2S) polarization CMS
  loglike += loglikePolarization(m_psi2S_CMS_pol, {psi2SXSecModel}, {psi2SPolModel}, {id}, {1.0}, M_PSI2S);

  // jpsi polarization CMS
  loglike += loglikePolarization(m_jpsi_CMS_pol, jpsiCSModels, jpsiPolModels, jpsiFDTrafos, jpsiFDFracs, M_JPSI);

  // chic2 / chic1 costh ratios
  const std::array<double, 3> normalizations = {p[IPAR("norm_costh_1")], p[IPAR("norm_costh_2")],
                                                p[IPAR("norm_costh_3")]};

  for (size_t iPtBin = 0; iPtBin < normalizations.size(); ++iPtBin) {
    const auto& ratioData = m_chic_ratios_CMS_pol[iPtBin];
    const auto [xs1, lambda1] =
        crossSecAndLambda(ratioData.first, 0.5 / M_JPSI, chi1CSModels, chi1PolModels, chi1FDTrafos, chi1FDFracs);

    const auto [xs2, lambda2] =
        crossSecAndLambda(ratioData.first, 0.5 / M_JPSI, chi2CSModels, chi2PolModels, chi2FDTrafos, chi2FDFracs);

    CosthRatioModel costhRatioModel = [lambda2, lambda1, normalizations, iPtBin](double costh) {
      return costhRatio(costh, lambda2, lambda1, normalizations[iPtBin]);
    };

    loglike += loglikeCosthRatio(ratioData.second, costhRatioModel);
  }

  // nuissance parameters
  for (const auto& nuissPar : m_nuissParams) {
    loglike += nuissPar.second(p[nuissPar.first]);
  }

  return -loglike;
}

CSModel GlobalLikelihoodNRQCD::getPsi2SXSecModel(const double* p) const {
  const auto l_3S1_1_psip = p[IPAR("l_3S1_1_psip")];
  const auto l_1S0_8_psip = p[IPAR("l_1S0_8_psip")];
  const auto l_3PJ_8_psip = p[IPAR("l_rr_3PJ_8_1S0_8_psip_jpsi")] * p[IPAR("l_r_3PJ_8_1S0_8_jpsi")] * l_1S0_8_psip;
  const auto l_3S1_8_psip = p[IPAR("l_rr_3S1_8_1S0_8_psip_jpsi")] * p[IPAR("l_r_3S1_8_1S0_8_jpsi")] * l_1S0_8_psip;

  return [=](double ptm) {
    std::vector<double> ldmes(4, 0);
    ldmes[util::to_index(PsiSDCs::s3S1_1)] = l_3S1_1_psip;
    ldmes[util::to_index(PsiSDCs::s3S1_8)] = l_3S1_8_psip;
    ldmes[util::to_index(PsiSDCs::s3PJ_8)] = l_3PJ_8_psip;
    ldmes[util::to_index(PsiSDCs::s1S0_8)] = l_1S0_8_psip;

    return m_SDC_psi.tot(ptm, ldmes);
  };
}

CSModel GlobalLikelihoodNRQCD::getJpsiXSecModel(const double* p) const {
  const auto l_3S1_1_jpsi = p[IPAR("l_3S1_1_jpsi")];
  const auto l_1S0_8_jpsi = p[IPAR("l_1S0_8_jpsi")];
  const auto l_3PJ_8_jpsi = p[IPAR("l_r_3PJ_8_1S0_8_jpsi")] * l_1S0_8_jpsi;
  const auto l_3S1_8_jpsi = p[IPAR("l_r_3S1_8_1S0_8_jpsi")] * l_1S0_8_jpsi;

  return [=](double ptm) {
    std::vector<double> ldmes(4, 0);
    ldmes[util::to_index(PsiSDCs::s3S1_1)] = l_3S1_1_jpsi;
    ldmes[util::to_index(PsiSDCs::s3S1_8)] = l_3S1_8_jpsi;
    ldmes[util::to_index(PsiSDCs::s3PJ_8)] = l_3PJ_8_jpsi;
    ldmes[util::to_index(PsiSDCs::s1S0_8)] = l_1S0_8_jpsi;

    return m_SDC_psi.tot(ptm, ldmes);
  };
}

CSModel GlobalLikelihoodNRQCD::getChi1XSecModel(const double* p) const {
  const auto l_3S1_8_c1 = p[IPAR("l_3S1_8_c0")] * 3;
  const auto l_3P1_1_c1 = p[IPAR("l_3P0_1_c0")] * 3;

  return [=](double ptm) {
    std::vector<double> ldmes(2, 0);
    ldmes[util::to_index(Chic1SDCs::s3P1_1)] = l_3P1_1_c1;
    ldmes[util::to_index(Chic1SDCs::s3S1_8)] = l_3S1_8_c1;

    return m_SDC_chic1.tot(ptm, ldmes);
  };
}

CSModel GlobalLikelihoodNRQCD::getChi2XSecModel(const double* p) const {
  const auto l_3S1_8_c2 = p[IPAR("l_3S1_8_c0")] * 5;
  const auto l_3P2_1_c2 = p[IPAR("l_3P0_1_c0")] * 5;

  return [=](double ptm) {
    std::vector<double> ldmes(2, 0);
    ldmes[util::to_index(Chic2SDCs::s3P2_1)] = l_3P2_1_c2;
    ldmes[util::to_index(Chic2SDCs::s3S1_8)] = l_3S1_8_c2;

    return m_SDC_chic2.tot(ptm, ldmes);
  };
}

PolModel GlobalLikelihoodNRQCD::getChi1PolModel(const double* p) const {
  const auto l_3S1_8_c1 = p[IPAR("l_3S1_8_c0")] * 3;
  const auto l_3P1_1_c1 = p[IPAR("l_3P0_1_c0")] * 3;

  return [=](double ptm) {
    std::vector<double> ldmes(2, 0);
    ldmes[util::to_index(Chic1SDCs::s3P1_1)] = l_3P1_1_c1;
    ldmes[util::to_index(Chic1SDCs::s3S1_8)] = l_3S1_8_c1;

    const auto total = m_SDC_chic1.tot(ptm, ldmes);
    const auto longitudinal = m_SDC_chic1.lng(ptm, ldmes);
    return lambdath(longitudinal / total);
  };
}

PolModel GlobalLikelihoodNRQCD::getChi2PolModel(const double* p) const {
  const auto l_3S1_8_c2 = p[IPAR("l_3S1_8_c0")] * 5;
  const auto l_3P2_1_c2 = p[IPAR("l_3P0_1_c0")] * 5;

  return [=](double ptm) {
    std::vector<double> ldmes(2, 0);
    ldmes[util::to_index(Chic2SDCs::s3P2_1)] = l_3P2_1_c2;
    ldmes[util::to_index(Chic2SDCs::s3S1_8)] = l_3S1_8_c2;

    const auto total = m_SDC_chic2.tot(ptm, ldmes);
    const auto longitudinal = m_SDC_chic2.lng(ptm, ldmes);
    return lambdath(longitudinal / total);
  };
}

PolModel GlobalLikelihoodNRQCD::getPsi2SPolModel(const double* p) const {
  const auto l_3S1_1_psip = p[IPAR("l_3S1_1_psip")];
  const auto l_1S0_8_psip = p[IPAR("l_1S0_8_psip")];
  const auto l_3PJ_8_psip = p[IPAR("l_rr_3PJ_8_1S0_8_psip_jpsi")] * p[IPAR("l_r_3PJ_8_1S0_8_jpsi")] * l_1S0_8_psip;
  const auto l_3S1_8_psip = p[IPAR("l_rr_3S1_8_1S0_8_psip_jpsi")] * p[IPAR("l_r_3S1_8_1S0_8_jpsi")] * l_1S0_8_psip;

  return [=](double ptm) {
    std::vector<double> ldmes(4, 0);
    ldmes[util::to_index(PsiSDCs::s3S1_1)] = l_3S1_1_psip;
    ldmes[util::to_index(PsiSDCs::s3S1_8)] = l_3S1_8_psip;
    ldmes[util::to_index(PsiSDCs::s3PJ_8)] = l_3PJ_8_psip;
    ldmes[util::to_index(PsiSDCs::s1S0_8)] = l_1S0_8_psip;

    const auto total = m_SDC_psi.tot(ptm, ldmes);
    const auto longitudinal = m_SDC_psi.lng(ptm, ldmes);
    return lambdath(longitudinal / total);
  };
}

PolModel GlobalLikelihoodNRQCD::getJpsiPolModel(const double* p) const {
  const auto l_3S1_1_jpsi = p[IPAR("l_3S1_1_jpsi")];
  const auto l_1S0_8_jpsi = p[IPAR("l_1S0_8_jpsi")];
  const auto l_3PJ_8_jpsi = p[IPAR("l_r_3PJ_8_1S0_8_jpsi")] * l_1S0_8_jpsi;
  const auto l_3S1_8_jpsi = p[IPAR("l_r_3S1_8_1S0_8_jpsi")] * l_1S0_8_jpsi;

  return [=](double ptm) {
    std::vector<double> ldmes(4, 0);
    ldmes[util::to_index(PsiSDCs::s3S1_1)] = l_3S1_1_jpsi;
    ldmes[util::to_index(PsiSDCs::s3S1_8)] = l_3S1_8_jpsi;
    ldmes[util::to_index(PsiSDCs::s3PJ_8)] = l_3PJ_8_jpsi;
    ldmes[util::to_index(PsiSDCs::s1S0_8)] = l_1S0_8_jpsi;

    const auto total = m_SDC_psi.tot(ptm, ldmes);
    const auto longitudinal = m_SDC_psi.lng(ptm, ldmes);
    return lambdath(longitudinal / total);
  };
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
  // TODO: Get reasonable start values
  setParam("l_3S1_8_c0", 0.002, 1);
  setParam("l_3P0_1_c0", 0.007, 1);

  setParam("l_3S1_1_jpsi", 10, 1);   // TODO: fix these as nuisances to theory?
  setParam("l_3S1_1_psip", 0.05, 1); // TODO: fix these as nuisances to theory?

  setParam("l_1S0_8_jpsi", 1e-6, 0.1);
  setParam("l_1S0_8_psip", 1e-6, 0.1);

  setParam("l_r_3PJ_8_1S0_8_jpsi", 30, 5);
  setParam("l_r_3S1_8_1S0_8_jpsi", 30, 5);
  setParam("l_rr_3PJ_8_1S0_8_psip_jpsi", 1, 0.1);
  setParam("l_rr_3S1_8_1S0_8_psip_jpsi", 1, 0.1);

  setParam("br_psip_dp", 1, 0.01);
  setParam("br_psip_mm", 1, 0.01);
  setParam("br_psip_c2", 1, 0.01);
  setParam("br_psip_c1", 1, 0.01);
  setParam("br_psip_jpsi", 1, 0.01);
  setParam("br_c2_jpsi", 1, 0.01);
  setParam("br_c1_jpsi", 1, 0.01);
  setParam("br_jpsi_mm", 1, 0.01);

  setParam("L_CMS", 1, 0.01);
  setParam("L_ATLAS", 1, 0.01);

  setParam("norm_costh_1", 0.45, 0.1);
  setParam("norm_costh_2", 0.45, 0.1);
  setParam("norm_costh_3", 0.45, 0.1);
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

std::vector<TGraphAsymmErrors> GlobalLikelihoodNRQCD::getDataGraphs(const ROOT::Fit::FitResult& fitResult) const {
  const auto parValues = fitResult.Parameters();
  std::vector<TGraphAsymmErrors> graphs;

  const double L_CMS = parValues[IPAR("L_CMS")];
  const double br_psip_mm = parValues[IPAR("br_psip_mm")];
  const auto psi2SXSecModel = getPsi2SXSecModel(parValues.data());
  const auto psi2SPolModel = getPsi2SPolModel(parValues.data());

  Identity<double> id;

  graphs.push_back(correctedCSGraph(m_psi2S_CMS_cs, {psi2SXSecModel}, {psi2SPolModel}, {id}, {1.0}, L_CMS * br_psip_mm,
                                    M_PSI2S, "psi2S_CMS_cs"));

  const double L_ATLAS = parValues[IPAR("L_ATLAS")];
  const double br_psip_dp = parValues[IPAR("br_psip_dp")];
  const double br_jpsi_mm = parValues[IPAR("br_jpsi_mm")];

  graphs.push_back(correctedCSGraph(m_psi2S_ATLAS_cs, {psi2SXSecModel}, {psi2SPolModel}, {id}, {1.0},
                                    L_ATLAS * br_psip_dp * br_jpsi_mm, M_PSI2S, "psi2S_ATLAS_cs"));

  const double br_c1_jpsi = parValues[IPAR("br_c1_jpsi")];
  const double br_psip_c1 = parValues[IPAR("br_psip_c1")];

  const auto chi1XSecModel = getChi1XSecModel(parValues.data());
  const auto chi1PolModel = getChi1PolModel(parValues.data());

  const std::vector<CSModel> chi1CSModels = {chi1XSecModel, psi2SXSecModel};
  const std::vector<PolModel> chi1PolModels = {chi1PolModel, psi2SPolModel};
  const std::vector<PolFeedDownTrafo> chi1FDTrafos = {id, lambdaPsiToChi1};
  const std::vector<double> chi1FDFracs = {1.0, B_PSIP_CHIC1[0] / br_psip_c1};

  graphs.push_back(correctedCSGraph(m_chic1_ATLAS_cs, chi1CSModels, chi1PolModels, chi1FDTrafos, chi1FDFracs,
                                    L_ATLAS * br_c1_jpsi * br_jpsi_mm, M_CHIC1, "chic1_ATLAS_cs"));

  const double br_c2_jpsi = parValues[IPAR("br_c2_jpsi")];
  const double br_psip_c2 = parValues[IPAR("br_psip_c2")];

  const auto chi2XSecModel = getChi2XSecModel(parValues.data());
  const auto chi2PolModel = getChi2PolModel(parValues.data());

  const std::vector<CSModel> chi2CSModels = {chi2XSecModel, psi2SXSecModel};
  const std::vector<PolModel> chi2PolModels = {chi2PolModel, psi2SPolModel};
  const std::vector<PolFeedDownTrafo> chi2FDTrafos = {id, lambdaPsiToChi2};
  const std::vector<double> chi2FDFracs = {1.0, B_PSIP_CHIC2[0] / br_psip_c2};

  graphs.push_back(correctedCSGraph(m_chic2_ATLAS_cs, chi2CSModels, chi2PolModels, chi2FDTrafos, chi2FDFracs,
                                    L_ATLAS * br_c2_jpsi * br_jpsi_mm, M_CHIC2, "chic2_ATLAS_cs"));
  //
  // chic ratios first have to compute the corrections individually

  const auto chi1Corrections = getCorrectionFactors(m_chic_ratio_CMS_cs, chi1CSModels, chi1PolModels, chi1FDTrafos,
                                                    chi1FDFracs, L_ATLAS * br_c1_jpsi * br_jpsi_mm, M_JPSI);

  const auto chi2Corrections = getCorrectionFactors(m_chic_ratio_CMS_cs, chi2CSModels, chi2PolModels, chi2FDTrafos,
                                                    chi2FDFracs, L_ATLAS * br_c2_jpsi * br_jpsi_mm, M_JPSI);

  std::vector<double> ratioCorrections;
  for (size_t i = 0; i < chi1Corrections.size(); ++i) {
    ratioCorrections.push_back(chi2Corrections[i] / chi1Corrections[i]);
  }

  graphs.push_back(correctGraph(m_chic_ratio_CMS_cs, "chic_ratio_CMS_cs", ratioCorrections));

  const double br_psip_jpsi = parValues[IPAR("br_psip_jpsi")];
  const auto jpsiXSecModel = getJpsiXSecModel(parValues.data());
  const auto jpsiPolModel = getJpsiPolModel(parValues.data());

  const std::vector<CSModel> jpsiCSModels = {jpsiXSecModel,  chi1XSecModel,  chi2XSecModel,
                                             psi2SXSecModel, psi2SXSecModel, psi2SXSecModel};
  const std::vector<PolModel> jpsiPolModels = {jpsiPolModel,  chi1PolModel,  chi2PolModel,
                                               psi2SPolModel, psi2SPolModel, psi2SPolModel};
  const std::vector<PolFeedDownTrafo> jpsiFDTrafos = {id, id, id, id, lambdaPsiToChi1, lambdaPsiToChi2};
  const std::vector<double> jpsiFDFracs = {1.0,
                                           B_CHIC1_JPSI[0] / br_c1_jpsi,
                                           B_CHIC2_JPSI[0] / br_c2_jpsi,
                                           B_PSIP_JPSI[0] / br_psip_jpsi,
                                           B_PSIP_CHIC1[0] * B_CHIC1_JPSI[0] / br_c1_jpsi / br_psip_c1,
                                           B_PSIP_CHIC2[0] * B_CHIC2_JPSI[0] / br_c2_jpsi / br_psip_c2};

  graphs.push_back(correctedCSGraph(m_jpsi_CMS_cs, jpsiCSModels, jpsiPolModels, jpsiFDTrafos, jpsiFDFracs,
                                    L_CMS * br_jpsi_mm, M_JPSI, "jpsi_CMS_cs"));

  // Polarizations
  graphs.push_back(asTGraph(m_psi2S_CMS_pol));
  graphs.back().SetName("psi2S_CMS_pol");

  graphs.push_back(asTGraph(m_jpsi_CMS_pol));
  graphs.back().SetName("jpsi_CMS_pol");

  for (const auto& cthRatioMeas : m_chic_ratios_CMS_pol) {
    graphs.push_back(asTGraph(cthRatioMeas.second));
    std::stringstream gn;
    const int after_dec = (cthRatioMeas.first.ptM - int(cthRatioMeas.first.ptM)) * 100;
    gn << "chic_CMS_pol_ptM_" << int(cthRatioMeas.first.ptM) << "p" << after_dec;
    graphs.back().SetName(gn.str().c_str());
  }

  return graphs;
}

std::vector<TF1> GlobalLikelihoodNRQCD::getBestFitModels(const ROOT::Fit::FitResult& fitResult) const {
  const auto parValues = fitResult.Parameters();
  std::vector<TF1> models;

  const auto psi2SXSecModel = getPsi2SXSecModel(parValues.data());
  const auto psi2SPolModel = getPsi2SPolModel(parValues.data());
  models.push_back(modelAsTF1(psi2SXSecModel, "psip_cs_direct", 2, 40));
  models.push_back(modelAsTF1(psi2SPolModel, "psip_pol_direct", 2, 40));

  const auto chi1XSecModel = getChi1XSecModel(parValues.data());
  const auto chi1PolModel = getChi1PolModel(parValues.data());
  models.push_back(modelAsTF1(chi1XSecModel, "chic1_cs_direct", 2, 40));
  models.push_back(modelAsTF1(chi1PolModel, "chic1_pol_direct", 2, 40));

  const auto chi2XSecModel = getChi2XSecModel(parValues.data());
  const auto chi2PolModel = getChi2PolModel(parValues.data());
  models.push_back(modelAsTF1(chi2XSecModel, "chic2_cs_direct", 2, 40));
  models.push_back(modelAsTF1(chi2PolModel, "chic2_pol_direct", 2, 40));

  const auto jpsiXSecModel = getJpsiXSecModel(parValues.data());
  const auto jpsiPolModel = getJpsiPolModel(parValues.data());

  models.push_back(modelAsTF1(jpsiXSecModel, "jpsi_cs_direct", 2, 40));
  models.push_back(modelAsTF1(psi2SPolModel, "jpsi_pol_direct", 2, 40));

  const Identity<double> id;

  const double br_psip_c1 = parValues[IPAR("br_psip_c1")];
  const auto [chic1XSecFull, chic1PolFull] =
      combineModels({chi1XSecModel, psi2SXSecModel}, {chi1PolModel, psi2SPolModel}, {id, lambdaPsiToChi1},
                    {1.0, B_PSIP_CHIC1[0] / br_psip_c1});

  models.push_back(modelAsTF1(chic1XSecFull, "chic1_cs_full", 2, 40));
  models.push_back(modelAsTF1(chic1PolFull, "chic1_pol_full", 2, 40));

  const double br_psip_c2 = parValues[IPAR("br_psip_c2")];
  const auto [chic2XSecFull, chic2PolFull] =
      combineModels({chi2XSecModel, psi2SXSecModel}, {chi2PolModel, psi2SPolModel}, {id, lambdaPsiToChi2},
                    {1.0, B_PSIP_CHIC2[0] / br_psip_c2});
  models.push_back(modelAsTF1(chic2XSecFull, "chic2_cs_full", 2, 40));
  models.push_back(modelAsTF1(chic2PolFull, "chic2_pol_full", 2, 40));

  // chic ratio cross section
  const CSModel chicRatioModel = [chic1XSecFull, chic2XSecFull](double ptm) -> double {
    return chic2XSecFull(ptm) / chic1XSecFull(ptm);
  };
  models.push_back(modelAsTF1(chicRatioModel, "chic_ratio_cs_full", 2, 40));

  const double br_c1_jpsi = parValues[IPAR("br_c1_jpsi")];
  const double br_c2_jpsi = parValues[IPAR("br_c2_jpsi")];
  const double br_psip_jpsi = parValues[IPAR("br_psip_jpsi")];

  const auto [jpsiXSecFull, jpsiPolFull] =
      combineModels({jpsiXSecModel, chic1XSecFull, chic2XSecFull, psi2SXSecModel},
                    {jpsiPolModel, chic1PolFull, chic2PolFull, psi2SPolModel}, {id, id, id, id},
                    {1.0, B_CHIC1_JPSI[0] / br_c1_jpsi, B_CHIC2_JPSI[0] / br_c2_jpsi, B_PSIP_JPSI[0] / br_psip_jpsi});

  models.push_back(modelAsTF1(jpsiXSecFull, "jpsi_cs_full", 2, 40));
  models.push_back(modelAsTF1(jpsiPolFull, "jpsi_pol_full", 2, 40));

  const std::array<double, 3> normalizations = {parValues[IPAR("norm_costh_1")], parValues[IPAR("norm_costh_2")],
                                                parValues[IPAR("norm_costh_3")]};

  for (size_t iPtBin = 0; iPtBin < normalizations.size(); ++iPtBin) {
    const auto& ratioData = m_chic_ratios_CMS_pol[iPtBin];

    const auto [xs1, lambda1] =
        crossSecAndLambda(ratioData.first, 0.5 / M_JPSI, {chi1XSecModel, psi2SXSecModel}, {chi1PolModel, psi2SPolModel},
                          {id, lambdaPsiToChi1}, {1.0, B_PSIP_CHIC1[0] / br_psip_c1});

    const auto [xs2, lambda2] =
        crossSecAndLambda(ratioData.first, 0.5 / M_JPSI, {chi2XSecModel, psi2SXSecModel}, {chi2PolModel, psi2SPolModel},
                          {id, lambdaPsiToChi2}, {1.0, B_PSIP_CHIC2[0] / br_psip_c2});

    const CosthRatioModel costhRatioModel = [lambda2, lambda1, normalizations, iPtBin](double costh) {
      return costhRatio(costh, lambda2, lambda1, normalizations[iPtBin]);
    };

    std::stringstream gn;
    const int after_dec = (ratioData.first.ptM - int(ratioData.first.ptM)) * 100;
    gn << "chic_pol_ptM_" << int(ratioData.first.ptM) << "p" << after_dec;

    models.push_back(modelAsTF1(costhRatioModel, gn.str().c_str(), 0, 1));
  }

  return models;
}

#endif
