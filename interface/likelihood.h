#ifndef H_LIKLIEHOOD__
#define H_LIKLIEHOOD__

#include "data_structures.h"
#include "constants.h"
#include "likelihood_helpers.h"
#include "fit_parameters.h"

#include "Fit/ParameterSettings.h"
#include "Fit/FitResult.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TTree.h"

#include <vector>
#include <array>
#include <utility>
#include <string>
#include <sstream>


class GlobalLikelihood {
public:
  GlobalLikelihood() = default;

  GlobalLikelihood(const CrossSectionMeasurement& psi2S_ATLAS, const CrossSectionMeasurement& psi2S_CMS,
                   const CrossSectionMeasurement& chic2_ATLAS, const CrossSectionMeasurement& chic1_ATLAS,
                   const CrossSectionMeasurement& jpsi_CMS, const CrossSectionMeasurement& chic_ratio_CMS,
                   const PolarizationMeasurement& psi2S_CMS_p, const PolarizationMeasurement& jpsi_CMS_p,
                   const std::vector<PtCosthRatioMeasurement> chic_costh_ratios_CMS,
                   const bool usePsiPol=true, const bool clipCorrs=false) :
    m_psi2S_ATLAS_cs(psi2S_ATLAS), m_psi2S_CMS_cs(psi2S_CMS), m_chic2_ATLAS_cs(chic2_ATLAS),
    m_chic1_ATLAS_cs(chic1_ATLAS), m_jpsi_CMS_cs(jpsi_CMS), m_chic_ratio_CMS_cs(chic_ratio_CMS),
    m_psi2S_CMS_pol(psi2S_CMS_p), m_jpsi_CMS_pol(jpsi_CMS_p), m_chic_ratios_CMS_pol(chic_costh_ratios_CMS),
    m_usePsiPolarizations(usePsiPol), m_clipCorrectionsPhysicalRange(clipCorrs)
  {
    setupFit();
  }

  GlobalLikelihood(const CrossSectionMeasurement& psi2S_ATLAS, const CrossSectionMeasurement& psi2S_CMS,
                   const CrossSectionMeasurement& chic2_ATLAS, const CrossSectionMeasurement& chic1_ATLAS,
                   const CrossSectionMeasurement& jpsi_CMS, const CrossSectionMeasurement& chic_ratio_CMS,
                   const PolarizationMeasurement& psi2S_CMS_p, const PolarizationMeasurement& jpsi_CMS_p,
                   const bool clipCorrs=false) :
    m_psi2S_ATLAS_cs(psi2S_ATLAS), m_psi2S_CMS_cs(psi2S_CMS), m_chic2_ATLAS_cs(chic2_ATLAS),
    m_chic1_ATLAS_cs(chic1_ATLAS), m_jpsi_CMS_cs(jpsi_CMS), m_chic_ratio_CMS_cs(chic_ratio_CMS),
    m_psi2S_CMS_pol(psi2S_CMS_p), m_jpsi_CMS_pol(jpsi_CMS_p), m_clipCorrectionsPhysicalRange(clipCorrs)
  {
    // if no costh ratios are used, some parts of the setups have to be done
    // differently. Most importantly, the start parameters vector has to be
    // correctly resized, and the initialization of the normalizations can be
    // omitted.
    m_useCosthRatios = false;
    m_startParams.resize(NPARS() - 3);
    setupFit();
  }

  /**
   * Evaluate the likelihood with the given parameters
   */
  double operator()(const double* p) const;

  /**
   * Get the number of parameters
   */
  size_t nPars() const { return m_startParams.size(); }

  static int getParIdx(const std::string& name) { return IPAR(name.c_str()); }

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
  void fixParameter(const char* name, const double val)
  {
    auto& par = m_startParams[IPAR(name)];
    par.SetValue(val);
    par.Fix();
  }

  TTree* storeParameterIndices() const
  {
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
  template<class... Args>
  void setParam(const char* name, Args&&... args)
  { m_startParams[IPAR(name)] = ROOT::Fit::ParameterSettings(name, args...); }

  /**
   * add all the nuissance parameters
   */
  void addNuissances();

  template<typename T>
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
  bool m_useCosthRatios{true};
  bool m_usePsiPolarizations{true};
  bool m_clipCorrectionsPhysicalRange{false};

  ParamsSettings m_startParams{NPARS()}; // default initialize
  std::vector<std::pair<int, NuissanceParameter>> m_nuissParams;
};



CSModel GlobalLikelihood::getPsi2SXSecModel(const double* p)
{
  const double sigma_psip = p[IPAR("sigma_psip")];
  const double beta_total_psi = p[IPAR("beta_total_psi")];
  const double gamma = p[IPAR("gamma")];

  return  [sigma_psip, beta_total_psi, gamma] (double ptm) {
    return sig_dir(ptm, sigma_psip, beta_total_psi, gamma);
  };
}

CSModel GlobalLikelihood::getChi1XSecModel(const double* p)
{
  const double sigma_chic1 = p[IPAR("sigma_chic1")];
  const double beta_total_c1 = p[IPAR("beta_total_c1")];
  const double gamma = p[IPAR("gamma")];

  return [sigma_chic1, beta_total_c1, gamma] (double ptm) {
    return sig_dir(ptm, sigma_chic1, beta_total_c1, gamma);
  };
}

CSModel GlobalLikelihood::getChi2XSecModel(const double* p)
{
  const double sigma_chic2 = p[IPAR("sigma_chic2")];
  const double beta_total_c2 = p[IPAR("beta_total_c2")];
  const double gamma = p[IPAR("gamma")];

  return [sigma_chic2, beta_total_c2, gamma] (double ptm) {
    return sig_dir(ptm, sigma_chic2, beta_total_c2, gamma);
  };
}

CSModel GlobalLikelihood::getJpsiXSecModel(const double* p)
{
  const double sigma_jpsi = p[IPAR("sigma_jpsi")];
  const double beta_total_psi = p[IPAR("beta_total_psi")];
  const double gamma = p[IPAR("gamma")];

  return [sigma_jpsi, beta_total_psi, gamma] (double ptm) {
    return sig_dir(ptm, sigma_jpsi, beta_total_psi, gamma);
  };
}

PolModel GlobalLikelihood::getPsiPolModel(const double* p)
{
  const double lth_psi = p[IPAR("lth_psi")];
  const double beta_total_psi = p[IPAR("beta_total_psi")];
  const double beta_long_psi = p[IPAR("beta_long_psi")];
  const double gamma = p[IPAR("gamma")];

  return [lth_psi, beta_total_psi, beta_long_psi, gamma] (double ptm) {
    return lambdaThetaAtPtm(ptm, lth_psi, beta_long_psi, beta_total_psi, gamma);
  };
}

PolModel GlobalLikelihood::getChi1PolModel(const double* p)
{
  const double beta_long_c1 = p[IPAR("beta_long_c1")];
  const double beta_total_c1 = p[IPAR("beta_total_c1")];
  const double lth_c1 = p[IPAR("lth_c1")];
  const double gamma = p[IPAR("gamma")];

  return [lth_c1, beta_total_c1, beta_long_c1, gamma] (double ptm) {
    return lambdaThetaAtPtm(ptm, lth_c1, beta_long_c1, beta_total_c1, gamma);
  };
}

PolModel GlobalLikelihood::getChi2PolModel(const double* p)
{
  const double beta_long_c2 = p[IPAR("beta_long_c2")];
  const double beta_total_c2 = p[IPAR("beta_total_c2")];
  const double lth_c2 = p[IPAR("lth_c2")];
  const double gamma = p[IPAR("gamma")];

  return [lth_c2, beta_total_c2, beta_long_c2, gamma] (double ptm) {
    return lambdaThetaAtPtm(ptm, lth_c2, beta_long_c2, beta_total_c2, gamma);
  };
}

size_t GlobalLikelihood::nDataPoints() const
{
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


double GlobalLikelihood::operator()(const double* p) const
{
  double loglike = 0;

  // pulling the different parameters out of the p pointer here to guarantee, that the call
  // to IPAR is evaluated at compile time
  const double L_CMS = p[IPAR("L_CMS")];
  const double L_ATLAS = p[IPAR("L_ATLAS")];
  const double br_psip_mm = p[IPAR("br_psip_mm")];
  const double br_psip_dp = p[IPAR("br_psip_dp")];
  const double br_jpsi_mm = p[IPAR("br_jpsi_mm")];

  const auto psi2SXSecModel = getPsi2SXSecModel(p);
  const auto psiPolModel = getPsiPolModel(p);

  const Identity<double> id;

  // psi(2S) CMS
  loglike += loglikeCrossSection(m_psi2S_CMS_cs,
                                 {psi2SXSecModel}, {psiPolModel}, {id}, {1.0},
                                 L_CMS * br_psip_mm, M_PSI2S,
                                 m_clipCorrectionsPhysicalRange, {-1.0, 1.0});

  // psi(2S) ATLAS
  loglike += loglikeCrossSection(m_psi2S_ATLAS_cs,
                                 {psi2SXSecModel}, {psiPolModel}, {id}, {1.0},
                                 L_ATLAS * br_psip_dp * br_jpsi_mm, M_PSI2S,
                                 m_clipCorrectionsPhysicalRange, {-1.0, 1.0});

  const double br_c2_jpsi = p[IPAR("br_c2_jpsi")];
  const double br_psip_c2 = p[IPAR("br_psip_c2")];

  const auto chic2XSecModel = getChi2XSecModel(p);
  const auto chi2PolModel = getChi2PolModel(p);

  const std::vector<CSModel> chi2CSModels = {chic2XSecModel, psi2SXSecModel};
  const std::vector<PolModel> chi2PolModels = {chi2PolModel, psiPolModel};
  const std::vector<PolFeedDownTrafo> chi2FDTrafos = {id, lambdaPsiToChi2};
  const std::vector<double> chi2FDFracs = {1.0, B_PSIP_CHIC2[0] / br_psip_c2};

  loglike += loglikeCrossSection(m_chic2_ATLAS_cs, chi2CSModels, chi2PolModels, chi2FDTrafos, chi2FDFracs,
                                 L_ATLAS * br_c2_jpsi * br_jpsi_mm, M_CHIC2,
                                 m_clipCorrectionsPhysicalRange, {-0.6, 1.0});


  const double br_c1_jpsi = p[IPAR("br_c1_jpsi")];
  const double br_psip_c1 = p[IPAR("br_psip_c1")];
  const auto chic1XSecModel = getChi1XSecModel(p);
  const auto chi1PolModel = getChi1PolModel(p);

  const std::vector<CSModel> chi1CSModels = {chic1XSecModel, psi2SXSecModel};
  const std::vector<PolModel> chi1PolModels = {chi1PolModel, psiPolModel};
  const std::vector<PolFeedDownTrafo> chi1FDTrafos = {id, lambdaPsiToChi1};
  const std::vector<double> chi1FDFracs = {1.0, B_PSIP_CHIC1[0] / br_psip_c1};

  loglike += loglikeCrossSection(m_chic1_ATLAS_cs, chi1CSModels, chi1PolModels, chi1FDTrafos, chi1FDFracs,
                                 L_ATLAS * br_c1_jpsi * br_jpsi_mm, M_CHIC1,
                                 m_clipCorrectionsPhysicalRange, {-1./3., 1.0});


  const double br_psip_jpsi = p[IPAR("br_psip_jpsi")];

  const auto jpsiXSecModel = getJpsiXSecModel(p);

  // jpsi direct polarization is the same as psi(2S) direct polarization
  // TODO: make it possible to compose one function that reflects the chi
  // polarization dependening on pt/M. Currently, we it is necessary to be very
  // general here and list all the contributions to the j/psi individually,
  // considering double feed-down accordingly. NOTE: The lambdas are only
  // affected in the feed-down decay from psi(2S) -> chi

  const std::vector<CSModel> jpsiCSModels = {
    jpsiXSecModel, chic1XSecModel, chic2XSecModel, psi2SXSecModel,
    psi2SXSecModel, psi2SXSecModel
  };
  const std::vector<PolModel> jpsiPolModels = {
    psiPolModel, chi1PolModel, chi2PolModel, psiPolModel,
    psiPolModel, psiPolModel
  };
  const std::vector<PolFeedDownTrafo> jpsiFDTrafos = {
    id, id, id, id, lambdaPsiToChi1, lambdaPsiToChi2
  };
  const std::vector<double> jpsiFDFracs = {
    1.0, B_CHIC1_JPSI[0] / br_c1_jpsi, B_CHIC2_JPSI[0] / br_c2_jpsi, B_PSIP_JPSI[0] / br_psip_jpsi,
    B_PSIP_CHIC1[0] * B_CHIC1_JPSI[0] / br_c1_jpsi / br_psip_c1,
    B_PSIP_CHIC2[0] * B_CHIC2_JPSI[0] / br_c2_jpsi / br_psip_c2
  };

  loglike += loglikeCrossSection(m_jpsi_CMS_cs, jpsiCSModels, jpsiPolModels, jpsiFDTrafos, jpsiFDFracs,
                                 L_CMS * br_jpsi_mm, M_JPSI,
                                 m_clipCorrectionsPhysicalRange, {-1.0, 1.0});


  loglike += loglikeCrossSectionRatio(m_chic_ratio_CMS_cs,
                                      chi2CSModels, chi1CSModels, chi2PolModels, chi1PolModels,
                                      chi2FDTrafos, chi1FDTrafos, chi2FDFracs, chi1FDFracs,
                                      br_psip_c2 / br_psip_c1, M_JPSI,
                                      m_clipCorrectionsPhysicalRange, {-0.6, 1.0}, {-1./3., 1.0});


  if (m_usePsiPolarizations) {
    loglike += loglikePolarization(m_psi2S_CMS_pol, {psi2SXSecModel}, {psiPolModel}, {id}, {1.0}, M_PSI2S);

    // same as above for the cross-sections. All contributions have to be considered individually
    loglike += loglikePolarization(m_jpsi_CMS_pol, jpsiCSModels, jpsiPolModels, jpsiFDTrafos, jpsiFDFracs,
                                   M_JPSI);
  }


  if (m_useCosthRatios) {

    const std::array<double, 3> normalizations = {
      p[IPAR("norm_costh_1")],
      p[IPAR("norm_costh_2")],
      p[IPAR("norm_costh_3")]
    };

    // costh ratios CMS
    for (size_t iPtBin = 0; iPtBin < normalizations.size(); ++iPtBin) {
      double lambda1, lambda2;

      const auto& ratioData = m_chic_ratios_CMS_pol[iPtBin];

      std::tie(std::ignore, lambda1) = crossSecAndLambda(ratioData.first, 0.5 / M_JPSI,
                                                         chi1CSModels, chi1PolModels,
                                                         chi1FDTrafos, chi1FDFracs);

      std::tie(std::ignore, lambda2) = crossSecAndLambda(ratioData.first, 0.5 / M_JPSI,
                                                         chi2CSModels, chi2PolModels,
                                                         chi2FDTrafos, chi2FDFracs);

      CosthRatioModel costhRatioModel = [lambda2, lambda1, normalizations, iPtBin] (double costh) {
        return costhRatio(costh, lambda2, lambda1, normalizations[iPtBin]);
      };

      loglike += loglikeCosthRatio(ratioData.second, costhRatioModel);
    }
  }



  // nuissance parameters
  for (const auto& nuissPar : m_nuissParams) {
    loglike += nuissPar.second(p[nuissPar.first]);
  }

  return -loglike;
}

void GlobalLikelihood::setupFit()
{
  defineStartParams();
  addNuissances();

}

void GlobalLikelihood::defineStartParams()
{
  setParam("sigma_psip", 0.4, .3);
  setParam("sigma_chic2", 0.5, .7);
  setParam("sigma_chic1", 1.0, 1);
  setParam("sigma_jpsi", 1.5, 1);

  setParam("lth_psi", 0.4, 0.1, 0, 1);
  setParam("lth_c1", 0.11, 0.1);
  setParam("lth_c2", 0.5, 0.1);

  setParam("gamma", 0.65, 0.1);
  setParam("beta_long_psi", 3.5, 0.1);
  setParam("beta_total_psi", 3.3, 0.1);
  setParam("beta_long_c1", 3.5, 0.1);
  setParam("beta_total_c1", 3.35, 0.1);
  setParam("beta_long_c2", 3.5, 0.1);
  setParam("beta_total_c2", 3.35, 0.1);

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

  if (m_useCosthRatios) {
    setParam("norm_costh_1", 0.45, 0.1);
    setParam("norm_costh_2", 0.45, 0.1);
    setParam("norm_costh_3", 0.45, 0.1);
  }
}

void GlobalLikelihood::addNuissances()
{
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

template<typename T>
void GlobalLikelihood::addNuissanceParameter(const std::string& name, const T& par)
{
  const auto index = IPAR(name.c_str());
  m_nuissParams.push_back({index, {par}});
}


std::vector<TGraphAsymmErrors> GlobalLikelihood::getDataGraphs(const ROOT::Fit::FitResult& fitResult) const
{
  const auto parValues = fitResult.Parameters();
  std::vector<TGraphAsymmErrors> graphs;

  const double L_CMS = parValues[IPAR("L_CMS")];
  const double br_psip_mm = parValues[IPAR("br_psip_mm")];

  const auto psi2SXSecModel = getPsi2SXSecModel(parValues.data());
  const auto psiPolModel = getPsiPolModel(parValues.data());

  Identity<double> id;

  graphs.push_back(correctedCSGraph(m_psi2S_CMS_cs, {psi2SXSecModel}, {psiPolModel}, {id}, {1.0},
                                    L_CMS * br_psip_mm, M_PSI2S, "psi2S_CMS_cs"));


  const double L_ATLAS = parValues[IPAR("L_ATLAS")];
  const double br_psip_dp = parValues[IPAR("br_psip_dp")];
  const double br_jpsi_mm = parValues[IPAR("br_jpsi_mm")];

  graphs.push_back(correctedCSGraph(m_psi2S_ATLAS_cs, {psi2SXSecModel}, {psiPolModel}, {id}, {1.0},
                                    L_ATLAS * br_psip_dp * br_jpsi_mm, M_PSI2S, "psi2S_ATLAS_cs"));


  const double br_c1_jpsi = parValues[IPAR("br_c1_jpsi")];
  const double br_psip_c1 = parValues[IPAR("br_psip_c1")];

  const auto chi1XSecModel = getChi1XSecModel(parValues.data());
  const auto chi1PolModel = getChi1PolModel(parValues.data());

  const std::vector<CSModel> chi1CSModels = {chi1XSecModel, psi2SXSecModel};
  const std::vector<PolModel> chi1PolModels = {chi1PolModel, psiPolModel};
  const std::vector<PolFeedDownTrafo> chi1FDTrafos = {id, lambdaPsiToChi1};
  const std::vector<double> chi1FDFracs = {1.0, B_PSIP_CHIC1[0] / br_psip_c1};

  graphs.push_back(correctedCSGraph(m_chic1_ATLAS_cs, chi1CSModels, chi1PolModels, chi1FDTrafos, chi1FDFracs,
                                    L_ATLAS * br_c1_jpsi * br_jpsi_mm, M_CHIC1,
                                    "chic1_ATLAS_cs"));


  const double br_c2_jpsi = parValues[IPAR("br_c2_jpsi")];
  const double br_psip_c2 = parValues[IPAR("br_psip_c2")];

  const auto chi2XSecModel = getChi2XSecModel(parValues.data());
  const auto chi2PolModel = getChi2PolModel(parValues.data());

  const std::vector<CSModel> chi2CSModels = {chi2XSecModel, psi2SXSecModel};
  const std::vector<PolModel> chi2PolModels = {chi2PolModel, psiPolModel};
  const std::vector<PolFeedDownTrafo> chi2FDTrafos = {id, lambdaPsiToChi2};
  const std::vector<double> chi2FDFracs = {1.0, B_PSIP_CHIC2[0] / br_psip_c2};

  graphs.push_back(correctedCSGraph(m_chic2_ATLAS_cs, chi2CSModels, chi2PolModels, chi2FDTrafos, chi2FDFracs,
                                    L_ATLAS * br_c2_jpsi * br_jpsi_mm, M_CHIC2,
                                    "chic2_ATLAS_cs"));


  // chic ratios first have to compute the corrections individually

  const auto chi1Corrections = getCorrectionFactors(m_chic_ratio_CMS_cs,
                                                    chi1CSModels, chi1PolModels, chi1FDTrafos, chi1FDFracs,
                                                    L_ATLAS * br_c1_jpsi * br_jpsi_mm, M_JPSI,
                                                    m_clipCorrectionsPhysicalRange, {-1./3., 1.0});

  const auto chi2Corrections = getCorrectionFactors(m_chic_ratio_CMS_cs,
                                                    chi2CSModels, chi2PolModels, chi2FDTrafos, chi2FDFracs,
                                                    L_ATLAS * br_c2_jpsi * br_jpsi_mm, M_JPSI,
                                                    m_clipCorrectionsPhysicalRange, {-0.6, 1.0});

  std::vector<double> ratioCorrections;
  for (size_t i = 0; i < chi1Corrections.size(); ++i) {
    ratioCorrections.push_back(chi2Corrections[i] / chi1Corrections[i]);
  }

  graphs.push_back(correctGraph(m_chic_ratio_CMS_cs, "chic_ratio_CMS_cs", ratioCorrections));



  const double br_psip_jpsi = parValues[IPAR("br_psip_jpsi")];
  const auto jpsiXSecModel = getJpsiXSecModel(parValues.data());

  const std::vector<CSModel> jpsiCSModels = {jpsiXSecModel, chi1XSecModel, chi2XSecModel, psi2SXSecModel,
                                             psi2SXSecModel, psi2SXSecModel};
  const std::vector<PolModel> jpsiPolModels = {psiPolModel, chi1PolModel, chi2PolModel, psiPolModel,
                                               psiPolModel, psiPolModel};
  const std::vector<PolFeedDownTrafo> jpsiFDTrafos = {id, id, id, id, lambdaPsiToChi1, lambdaPsiToChi2};
  const std::vector<double> jpsiFDFracs = {1.0, B_CHIC1_JPSI[0] / br_c1_jpsi,
                                           B_CHIC2_JPSI[0] / br_c2_jpsi, B_PSIP_JPSI[0] / br_psip_jpsi,
                                           B_PSIP_CHIC1[0] * B_CHIC1_JPSI[0] / br_c1_jpsi / br_psip_c1,
                                           B_PSIP_CHIC2[0] * B_CHIC2_JPSI[0] / br_c2_jpsi / br_psip_c2};

  graphs.push_back(correctedCSGraph(m_jpsi_CMS_cs,
                                    jpsiCSModels, jpsiPolModels, jpsiFDTrafos, jpsiFDFracs,
                                    L_CMS * br_jpsi_mm, M_JPSI,
                                    "jpsi_CMS_cs"));


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


std::vector<TF1> GlobalLikelihood::getBestFitModels(const ROOT::Fit::FitResult& fitResult) const
{
  const auto parValues = fitResult.Parameters();
  std::vector<TF1> models;

  const auto psi2SXSecModel = getPsi2SXSecModel(parValues.data());
  const auto psiPolModel = getPsiPolModel(parValues.data());
  models.push_back(modelAsTF1(psi2SXSecModel, "psi_cs_direct", 2, 40));
  models.push_back(modelAsTF1(psiPolModel, "psi_pol_direct", 2, 40));

  const auto chi1XSecModel = getChi1XSecModel(parValues.data());
  const auto chi1PolModel = getChi1PolModel(parValues.data());
  models.push_back(modelAsTF1(chi1XSecModel, "chic1_cs_direct", 2, 40));
  models.push_back(modelAsTF1(chi1PolModel, "chic1_pol_direct", 2, 40));

  const auto chi2XSecModel = getChi2XSecModel(parValues.data());
  const auto chi2PolModel = getChi2PolModel(parValues.data());
  models.push_back(modelAsTF1(chi2XSecModel, "chic2_cs_direct", 2, 40));
  models.push_back(modelAsTF1(chi2PolModel, "chic2_pol_direct", 2, 40));

  const auto jpsiXSecModel = getJpsiXSecModel(parValues.data());
  models.push_back(modelAsTF1(jpsiXSecModel, "jpsi_cs_direct", 2, 40));

  CSModel chic1XSecFull, chic2XSecFull, jpsiXSecFull;
  PolModel chic1PolFull, chic2PolFull, jpsiPolFull;

  const Identity<double> id;

  const double br_psip_c1 = parValues[IPAR("br_psip_c1")];
  std::tie(chic1XSecFull, chic1PolFull) = combineModels({chi1XSecModel, psi2SXSecModel},
                                                        {chi1PolModel, psiPolModel},
                                                        {id, lambdaPsiToChi1},
                                                        {1.0, B_PSIP_CHIC1[0] / br_psip_c1});

  models.push_back(modelAsTF1(chic1XSecFull, "chic1_cs_full", 2, 40));
  models.push_back(modelAsTF1(chic1PolFull, "chic1_pol_full", 2, 40));

  const double br_psip_c2 = parValues[IPAR("br_psip_c2")];
  std::tie(chic2XSecFull, chic2PolFull) = combineModels({chi2XSecModel, psi2SXSecModel},
                                                        {chi2PolModel, psiPolModel},
                                                        {id, lambdaPsiToChi2},
                                                        {1.0, B_PSIP_CHIC2[0] / br_psip_c2});
  models.push_back(modelAsTF1(chic2XSecFull, "chic2_cs_full", 2, 40));
  models.push_back(modelAsTF1(chic2PolFull, "chic2_pol_full", 2, 40));

  // chic ratio cross section
  const CSModel chicRatioModel = [chic1XSecFull, chic2XSecFull] (double ptm) -> double {
    return chic2XSecFull(ptm) / chic1XSecFull(ptm);
  };
  models.push_back(modelAsTF1(chicRatioModel, "chic_ratio_cs_full", 2, 40));

  const double br_c1_jpsi = parValues[IPAR("br_c1_jpsi")];
  const double br_c2_jpsi = parValues[IPAR("br_c2_jpsi")];
  const double br_psip_jpsi = parValues[IPAR("br_psip_jpsi")];

  std::tie(jpsiXSecFull, jpsiPolFull) = combineModels({jpsiXSecModel, chic1XSecFull, chic2XSecFull,
        psi2SXSecModel},
    {psiPolModel, chic1PolFull, chic2PolFull, psiPolModel},
    {id, id, id, id},
    {1.0, B_CHIC1_JPSI[0] / br_c1_jpsi, B_CHIC2_JPSI[0] / br_c2_jpsi, B_PSIP_JPSI[0] / br_psip_jpsi});

  models.push_back(modelAsTF1(jpsiXSecFull, "jpsi_cs_full", 2, 40));
  models.push_back(modelAsTF1(jpsiPolFull, "jpsi_pol_full", 2, 40));


  if (m_useCosthRatios) {
    const std::array<double, 3> normalizations = {
      parValues[IPAR("norm_costh_1")],
      parValues[IPAR("norm_costh_2")],
      parValues[IPAR("norm_costh_3")]
    };

    // costh ratios CMS
    for (size_t iPtBin = 0; iPtBin < normalizations.size(); ++iPtBin) {
      double lambda1, lambda2;

      const auto& ratioData = m_chic_ratios_CMS_pol[iPtBin];

      std::tie(std::ignore, lambda1) = crossSecAndLambda(ratioData.first, 0.5 / M_JPSI,
                                                         {chi1XSecModel, psi2SXSecModel},
                                                         {chi1PolModel,psiPolModel},
                                                         {id, lambdaPsiToChi1},
                                                         {1.0, B_PSIP_CHIC1[0] / br_psip_c1});

      std::tie(std::ignore, lambda2) = crossSecAndLambda(ratioData.first, 0.5 / M_JPSI,
                                                         {chi2XSecModel, psi2SXSecModel},
                                                         {chi2PolModel, psiPolModel},
                                                         {id, lambdaPsiToChi2},
                                                         {1.0, B_PSIP_CHIC2[0] / br_psip_c2});

      const CosthRatioModel costhRatioModel = [lambda2, lambda1, normalizations, iPtBin] (double costh) {
        return costhRatio(costh, lambda2, lambda1, normalizations[iPtBin]);
      };

      std::stringstream gn;
      const int after_dec = (ratioData.first.ptM - int(ratioData.first.ptM)) * 100;
      gn << "chic_pol_ptM_" << int(ratioData.first.ptM) << "p" << after_dec;

      models.push_back(modelAsTF1(costhRatioModel, gn.str().c_str(), 0, 1));
    }
  }

  return models;
}


#endif
