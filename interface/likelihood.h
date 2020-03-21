#ifndef H_LIKLIEHOOD__
#define H_LIKLIEHOOD__

#include "data_structures.h"
#include "constants.h"
#include "likelihood_helpers.h"
#include "simple_compile_time_map.h"

#include "Fit/ParameterSettings.h"
#include "Fit/FitResult.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"

#include <vector>
#include <array>
#include <utility>
#include <string>

constexpr std::array<ParameterIndex, 24> PARAMETERS = {{
    {"sigma_psip",  0},
    {"sigma_chic2", 1},
    {"sigma_chic1", 2},
    {"sigma_jpsi", 3},

    {"f_long_psi", 4},
    {"f_long_c1", 5},
    {"f_long_c2", 6},

    {"gamma", 7},
    {"beta_long_psi", 8},
    {"beta_trans_psi", 9},
    {"beta_long_c1", 10},
    {"beta_trans_c1", 11},
    {"beta_long_c2", 12},
    {"beta_trans_c2", 13},

    {"br_psip_dp", 14},
    {"br_psip_mm", 15},
    {"br_psip_c2", 16},
    {"br_psip_c1", 17},
    {"br_psip_jpsi", 18},
    {"br_c2_jpsi", 19},
    {"br_c1_jpsi", 20},
    {"br_jpsi_mm", 21},

    {"L_CMS", 22},
    {"L_ATLAS", 23}
  }};

/**
 * Convenience overload
 */
static constexpr int IPAR(const char* name)
{
  return getParIdx(PARAMETERS, name);
}


using ParamsSettings = std::vector<ROOT::Fit::ParameterSettings>;

class GlobalLikelihood {
public:
  GlobalLikelihood(const CrossSectionMeasurement& psi2S_ATLAS, const CrossSectionMeasurement& psi2S_CMS,
                   const CrossSectionMeasurement& chic2_ATLAS, const CrossSectionMeasurement& chic1_ATLAS,
                   const CrossSectionMeasurement& jpsi_CMS, const CrossSectionMeasurement& chic_ratio_CMS,
                   const PolarizationMeasurement& psi2S_CMS_p, const PolarizationMeasurement& jpsi_CMS_p) :
    m_psi2S_ATLAS_cs(psi2S_ATLAS), m_psi2S_CMS_cs(psi2S_CMS), m_chic2_ATLAS_cs(chic2_ATLAS),
    m_chic1_ATLAS_cs(chic1_ATLAS), m_jpsi_CMS_cs(jpsi_CMS), m_chic_ratio_CMS_cs(chic_ratio_CMS),
    m_psi2S_CMS_pol(psi2S_CMS_p), m_jpsi_CMS_pol(jpsi_CMS_p)
  {
    setupFit();
  }

  /**
   * Evaluate the likelihood with the given parameters
   */
  double operator()(const double* p) const;

  /**
   * Get the number of parameters
   */
  int nPars() const { return m_startParams.size(); }

  /**
   * Get the start parameter settings (necessary to pass them to the fitter)
   */
  const ParamsSettings& getStartParams() const { return m_startParams; }

  /**
   * Get the models (and sub-models) using the parameters of the fit result
   */
  std::vector<TF1> getPsi2SCSModel(const ROOT::Fit::FitResult& fitResult) const;
  std::vector<TF1> getChicCSModel(const ROOT::Fit::FitResult& fitResult) const;
  std::vector<TF1> getJpsiCSModel(const ROOT::Fit::FitResult& fitResult) const;

  /**
   * Get the data graphs (with applied corrections, etc.)
   */
  std::vector<TGraphAsymmErrors> getDataGraphs(const ROOT::Fit::FitResult& fitResult) const;

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

  ParamsSettings m_startParams{mapSize(PARAMETERS)}; // default initialize
  std::vector<std::pair<int, NuissanceParameter>> m_nuissParams;
};


double GlobalLikelihood::operator()(const double* p) const
{
  double loglike = 0;

  // pulling the different parameters out of the p pointer here to guarantee, that the call
  // to IPAR is evaluated at compile time
  const double sigma_psip = p[IPAR("sigma_psip")];
  const double f_long_psi = p[IPAR("f_long_psi")];
  const double beta_trans_psi = p[IPAR("beta_trans_psi")];
  const double beta_long_psi = p[IPAR("beta_long_psi")];
  const double gamma = p[IPAR("gamma")];
  const double L_CMS = p[IPAR("L_CMS")];
  const double L_ATLAS = p[IPAR("L_ATLAS")];
  const double br_psip_mm = p[IPAR("br_psip_mm")];
  const double br_psip_dp = p[IPAR("br_psip_dp")];
  const double br_jpsi_mm = p[IPAR("br_jpsi_mm")];

  CSModel psi2SXSecModel = [sigma_psip, f_long_psi, beta_trans_psi, beta_long_psi, gamma] (double ptm) {
    return sig_dir(ptm, sigma_psip, f_long_psi, beta_long_psi, beta_trans_psi, gamma);
  };
  PolModel psi2SPolModel = [f_long_psi, beta_trans_psi, beta_long_psi, gamma] (double ptm) {
    return lambdaTheta(ptm, f_long_psi, beta_long_psi, beta_trans_psi, gamma);
  };

  const Identity<double> id;

  // psi(2S) CMS
  loglike += loglikeCrossSection(m_psi2S_CMS_cs,
                                 {psi2SXSecModel}, {psi2SPolModel}, {id}, {1.0},
                                 L_CMS * br_psip_mm, M_PSI2S);

  // psi(2S) ATLAS
  loglike += loglikeCrossSection(m_psi2S_ATLAS_cs,
                                 {psi2SXSecModel}, {psi2SPolModel}, {id}, {1.0},
                                 L_ATLAS * br_psip_dp * br_jpsi_mm, M_PSI2S);

  const double sigma_chic2 = p[IPAR("sigma_chic2")];
  const double beta_long_c2 = p[IPAR("beta_long_c2")];
  const double beta_trans_c2 = p[IPAR("beta_trans_c2")];
  const double f_long_c2 = p[IPAR("f_long_c2")];
  const double br_c2_jpsi = p[IPAR("br_c2_jpsi")];
  const double br_psip_c2 = p[IPAR("br_psip_c2")];

  CSModel chic2XSecModel = [sigma_chic2, f_long_c2, beta_long_c2, beta_trans_c2, gamma] (double ptm) {
    return sig_dir(ptm, sigma_chic2, f_long_c2, beta_long_c2, beta_trans_c2, gamma);
  };
  PolModel chi2PolModel = [f_long_c2, beta_long_c2, beta_trans_c2, gamma] (double ptm) {
    return lambdaTheta(ptm, f_long_c2, beta_long_c2, beta_trans_c2, gamma);
  };

  loglike += loglikeCrossSection(m_chic2_ATLAS_cs,
                                 {chic2XSecModel, psi2SXSecModel}, {chi2PolModel, psi2SPolModel},
                                 {id, lambdaPsiToChi2}, {1.0, B_PSIP_CHIC2[0] / br_psip_c2},
                                 L_ATLAS * br_c2_jpsi * br_jpsi_mm, M_CHIC2);



  const double sigma_chic1 = p[IPAR("sigma_chic1")];
  const double beta_long_c1 = p[IPAR("beta_long_c1")];
  const double beta_trans_c1 = p[IPAR("beta_trans_c1")];
  const double f_long_c1 = p[IPAR("f_long_c1")];
  const double br_c1_jpsi = p[IPAR("br_c1_jpsi")];
  const double br_psip_c1 = p[IPAR("br_psip_c1")];

  CSModel chic1XSecModel = [sigma_chic1, f_long_c1, beta_long_c1, beta_trans_c1, gamma] (double ptm) {
    return sig_dir(ptm, sigma_chic1, f_long_c1, beta_long_c1, beta_trans_c1, gamma);
  };
  PolModel chi1PolModel = [f_long_c1, beta_trans_c1, beta_long_c1, gamma] (double ptm) {
    return lambdaTheta(ptm, f_long_c1, beta_long_c1, beta_trans_c1, gamma);
  };

  loglike += loglikeCrossSection(m_chic1_ATLAS_cs,
                                 {chic1XSecModel, psi2SXSecModel}, {chi1PolModel, psi2SPolModel},
                                 {id, lambdaPsiToChi1}, {1.0, B_PSIP_CHIC1[0] / br_psip_c1},
                                 L_ATLAS * br_c1_jpsi * br_jpsi_mm, M_CHIC1);


  const double sigma_jpsi = p[IPAR("sigma_jpsi")];
  const double br_psip_jpsi = p[IPAR("br_psip_jpsi")];

  CSModel jpsiXSecModel = [sigma_jpsi, f_long_psi, beta_trans_psi, beta_long_psi, gamma] (double ptm) {
    return sig_dir(ptm, sigma_jpsi, f_long_psi, beta_trans_psi, beta_long_psi, gamma);
  };

  // jpsi direct polarization is the same as psi(2S) direct polarization
  // TODO: make it possible to compose one function that reflects the chi
  // polarization dependening on pt/M. Currently, we it is necessary to be very
  // general here and list all the contributions to the j/psi individually,
  // considering double feed-down accordingly. NOTE: The lambdas are only
  // affected in the feed-down decay from psi(2S) -> chi

  loglike += loglikeCrossSection(m_jpsi_CMS_cs,
                                 {jpsiXSecModel, chic1XSecModel, chic2XSecModel, psi2SXSecModel,
                                     psi2SXSecModel, psi2SXSecModel},
                                 {psi2SPolModel, chi1PolModel, chi2PolModel, psi2SPolModel,
                                     psi2SPolModel, psi2SPolModel},
                                 {id, id, id, id, lambdaPsiToChi1, lambdaPsiToChi2},
                                 {1.0, B_CHIC1_JPSI[0] / br_c1_jpsi, B_CHIC2_JPSI[0] / br_c2_jpsi,
                                     B_PSIP_JPSI[0] / br_psip_jpsi,
                                     B_PSIP_CHIC1[0] * B_CHIC1_JPSI[0] / br_c1_jpsi / br_psip_c1,
                                     B_PSIP_CHIC2[0] * B_CHIC2_JPSI[0] / br_c2_jpsi / br_psip_c2},
                                 L_CMS * br_jpsi_mm, M_JPSI);


  loglike += loglikeCrossSectionRatio(m_chic_ratio_CMS_cs,
                                      {chic2XSecModel, psi2SXSecModel}, {chic1XSecModel, psi2SXSecModel},
                                      {chi2PolModel, psi2SPolModel}, {chi1PolModel, psi2SPolModel},
                                      {id, lambdaPsiToChi2}, {id, lambdaPsiToChi1},
                                      {1.0, B_PSIP_CHIC2[0] / br_psip_c2},
                                      {1.0, B_PSIP_CHIC1[0] / br_psip_c1},
                                      br_psip_c2 / br_psip_c1, M_JPSI);



  loglike += loglikePolarization(m_psi2S_CMS_pol, {psi2SXSecModel}, {psi2SPolModel}, {id}, {1.0}, M_PSI2S);

  // same as above for the cross-sections. All contributions have to be considered individually
  loglike += loglikePolarization(m_jpsi_CMS_pol,
                                 {jpsiXSecModel, chic1XSecModel, chic2XSecModel, psi2SXSecModel,
                                     psi2SXSecModel, psi2SXSecModel},
                                 {psi2SPolModel, chi1PolModel, chi2PolModel, psi2SPolModel,
                                     psi2SPolModel, psi2SPolModel},
                                 {id, id, id, id, lambdaPsiToChi1, lambdaPsiToChi2},
                                 {1.0, B_CHIC1_JPSI[0] / br_c1_jpsi, B_CHIC2_JPSI[0] / br_c2_jpsi,
                                     B_PSIP_JPSI[0] / br_psip_jpsi,
                                     B_PSIP_CHIC1[0] * B_CHIC1_JPSI[0] / br_c1_jpsi / br_psip_c1,
                                     B_PSIP_CHIC2[0] * B_CHIC2_JPSI[0] / br_c2_jpsi / br_psip_c2},
                                 M_JPSI);

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
  setParam("sigma_psip", 4.0, 3);
  setParam("sigma_chic2", 14.0, 7);
  setParam("sigma_chic1", 17.5, 10);
  setParam("sigma_jpsi", 19.8, 10);

  setParam("f_long_psi", 0.3, 0.1, 0, 1);
  setParam("f_long_c1", 0.11, 0.1, 0, 1);
  setParam("f_long_c2", 0.5, 0.1, 0, 1);

  setParam("gamma", 0.74, 0.1);
  setParam("beta_long_psi", 3.3, 0.1);
  setParam("beta_trans_psi", 3.3, 0.1);
  setParam("beta_long_c1", 3.3, 0.1);
  setParam("beta_trans_c1", 3.3, 0.1);
  setParam("beta_long_c2", 3.3, 0.1);
  setParam("beta_trans_c2", 3.3, 0.1);

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


std::vector<TF1> GlobalLikelihood::getPsi2SCSModel(const ROOT::Fit::FitResult& fitResult) const
{
  const auto& parValues = fitResult.Parameters();
  const double sigma = parValues[IPAR("sigma_psip")];
  const double fLong = parValues[IPAR("f_long_psi")];
  const double betaLong = parValues[IPAR("beta_long_psi")];
  const double betaTrans = parValues[IPAR("beta_trans_psi")];
  const double gamma = parValues[IPAR("gamma")];

  return psi2SCSModel(sigma, fLong, betaLong, betaTrans, gamma);
}

std::vector<TF1> GlobalLikelihood::getChicCSModel(const ROOT::Fit::FitResult& fitResult) const
{
  const auto& parValues = fitResult.Parameters();

  // psi2S
  const double sigmaPsi = parValues[IPAR("sigma_psip")];
  const double fLongPsi = parValues[IPAR("f_long_psi")];
  const double betaLongPsi = parValues[IPAR("beta_long_psi")];
  const double betaTransPsi = parValues[IPAR("beta_trans_psi")];
  const double gamma = parValues[IPAR("gamma")];

  // chic2
  const double sigma2 = parValues[IPAR("sigma_chic2")];
  const double fLong2 = parValues[IPAR("f_long_c2")];
  const double betaLong2 = parValues[IPAR("beta_long_c2")];
  const double betaTrans2 = parValues[IPAR("beta_trans_c2")];
  const double brPsiPChi2 = B_PSIP_CHIC2[0] / parValues[IPAR("br_psip_c2")];

  auto models = chicCSModel("chic2", sigma2, fLong2, betaLong2, betaTrans2, gamma,
                            sigmaPsi, fLongPsi, betaLongPsi, betaTransPsi, brPsiPChi2);

  // chic1
  const double sigma1 = parValues[IPAR("sigma_chic1")];
  const double fLong1 = parValues[IPAR("f_long_c1")];
  const double betaLong1 = parValues[IPAR("beta_long_c1")];
  const double betaTrans1 = parValues[IPAR("beta_trans_c1")];
  const double brPsiPChi1 =  B_PSIP_CHIC1[0] / parValues[IPAR("br_psip_c1")];

  auto chi1Models = chicCSModel("chic1", sigma1, fLong1, betaLong1, betaTrans1, gamma,
                                sigmaPsi, fLongPsi, betaLongPsi, betaTransPsi, brPsiPChi1);

  for (const auto& m : chi1Models) {
    models.push_back(m);
  }

  return models;
}

std::vector<TF1> GlobalLikelihood::getJpsiCSModel(const ROOT::Fit::FitResult& fitResult) const
{
  const auto& parValues = fitResult.Parameters();

  // psi2S
  const double sigmaPsi = parValues[IPAR("sigma_psip")];
  const double fLongPsi = parValues[IPAR("f_long_psi")];
  const double betaLongPsi = parValues[IPAR("beta_long_psi")];
  const double betaTransPsi = parValues[IPAR("beta_trans_psi")];
  const double gamma = parValues[IPAR("gamma")];

  // chic2
  const double sigma2 = parValues[IPAR("sigma_chic2")];
  const double fLong2 = parValues[IPAR("f_long_c2")];
  const double betaLong2 = parValues[IPAR("beta_long_c2")];
  const double betaTrans2 = parValues[IPAR("beta_trans_c2")];
  const double brPsiPChi2 = B_PSIP_CHIC2[0] / parValues[IPAR("br_psip_c2")];
  const double brChi2Jpsi = B_CHIC2_JPSI[0] / parValues[(IPAR("br_c2_jpsi"))];

  // chic1
  const double sigma1 = parValues[IPAR("sigma_chic1")];
  const double fLong1 = parValues[IPAR("f_long_c1")];
  const double betaLong1 = parValues[IPAR("beta_long_c1")];
  const double betaTrans1 = parValues[IPAR("beta_trans_c1")];
  const double brPsiPChi1 =  B_PSIP_CHIC1[0] / parValues[IPAR("br_psip_c1")];
  const double brChi1Jpsi = B_CHIC1_JPSI[0] / parValues[IPAR("br_c1_jpsi")];

  // jpsi
  const double sigmaJ = parValues[IPAR("sigma_jpsi")];
  const double brPsiJpsi = B_PSIP_JPSI[0] / parValues[IPAR("br_psip_jpsi")];

  return jpsiCSModel(sigmaJ, fLongPsi, betaLongPsi, betaTransPsi, gamma,
                     sigma1, fLong1, betaLong1, betaTrans1, sigma2, fLong2, betaLong2, betaTrans2,
                     sigmaPsi, brPsiPChi1, brPsiPChi2, brChi1Jpsi, brChi2Jpsi, brPsiJpsi);

}

std::vector<TGraphAsymmErrors> GlobalLikelihood::getDataGraphs(const ROOT::Fit::FitResult& fitResult) const
{
  const auto& parValues = fitResult.Parameters();
  std::vector<TGraphAsymmErrors> graphs;

  const double f_long_psi = parValues[IPAR("f_long_psi")];
  const double beta_long_psi = parValues[IPAR("beta_long_psi")];
  const double beta_trans_psi = parValues[IPAR("beta_trans_psi")];
  const double gamma = parValues[IPAR("gamma")];
  const double br_psip_mm = parValues[IPAR("br_psip_mm")];
  const double L_CMS = parValues[IPAR("L_CMS")];

  auto psi2SPolModel = [f_long_psi, beta_trans_psi, beta_long_psi, gamma] (double ptm) {
    return lambdaTheta(ptm, f_long_psi, beta_long_psi, beta_trans_psi, gamma);
  };

  graphs.push_back(correctCSData(m_psi2S_CMS_cs, "psi2S_CMS_cs",
                                 psi2SCorrections(m_psi2S_CMS_cs, psi2SPolModel, br_psip_mm * L_CMS)));

  const double L_ATLAS = parValues[IPAR("L_ATLAS")];
  const double br_psip_jpsi = parValues[IPAR("br_psip_jpsi")];
  const double br_psip_dp = parValues[IPAR("br_psip_dp")];
  graphs.push_back(correctCSData(m_psi2S_ATLAS_cs, "psi2S_ATLAS_cs",
                                 psi2SCorrections(m_psi2S_ATLAS_cs, psi2SPolModel,
                                                  br_psip_dp * br_psip_jpsi * L_ATLAS)));

  const double sigma_psip = parValues[IPAR("sigma_psip")];

  CSModel psi2SXSecModel = [sigma_psip, f_long_psi, beta_trans_psi, beta_long_psi, gamma] (double ptm) {
    return sig_dir(ptm, sigma_psip, f_long_psi, beta_long_psi, beta_trans_psi, gamma);
  };

  const double sigma_chic2 = parValues[IPAR("sigma_chic2")];
  const double beta_long_c2 = parValues[IPAR("beta_long_c2")];
  const double beta_trans_c2 = parValues[IPAR("beta_trans_c2")];
  const double f_long_c2 = parValues[IPAR("f_long_c2")];
  const double br_c2_jpsi = parValues[IPAR("br_c2_jpsi")];
  const double br_psip_c2 = parValues[IPAR("br_psip_c2")];
  const double br_jpsi_mm = parValues[IPAR("br_jpsi_mm")];

  CSModel chic2XSecModel = [sigma_chic2, f_long_c2, beta_long_c2, beta_trans_c2, gamma] (double ptm) {
    return sig_dir(ptm, sigma_chic2, f_long_c2, beta_long_c2, beta_trans_c2, gamma);
  };
  auto chic2PolModel = [f_long_c2, beta_long_c2, beta_trans_c2, gamma, psi2SPolModel]
    (double ptm, double fdir) {
    return lambdaChic(ptm, f_long_c2, beta_long_c2, beta_trans_c2, gamma, psi2SPolModel(ptm), fdir);
  };

  graphs.push_back(correctCSData(m_chic2_ATLAS_cs, "chic2_ATLAS_cs",
                                 chicCorrections(m_chic2_ATLAS_cs, psi2SXSecModel, chic2XSecModel,
                                                 chic2PolModel, B_PSIP_CHIC2[0], M_CHIC2, br_psip_c2,
                                                 L_ATLAS * br_c2_jpsi * br_jpsi_mm)));

  const double sigma_chic1 = parValues[IPAR("sigma_chic1")];
  const double beta_long_c1 = parValues[IPAR("beta_long_c1")];
  const double beta_trans_c1 = parValues[IPAR("beta_trans_c1")];
  const double f_long_c1 = parValues[IPAR("f_long_c1")];
  const double br_c1_jpsi = parValues[IPAR("br_c1_jpsi")];
  const double br_psip_c1 = parValues[IPAR("br_psip_c1")];

  CSModel chic1XSecModel = [sigma_chic1, f_long_c1, beta_long_c1, beta_trans_c1, gamma] (double ptm) {
    return sig_dir(ptm, sigma_chic1, f_long_c1, beta_long_c1, beta_trans_c1, gamma);
  };
  auto chic1PolModel = [f_long_c1, beta_long_c1, beta_trans_c1, gamma, psi2SPolModel]
    (double ptm, double fdir) {
    return lambdaChic(ptm, f_long_c1, beta_long_c1, beta_trans_c1, gamma, psi2SPolModel(ptm), fdir);
  };

  graphs.push_back(correctCSData(m_chic1_ATLAS_cs, "chic1_ATLAS_cs",
                                 chicCorrections(m_chic1_ATLAS_cs, psi2SXSecModel, chic1XSecModel,
                                                 chic1PolModel, B_PSIP_CHIC1[0], M_CHIC1, br_psip_c1,
                                                 L_ATLAS * br_c1_jpsi * br_jpsi_mm)));

  const double sigma_jpsi = parValues[IPAR("sigma_jpsi")];

  CSModel jpsiXSecModel = [sigma_jpsi, f_long_psi, beta_trans_psi, beta_long_psi, gamma] (double ptm) {
    return sig_dir(ptm, sigma_jpsi, f_long_psi, beta_trans_psi, beta_long_psi, gamma);
  };
  auto jpsiPolModel = [psi2SPolModel, chic1PolModel, chic2PolModel]
    (double ptm, double fdirPsiChi1, double fdirPsiChi2, double fpsi, double fchi1, double fchi2) {
    return lambdaJpsi(psi2SPolModel(ptm), chic1PolModel(ptm, fdirPsiChi1), chic2PolModel(ptm, fdirPsiChi2),
                      fpsi, fchi1, fchi2);
  };

  graphs.push_back(correctCSData(m_jpsi_CMS_cs, "jpsi_CMS_cs",
                                 jpsiCorrections(m_jpsi_CMS_cs, psi2SXSecModel, chic1XSecModel,
                                                 chic2XSecModel, jpsiXSecModel, jpsiPolModel, br_psip_c1,
                                                 br_psip_c2, br_c1_jpsi, br_c2_jpsi, br_psip_jpsi,
                                                 L_CMS * br_jpsi_mm)));

  return graphs;
}

#endif
