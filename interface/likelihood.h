#ifndef H_LIKLIEHOOD__
#define H_LIKLIEHOOD__

#include "data_structures.h"
#include "constants.h"
#include "likelihood_helpers.h"
#include "simple_compile_time_map.h"

#include "Fit/ParameterSettings.h"

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
  size_t nPars() const { return m_startParams.size(); }

  int getParIdx(const std::string& name) const { return IPAR(name.c_str()); }

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

#endif
