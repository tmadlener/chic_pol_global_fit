#ifndef H_LIKLIEHOOD__
#define H_LIKLIEHOOD__

#include "data_structures.h"
#include "constants.h"
#include "likelihood_helpers.h"

#include "Fit/ParameterSettings.h"

#include <vector>
#include <array>
#include <utility>
#include <algorithm>
#include <iostream>
#include <functional>
#include <string>
#include <unordered_map>

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
  ParamsSettings getStartParams() const { return m_startParams; }

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
   * add all the nuissance parameters
   */
  void addNuissances();

  /**
   * Setup the partial contributions to the total likelihood
   */
  void setupContributions();

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

  ParamsSettings m_startParams;
  std::vector<std::pair<int, NuissanceParameter>> m_nuissParams;
  std::unordered_map<std::string, int> m_parNamesIdcs;
};


double GlobalLikelihood::operator()(const double* p) const
{
  double loglike = 0;

  CSModel psi2SXSecModel = [&p] (double ptm) {
    return sig_dir(ptm, p[0], p[4], p[7], p[6], p[5]);
  };
  auto psi2SPolModel = [&p] (double ptm) {
    return lambdathPsi(ptm, p[4], p[6], p[7], p[5]);
  };

  loglike += psi2SCrossSection(m_psi2S_CMS_cs, psi2SXSecModel, psi2SPolModel,
                               p[18], p[11]);

  loglike += psi2SCrossSection(m_psi2S_ATLAS_cs, psi2SXSecModel, psi2SPolModel,
                               p[19], p[10] * p[17]);

  CSModel chic2XSecModel = [&p] (double ptm) {
    return sig_dir(ptm, p[1], 1.0, p[9], p[6], p[5]);
  };
  ChiPolModel chic2PolModel = [&p, psi2SPolModel] (double ptm, double fdir) {
    return lambdathChic2(ptm, psi2SPolModel(ptm), fdir, 1.0, p[6], p[9], p[5]);
  };

  loglike += chicCrossSection(m_chic2_ATLAS_cs, psi2SXSecModel, chic2XSecModel, chic2PolModel,
                              p[19], p[12], p[15] * p[17], B_PSIP_CHIC2[0], M_CHIC2);


  CSModel chic1XSecModel = [&p] (double ptm) {
    return sig_dir(ptm, p[2], 1.0, p[8], p[6], p[5]);
  };
  ChiPolModel chic1PolModel = [&p, psi2SPolModel] (double ptm, double fdir) {
    return lambdathChic1(ptm, psi2SPolModel(ptm), fdir, 1.0, p[6], p[8], p[5]);
  };

  loglike += chicCrossSection(m_chic1_ATLAS_cs, psi2SXSecModel, chic1XSecModel, chic1PolModel,
                              p[19], p[13], p[16] * p[17], B_PSIP_CHIC1[0], M_CHIC1);


  CSModel jpsiXSecModel = [&p] (double ptm) {
    return sig_dir(ptm, p[3], p[4], p[7], p[6], p[5]);
  };
  JpsiPolModel jpsiPolModel = [&p, psi2SPolModel, chic1PolModel, chic2PolModel]
    (double ptm, double fdirPsiChi1, double fdirPsiChi2, double fpsi, double fchi1, double fchi2) {
    return lambdaJpsi(psi2SPolModel(ptm), chic1PolModel(ptm, fdirPsiChi1), chic2PolModel(ptm, fdirPsiChi2),
                      fpsi, fchi1, fchi2);
  };

  loglike += jpsiCrossSection(m_jpsi_CMS_cs, psi2SXSecModel, chic1XSecModel, chic2XSecModel, jpsiXSecModel, jpsiPolModel,
                              p[13], p[12], p[14], p[16], p[15], p[17], p[18]);


  loglike += chicCrossSectionRatio(m_chic_ratio_CMS_cs, psi2SXSecModel, chic1XSecModel, chic2XSecModel,
                                   chic1PolModel, chic2PolModel, p[13], p[12]);


  loglike += psi2SPolarization(m_psi2S_CMS_pol, psi2SPolModel);

  loglike += jpsiPolarization(m_jpsi_CMS_pol, jpsiPolModel, psi2SXSecModel, chic1XSecModel, chic2XSecModel, jpsiXSecModel,
                              p[13], p[12], p[14], p[16], p[15]);

  for (const auto& nuissPar : m_nuissParams) {
    loglike += nuissPar.second(p[nuissPar.first]);
  }

  return -2 * loglike;
}

void GlobalLikelihood::setupFit()
{
  defineStartParams();
  addNuissances();

}

void GlobalLikelihood::defineStartParams()
{
  m_startParams.emplace_back("sigma_psip", 24.81, 3);
  m_startParams.emplace_back("sigma_chic2", 84.25, 7);
  m_startParams.emplace_back("sigma_chic1", 104.48, 10);
  m_startParams.emplace_back("sigma_jpsi", 116.9, 10);

  m_startParams.emplace_back("f_ppsi", 0.02, 0.5, 0, 1);

  m_startParams.emplace_back("gamma", 0.74, 0.1);
  m_startParams.emplace_back("beta_u", 3.3, 0.01);
  m_startParams.emplace_back("beta_p", 2.7, 0.01);
  m_startParams.emplace_back("beta_c1", 3.3, 0.1);
  m_startParams.emplace_back("beta_c2", 3.3, 0.1);

  m_startParams.emplace_back("br_psip_dp", 1, 0.01);
  m_startParams.emplace_back("br_psip_mm", 1, 0.01);
  m_startParams.emplace_back("br_psip_c2", 1, 0.01);
  m_startParams.emplace_back("br_psip_c1", 1, 0.01);
  m_startParams.emplace_back("br_psip_jpsi", 1, 0.01);
  m_startParams.emplace_back("br_c2_jpsi", 1, 0.01);
  m_startParams.emplace_back("br_c1_jpsi", 1, 0.01);
  m_startParams.emplace_back("br_jpsi_mm", 1, 0.01);

  m_startParams.emplace_back("L_CMS", 1, 0.01);
  m_startParams.emplace_back("L_ATLAS", 1, 0.01);

  for (size_t i = 0; i < m_startParams.size(); ++i) {
    // m_parNamesIdcs
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

void GlobalLikelihood::setupContributions()
{
  // m_contributions.
}

template<typename T>
void GlobalLikelihood::addNuissanceParameter(const std::string& name, const T& par)
{
  const auto it = std::find_if(m_startParams.cbegin(), m_startParams.cend(),
                               [&name] (const ROOT::Fit::ParameterSettings& param) {
                                 return param.Name() == name;
                               });

  if (it != m_startParams.cend()) {
    auto index = std::distance(m_startParams.cbegin(), it);
    m_nuissParams.push_back({index, {par}});
  } else {
    std::cerr << "Problem when adding nuissance parameter: " << name << '\n';
  }
}

#endif
