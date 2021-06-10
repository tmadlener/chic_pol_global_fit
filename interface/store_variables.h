#ifndef STORE_VARIABLES_H__
#define STORE_VARIABLES_H__

#include "likelihood_helpers.h"

#ifndef NRQCD_FIT
#define NRQCD_FIT 0
#endif

#if NRQCD_FIT
#include "fit_parameters_nrqcd.h"
#else
#include "fit_parameters.h"
#endif

#include "TTree.h"

#include <tuple> // std::tie
#include <cmath> // std::isnan

template<typename LLH>
struct StoreVariables {
  StoreVariables(const LLH& llh) : m_llh(llh) {}

  TTree* create(bool storeParams=true);

  bool operator()(const double ptm, const std::vector<double>& p);

  bool valid() const;

  double ptM;
  double psip_cs;
  double jpsi_cs;
  double chic1_cs;
  double chic2_cs;
  double lth_psip;
  double lth_jpsi;
  double lth_chic1;
  double lth_chic2;

  double jpsi_cs_dir;
  double chic1_cs_dir;
  double chic2_cs_dir;
  double lth_chic1_dir;
  double lth_chic2_dir;
  double lth_jpsi_chic;

  double dlth_jpsi_psip;
  double dlth_chic2_chic1;

  double r_psip_chic2;
  double r_psip_chic1;
  double r_chic2_jpsi;
  double r_chic1_jpsi;
  double r_psip_jpsi;
  double r_chic_jpsi;

  double chic2_chic1_cs;
  double chic2_chic1_cs_br;

#if NRQCD_FIT
  double l_3PJ_8_jpsi;
  double l_3S1_8_jpsi;
  double l_3PJ_8_psip;
  double l_3S1_8_psip;
#endif

  std::array<double, NPARS()> pVals{};
  const LLH& m_llh;
};

template<typename LLH>
TTree* StoreVariables<LLH>::create(bool storeParams)
{
  TTree* tree = new TTree("ptm_dep_scan", "values of random scans as a function of pT/M");

  tree->Branch("ptM", &ptM);
  tree->Branch("psip_cs", &psip_cs);
  tree->Branch("jpsi_cs", &jpsi_cs);
  tree->Branch("chic1_cs", &chic1_cs);
  tree->Branch("chic2_cs", &chic2_cs);
  tree->Branch("lth_psip", &lth_psip);
  tree->Branch("lth_jpsi", &lth_jpsi);
  tree->Branch("lth_chic1", &lth_chic1);
  tree->Branch("lth_chic2", &lth_chic2);

  tree->Branch("jpsi_cs_dir", &jpsi_cs_dir);
  tree->Branch("chic1_cs_dir", &chic1_cs_dir);
  tree->Branch("chic2_cs_dir", &chic2_cs_dir);
  tree->Branch("lth_chic1_dir", &lth_chic1_dir);
  tree->Branch("lth_chic2_dir", &lth_chic2_dir);
  tree->Branch("lth_jpsi_chic", &lth_jpsi_chic);

  tree->Branch("dlth_jpsi_psip", &dlth_jpsi_psip);

  tree->Branch("r_psip_chic2", &r_psip_chic2);
  tree->Branch("r_psip_chic1", &r_psip_chic1);
  tree->Branch("r_chic2_jpsi", &r_chic2_jpsi);
  tree->Branch("r_chic1_jpsi", &r_chic1_jpsi);
  tree->Branch("r_psip_jpsi", &r_psip_jpsi);
  tree->Branch("r_chic_jpsi", &r_chic_jpsi);

  tree->Branch("chic2_chic1_cs", &chic2_chic1_cs);
  tree->Branch("chic2_chic1_cs_br", &chic2_chic1_cs_br);
  tree->Branch("dlth_chic2_chic1", &dlth_chic2_chic1);

#if NRQCD_FIT
  tree->Branch("l_3PJ_8_jpsi", &l_3PJ_8_jpsi);
  tree->Branch("l_3S1_8_jpsi", &l_3S1_8_jpsi);
  tree->Branch("l_3PJ_8_psip", &l_3PJ_8_psip);
  tree->Branch("l_3S1_8_psip", &l_3S1_8_psip);
#endif


  if (storeParams) {
    for (size_t iPar = 0; iPar < pVals.size(); ++iPar) {
      tree->Branch(PARNAME(iPar), &pVals[iPar]);
    }
  }

  return tree;
}

template<typename LLH>
bool StoreVariables<LLH>::operator()(const double ptm, const std::vector<double>& p)
{
  ptM = ptm;
  for (size_t i = 0; i < pVals.size(); ++i) {
    pVals[i] = p[i];
  }

  const auto psi2SXSecModel = m_llh.getPsi2SXSecModel(p.data());
  const auto chi2XSecModel = m_llh.getChi2XSecModel(p.data());
  const auto chi1XSecModel = m_llh.getChi1XSecModel(p.data());
  const auto jpsiXSecModel = m_llh.getJpsiXSecModel(p.data());

#if NRQCD_FIT
  const auto psi2SPolModel = m_llh.getPsi2SPolModel(p.data());
  const auto jpsiPolModel = m_llh.getJpsiPolModel(p.data());
#else
  const auto psi2SPolModel = m_llh.getPsiPolModel(p.data());
  const auto jpsiPolModel = m_llh.getPsiPolModel(p.data());
#endif
  const auto chi2PolModel = m_llh.getChi2PolModel(p.data());
  const auto chi1PolModel = m_llh.getChi1PolModel(p.data());

  jpsi_cs_dir = jpsiXSecModel(ptM);
  chic1_cs_dir = chi1XSecModel(ptM);
  chic2_cs_dir = chi2XSecModel(ptM);

  lth_chic1_dir = chi1PolModel(ptM);
  lth_chic2_dir = chi2PolModel(ptM);

  const Identity<double> id;

  std::tie(psip_cs, lth_psip) = crossSecAndLambda(ptM, {psi2SXSecModel}, {psi2SPolModel}, {id}, {1.0});

  const double br_psip_c2 = p[IPAR("br_psip_c2")];
  std::tie(chic2_cs, lth_chic2) = crossSecAndLambda(ptM, {chi2XSecModel, psi2SXSecModel},
                                                    {chi2PolModel, psi2SPolModel}, {id, lambdaPsiToChi2},
                                                    {1.0, B_PSIP_CHIC2[0] / br_psip_c2});

  const double br_psip_c1 = p[IPAR("br_psip_c1")];
  std::tie(chic1_cs, lth_chic1) = crossSecAndLambda(ptM, {chi1XSecModel, psi2SXSecModel},
                                                    {chi1PolModel, psi2SXSecModel}, {id, lambdaPsiToChi1},
                                                    {1.0, B_PSIP_CHIC1[0] / br_psip_c1});

  const double br_c2_jpsi = p[IPAR("br_c2_jpsi")];
  const double br_c1_jpsi = p[IPAR("br_c1_jpsi")];
  const double br_psip_jpsi = p[IPAR("br_psip_jpsi")];

  std::tie(jpsi_cs, lth_jpsi) = crossSecAndLambda(ptM,
                                                  {jpsiXSecModel, chi1XSecModel, chi2XSecModel,
                                                      psi2SXSecModel, psi2SXSecModel, psi2SXSecModel},
                                                  {jpsiPolModel, chi1PolModel, chi2PolModel,
                                                      psi2SPolModel, psi2SPolModel, psi2SPolModel},
                                                  {id, id, id, id, lambdaPsiToChi1, lambdaPsiToChi2},
                                                  {1.0, B_CHIC1_JPSI[0] / br_c1_jpsi,
                                                      B_CHIC2_JPSI[0] / br_c2_jpsi,
                                                      B_PSIP_JPSI[0] / br_psip_jpsi,
                                                      B_PSIP_CHIC1[0] * B_CHIC1_JPSI[0] / br_c1_jpsi / br_psip_c1,
                                                      B_PSIP_CHIC2[0] * B_CHIC2_JPSI[0] / br_c2_jpsi / br_psip_c2});

  // cross section ratio and delta lambda for the chic
  chic2_chic1_cs = chic2_cs / chic1_cs;
  chic2_chic1_cs_br = chic2_chic1_cs * B_CHIC2_JPSI[0] / B_CHIC1_JPSI[0] / br_c2_jpsi * br_c1_jpsi;
  dlth_chic2_chic1 = lth_chic2 - lth_chic1;


  // feed-down fractions
  r_psip_chic1 = psip_cs * B_PSIP_CHIC1[0] / br_psip_c1 / chic1_cs;
  r_psip_chic2 = psip_cs * B_PSIP_CHIC2[0] / br_psip_c2 / chic2_cs;

  r_chic1_jpsi = chic1_cs * B_CHIC1_JPSI[0] / br_c1_jpsi / jpsi_cs;
  r_chic2_jpsi = chic2_cs * B_CHIC2_JPSI[0] / br_c2_jpsi / jpsi_cs;
  r_psip_jpsi = psip_cs * B_PSIP_JPSI[0] / br_psip_jpsi / jpsi_cs;
  r_chic_jpsi = r_chic1_jpsi + r_chic2_jpsi;

  lth_jpsi_chic = weightedLambda({lth_chic1, lth_chic2}, {r_chic1_jpsi / r_chic_jpsi, r_chic2_jpsi / r_chic_jpsi});
  dlth_jpsi_psip = lth_jpsi - lth_psip;

#if NRQCD_FIT
  const double l_1S0_8_jpsi = p[IPAR("l_1S0_8_jpsi")];
  l_3PJ_8_jpsi = p[IPAR("l_r_3PJ_8_1S0_8_jpsi")] * l_1S0_8_jpsi;
  l_3S1_8_jpsi = p[IPAR("l_r_3S1_8_1S0_8_jpsi")] * l_1S0_8_jpsi;

  const double l_1S0_8_psip = p[IPAR("l_1S0_8_psip")];
  l_3PJ_8_psip = p[IPAR("l_rr_3PJ_8_1S0_8_psip_jpsi")] * p[IPAR("l_r_3PJ_8_1S0_8_jpsi")] * l_1S0_8_psip;
  l_3S1_8_psip = p[IPAR("l_rr_3S1_8_1S0_8_psip_jpsi")] * p[IPAR("l_r_3S1_8_1S0_8_jpsi")] * l_1S0_8_psip;
#endif

  return valid();
}

template<typename LLH>
bool StoreVariables<LLH>::valid() const
{
  return !(std::isnan(psip_cs) ||
           std::isnan(jpsi_cs) ||
           std::isnan(chic1_cs) ||
           std::isnan(chic2_cs) ||
           std::isnan(lth_psip) ||
           std::isnan(lth_jpsi) ||
           std::isnan(lth_chic1) ||
           std::isnan(lth_chic2) ||
           std::isnan(jpsi_cs_dir) ||
           std::isnan(chic1_cs_dir) ||
           std::isnan(chic2_cs_dir) ||
           std::isnan(lth_chic1_dir) ||
           std::isnan(lth_chic2_dir));
}


#endif
