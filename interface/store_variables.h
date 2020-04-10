#ifndef STORE_VARIABLES_H__
#define STORE_VARIABLES_H__

#include "likelihood.h"
#include "likelihood_helpers.h"

#include "TTree.h"

#include <tuple> // std::tie
#include <cmath> // std::isnan

struct StoreVariables {

  TTree* create();

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

  std::array<double, mapSize(PARAMETERS)> pVals{};
};

TTree* StoreVariables::create()
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


  for (size_t iPar = 0; iPar < pVals.size(); ++iPar) {
    tree->Branch(PARAMETERS[iPar].first, &pVals[iPar]);
  }

  return tree;
}

bool StoreVariables::operator()(const double ptm, const std::vector<double>& p)
{
  ptM = ptm;
  for (size_t i = 0; i < pVals.size(); ++i) {
    pVals[i] = p[i];
  }

  const auto psi2SXSecModel = GlobalLikelihood::getPsi2SXSecModel(p.data());
  const auto chi2XSecModel = GlobalLikelihood::getChi2XSecModel(p.data());
  const auto chi1XSecModel = GlobalLikelihood::getChi1XSecModel(p.data());
  const auto jpsiXSecModel = GlobalLikelihood::getJpsiXSecModel(p.data());

  const auto psiPolModel = GlobalLikelihood::getPsiPolModel(p.data());
  const auto chi2PolModel = GlobalLikelihood::getChi2PolModel(p.data());
  const auto chi1PolModel = GlobalLikelihood::getChi1PolModel(p.data());

  jpsi_cs_dir = jpsiXSecModel(ptM);
  chic1_cs_dir = chi1XSecModel(ptM);
  chic2_cs_dir = chi2XSecModel(ptM);

  lth_chic1_dir = chi1PolModel(ptM);
  lth_chic2_dir = chi2PolModel(ptM);

  const Identity<double> id;

  std::tie(psip_cs, lth_psip) = crossSecAndLambda(ptM, {psi2SXSecModel}, {psiPolModel}, {id}, {1.0});

  const double br_psip_c2 = p[IPAR("br_psip_c2")];
  std::tie(chic2_cs, lth_chic2) = crossSecAndLambda(ptM, {chi2XSecModel, psi2SXSecModel},
                                                    {chi2PolModel, psiPolModel}, {id, lambdaPsiToChi2},
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
                                                  {psiPolModel, chi1PolModel, chi2PolModel,
                                                      psiPolModel, psiPolModel, psiPolModel},
                                                  {id, id, id, id, lambdaPsiToChi1, lambdaPsiToChi2},
                                                  {1.0, B_CHIC1_JPSI[0] / br_c1_jpsi,
                                                      B_CHIC2_JPSI[0] / br_c2_jpsi,
                                                      B_PSIP_JPSI[0] / br_psip_jpsi,
                                                      B_PSIP_CHIC1[0] * B_CHIC1_JPSI[0] / br_c1_jpsi / br_psip_c1,
                                                      B_PSIP_CHIC2[0] * B_CHIC2_JPSI[0] / br_c2_jpsi / br_psip_c2});


  return valid();
}

bool StoreVariables::valid() const
{
  return !(  std::isnan(ptM) ||
             std::isnan(psip_cs) ||
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
