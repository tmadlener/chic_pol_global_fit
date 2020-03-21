#ifndef H_LIKELIHOOD_HELPERS__
#define H_LIKELIHOOD_HELPERS__

#include "data_structures.h"
#include "constants.h"

#include "TF1.h"

#include <cmath>
#include <functional>
#include <algorithm>
#include <utility>

using CSModel = std::function<double(double)>;
using PolModel = std::function<double(double)>;
using PolFeedDownTrafo = std::function<double(double)>;

double power_law(double x, double beta, double gamma)
{
  return x * std::pow(1 + x * x / (beta - 2) / gamma, -beta);
}

double norm_plaw(double x, double beta, double gamma, double xnorm=PTMNORM)
{
  return power_law(x, beta, gamma) / power_law(xnorm, beta, gamma);
}

double sig_dir(double ptm, double sig, double fLong, double betaLong, double betaTrans, double gamma)
{
  return sig * ((1 - fLong) * norm_plaw(ptm, betaTrans, gamma) +
                fLong * norm_plaw(ptm, betaLong, gamma));
}

/**
 * Calculate the weighted lambda from an arbitrary number of lambdas and weights
 */
double weightedLambda(const std::vector<double>& lambdas, const std::vector<double>& weights)
{
  double num = 0;
  double denom = 0;

  for (size_t i = 0; i < lambdas.size(); ++i) {
    const double lambda = lambdas[i];
    const double weight = weights[i];

    num += weight * lambda / (3 + lambda);
    denom += weight / (3 + lambda);
  }

  return num / denom;
}

/**
 * Convenience wrapper for two lambdas and weights
 */
double weightedLambda(const double lambda1, const double lambda2, const double w1, const double w2)
{
  return weightedLambda({lambda1, lambda2}, {w1, w2});
}

/**
 * lambda theta assuming a fully transverse and a fully longitudinal component
 * given the fraction of the longitudinal component
 */
double lambdath(const double fLong)
{
  return (1 - 3 * fLong) / (1 + fLong);
}

/**
 * ptm dependent lambda theta from the shape parameters of the two components and the longitudinal
 * fraction at a given reference point
 */
double lambdaTheta(double ptm, double fLong, double betaLong, double betaTrans, double gamma)
{
  const double contLong = fLong * norm_plaw(ptm, betaLong, gamma);
  const double contTrans = (1 - fLong) * norm_plaw(ptm, betaTrans, gamma);

  return lambdath(contLong / (contLong + contTrans));
}

double lambdaChic(double ptm, double fLong, double betaLong, double betaTrans, double gamma, double lambdaPsi, double fDir)
{
  const double lambdaChiDir = lambdaTheta(ptm, fLong, betaLong, betaTrans, gamma);
  return weightedLambda(lambdaChiDir, lambdaPsi / (4 + lambdaPsi), fDir, 1 - fDir);
}

double lambdaJpsi(double lambdaPsi, double lambdaChic1, double lambdaChic2,
                  double fdir, double fchic1, double fchic2)
{
  const double fpsi = 1 - (fdir + fchic1 + fchic2);
  return weightedLambda({lambdaPsi, lambdaPsi, lambdaChic1, lambdaChic2},
                        {fdir, fpsi, fchic1, fchic2});
}


/**
 * Calculate the polarization dependent correction for the measured cross section
 */
double acceptanceCorrection(const double kFactor, const double lambda)
{
  return (1 + 1. / 3. * kFactor) / (1 + (1 - lambda) / (3 + lambda) * kFactor);
}

/**
 * Sample the passed function in the passed range using the provided step size.
 */
template<typename Func>
double sampleFunc(const double min, const double max, const double step, Func func)
{
  double val = 0;
  int nsteps = 0;
  for (double x = min; x < max; x += step) {
    val += func(x);
    nsteps++;
  }
  return val / nsteps;
}


/**
 * Calculate all the partial cross sections at the given DataPoint from the
 * different cross section models that contribute with given branching
 * fractions.
 */
template<typename DataType>
std::vector<double> partialCrossSections(const DataType& point, double step,
                                         const std::vector<CSModel>& csModels,
                                         const std::vector<double>& brFracs)
{
  std::vector<double> partCS;
  for (size_t i = 0; i < csModels.size(); ++i) {


    partCS.push_back(sampleFunc(point.low, point.high + 1e-5, step, csModels[i]) * brFracs[i]);
  }
  return partCS;
}

/**
 * Calculate the lambda of a state that is affected by feed-down contributions.
 * Inputs are the models that describe the lambdas of the different
 * contributions and how they are affected in the feed-down process.
 * Additionally the fraction of each contribution is necessary.
 */
double lambdaFeedDown(double ptm, const std::vector<PolModel>& polModels,
                      const std::vector<PolFeedDownTrafo>& feedDownTrafos,
                      const std::vector<double>& weights)
{
  std::vector<double> lambdas;
  for (size_t i = 0; i < polModels.size(); ++i) {
    const double lam = polModels[i](ptm);
    lambdas.push_back(feedDownTrafos[i](lam));
  }
  return weightedLambda(lambdas, weights);
}

/**
 * Calculate the predicted cross section and lambda parameter for a given data
 * point, given the cross section and polarization models as well as the
 * functions describing how the lambdas transform in each feed-down decay and
 * the corresponding branching fractions for each process
 */
template<typename DataType>
std::pair<double, double> crossSecAndLambda(const DataType& point, double step,
                                            const std::vector<CSModel>& csModels,
                                            const std::vector<PolModel>& polModels,
                                            const std::vector<PolFeedDownTrafo>& feedDownTrafos,
                                            const std::vector<double>& brFracs)
{
  const auto partialCrossSecs = partialCrossSections(point, step, csModels, brFracs);

  const double totalCS = std::accumulate(partialCrossSecs.cbegin(), partialCrossSecs.cend(), 0.0);
  std::vector<double> polWeights;
  for (const auto& ps : partialCrossSecs) { polWeights.push_back(ps / totalCS); }

  const double lambda = lambdaFeedDown(point.ptM, polModels, feedDownTrafos, polWeights);
  return {totalCS, lambda};
}

/**
 * Calculate the contribution to the log-likelihood given a set of cross-section
 * measurements and models describing the cross section and polarization of each
 * state contributing to the measured state.
 */
double loglikeCrossSection(const CrossSectionMeasurement& data,
                           const std::vector<CSModel>& csModels,
                           const std::vector<PolModel>& polModels,
                           const std::vector<PolFeedDownTrafo>& feedDownTrafos,
                           const std::vector<double>& brFracs,
                           const double globNuiss,
                           const double mass)
{
  double loglike = 0;
  for (const auto& point : data) {
    double cs, lambda;
    std::tie(cs, lambda) = crossSecAndLambda(point, 0.5 / mass,
                                             csModels, polModels, feedDownTrafos, brFracs);

    const double accCorr = acceptanceCorrection(point.K, lambda);
    const double crossSec = point.xSec * accCorr * globNuiss;
    const double uncer = point.uncer * accCorr * globNuiss;

    const double relDiff = (cs - crossSec) / uncer;
    loglike += -0.5 * relDiff * relDiff;
  }

  return loglike;
}

/**
 * Calculate the contribution to the log-likelihood for a cross-section ratio
 * measurement
 */
double loglikeCrossSectionRatio(const CrossSectionMeasurement& data,
                                const std::vector<CSModel>& csModelsN,
                                const std::vector<CSModel>& csModelsD,
                                const std::vector<PolModel>& polModelsN,
                                const std::vector<PolModel>& polModelsD,
                                const std::vector<PolFeedDownTrafo>& fdTrafosN,
                                const std::vector<PolFeedDownTrafo>& fdTrafosD,
                                const std::vector<double>& brFracsN,
                                const std::vector<double>& brFracsD,
                                const double globNuiss, const double mass)
{
  double loglike = 0;
  for (const auto& point : data) {
    double csN, lambdaN, csD, lambdaD;
    std::tie(csN, lambdaN) = crossSecAndLambda(point, 0.5 / mass,
                                               csModelsN, polModelsN, fdTrafosN, brFracsN);
    std::tie(csD, lambdaD) = crossSecAndLambda(point, 0.5 / mass,
                                               csModelsD, polModelsD, fdTrafosD, brFracsD);

    const double accCorr = acceptanceCorrection(point.K, lambdaN) / acceptanceCorrection(point.K, lambdaD);
    const double crossSecRatio = point.xSec * accCorr * globNuiss;
    const double uncer = point.uncer * accCorr * globNuiss;

    const double relDiff = (csN / csD - crossSecRatio) / uncer;
    loglike += -0.5 * relDiff * relDiff;
  }

  return loglike;
}

/**
 * Calculate the contribution to the log-likelihood for a polarization
 * measurement. The cross section models are necessary to calculate the
 * contributions to the measured state.
 */
double loglikePolarization(const PolarizationMeasurement& data,
                           const std::vector<CSModel>& csModels,
                           const std::vector<PolModel>& polModels,
                           const std::vector<PolFeedDownTrafo>& fdTrafos,
                           const std::vector<double>& brFracs,
                           const double mass)
{
  double loglike = 0;
  for (const auto& point : data) {
    double lambda;
    std::tie(std::ignore, lambda) = crossSecAndLambda(point, 0.5 / mass,
                                                      csModels, polModels, fdTrafos, brFracs);

    const double relDiff = (lambda - point.lth) / point.uncer;
    loglike += -0.5 * relDiff * relDiff;
  }

  return loglike;
}


/**
 * Helper functor representing the identity function. In this case necessary
 * because in many feed-down decays the polarization is unaffected
 */
template<typename T>
struct Identity {
  T operator()(T x) const {return x;}
};

/**
 * Transformation of the psi(2S) lambda in the psi(2S) -> chi1 feed-down decay
 */
double lambdaPsiToChi1(const double lambdaPsi) {
  return lambdaPsi / (4 + lambdaPsi);
}

/**
 * Transformation of the psi(2S) lambda in the psi(2S) -> chi2 feed-down decay
 */
double lambdaPsiToChi2(const double lambdaPsi) {
  return 21 * lambdaPsi / (60 + 13 * lambdaPsi);
}





std::vector<TF1> psi2SCSModel(double sigma, double fLong, double betaLong, double betaTrans, double gamma)
{
  std::vector<TF1> models;

  TF1 fullModel("psi2S_cs_model", [](double* x, double* p)
                { return sig_dir(x[0], p[0], p[1], p[2], p[3], p[4]); },
                MIN_PTM, 40, 5);
  fullModel.SetParameters(sigma, fLong, betaLong, betaTrans, gamma);
  models.push_back(fullModel);

  TF1 longContr("psi2S_cs_long", [](double* x, double* p)
                { return p[0] * p[1] * norm_plaw(x[0], p[2], p[3]); },
                MIN_PTM, 40, 4);
  longContr.SetParameters(sigma, fLong, betaLong, gamma);
  models.push_back(longContr);

  TF1 transContr("psi2S_cs_trans", [](double* x, double* p)
                 { return p[0] * p[1] * norm_plaw(x[0], p[2], p[3]); },
                 MIN_PTM, 40, 4);
  transContr.SetParameters(sigma, 1 - fLong, betaTrans, gamma);
  models.push_back(transContr);

  return models;
}

std::vector<TF1> chicCSModel(const std::string& state,
                             double sigmaChi, double fLongChi, double betaLongChi, double betaTransChi, double gamma,
                             double sigmaPsi, double fLongPsi, double betaLongPsi, double betaTransPsi,
                             double brPsiChi)
{
  std::vector<TF1> models;
  TF1 directModel((state + "_cs_direct").c_str(),
                  [] (double* x, double* p) {return sig_dir(x[0], p[0], p[1], p[2], p[3], p[4]); },
                  MIN_PTM, 40, 5);
  directModel.SetParameters(sigmaChi, fLongChi, betaLongChi, betaTransChi, gamma);
  models.push_back(directModel);

  TF1 fullModel((state + "_cs_model").c_str(),
                [](double* x, double* p) {
                  return sig_dir(x[0], p[0], p[1], p[2], p[3], p[4]) +
                    p[5] * sig_dir(x[0], p[6], p[7], p[8], p[9], p[4]); },
                MIN_PTM, 40, 10);
  fullModel.SetParameters(sigmaChi, fLongChi, betaLongChi, betaTransChi, gamma, brPsiChi,
                          sigmaPsi, fLongPsi, betaLongPsi, betaTransPsi);
  models.push_back(fullModel);

  return models;
}

std::vector<TF1> jpsiCSModel(double sigmaJ, double fLongPsi, double betaLongPsi, double betaTransPsi,
                             double gamma, double sigma1, double fLong1, double betaLong1, double betaTrans1,
                             double sigma2, double fLong2, double betaLong2, double betaTrans2,
                             double sigmaPsi, double brPsiChi1, double brPsiChi2, double brChi1Jpsi,
                             double brChi2Jpsi, double brPsiJpsi)
{
  std::vector<TF1> models;

  TF1 directModel("jpsi_cs_direct",
                  [] (double* x, double* p) { return sig_dir(x[0], p[0], p[1], p[2], p[3], p[4]); },
                  MIN_PTM, 40, 5);
  directModel.SetParameters(sigmaJ, fLongPsi, betaLongPsi, betaTransPsi, gamma);
  models.push_back(directModel);

  TF1 longContr("jpsi_cs_direct_long", [](double* x, double* p)
                { return p[0] * p[1] * norm_plaw(x[0], p[2], p[3]); },
                MIN_PTM, 40, 4);
  longContr.SetParameters(sigmaJ, fLongPsi, betaLongPsi, gamma);
  models.push_back(longContr);

  TF1 transContr("jpsi_cs_direct_trans", [](double* x, double* p)
                 { return p[0] * p[1] * norm_plaw(x[0], p[2], p[3]); },
                 MIN_PTM, 40, 4);
  transContr.SetParameters(sigmaJ, 1 - fLongPsi, betaTransPsi, gamma);
  models.push_back(transContr);


  TF1 fullModel("jpsi_cs_model",
                [] (double* x, double* p) {
                  const double dir = sig_dir(x[0], p[0], p[1], p[2], p[3], p[4]);
                  const double dirChi1 = sig_dir(x[0], p[5], p[6], p[7], p[8], p[4]);
                  const double dirChi2 = sig_dir(x[0], p[9], p[10], p[11], p[12], p[4]);
                  const double dirPsi = sig_dir(x[0], p[13], p[1], p[2], p[3], p[4]);

                  const double chi1 = dirChi1 + dirPsi * p[14];
                  const double chi2 = dirChi2 + dirPsi * p[15];

                  const double fdChi1 = chi1 * p[16];
                  const double fdChi2 = chi2 * p[17];

                  const double jpsi = dir + fdChi1 + fdChi2 + dirPsi * p[18];

                  return jpsi;
                }, MIN_PTM, 40, 19);

  fullModel.SetParameters(sigmaJ, fLongPsi, betaLongPsi, betaTransPsi, gamma,
                          sigma1, fLong1, betaLong1, betaTrans1,
                          sigma2, fLong2);
  fullModel.SetParameter(11, betaLong2);
  fullModel.SetParameter(12, betaTrans2);
  fullModel.SetParameter(13, sigmaPsi);
  fullModel.SetParameter(14, brPsiChi1);
  fullModel.SetParameter(15, brPsiChi2);
  fullModel.SetParameter(16, brChi1Jpsi);
  fullModel.SetParameter(17, brChi2Jpsi);
  fullModel.SetParameter(18, brPsiJpsi);
  models.push_back(fullModel);

  TF1 feedDown("jpsi_cs_feed_down",
               [] (double* x, double* p) {
                 const double psi = sig_dir(x[0], p[0], p[1], p[2], p[3], p[4]);
                 const double dirChi1 = sig_dir(x[0], p[5], p[6], p[7], p[8], p[4]);
                 const double dirChi2 = sig_dir(x[0], p[9], p[10], p[11], p[12], p[4]);

                 const double chi1 = dirChi1 + psi * p[13];
                 const double chi2 = dirChi2 + psi * p[14];

                 return chi1 * p[15] + chi2 * p[16] + psi * p[17];
               }, MIN_PTM, 40, 18);
  feedDown.SetParameters(sigmaPsi, fLongPsi, betaLongPsi, betaTransPsi, gamma,
                         sigma1, fLong1, betaLong1, betaTrans1,
                         sigma2, fLong2);
  feedDown.SetParameter(11, betaLong2);
  feedDown.SetParameter(12, betaTrans2);
  feedDown.SetParameter(13, brPsiChi1);
  feedDown.SetParameter(14, brPsiChi2);
  feedDown.SetParameter(15, brChi1Jpsi);
  feedDown.SetParameter(16, brChi2Jpsi);
  feedDown.SetParameter(17, brPsiJpsi);

  models.push_back(feedDown);

  return models;
}

/**
 * point-by-point corrections for the data points for the psi2S cross sections
 */
template<typename PsiPolModel>
std::vector<double> psi2SCorrections(const CrossSectionMeasurement& data, PsiPolModel polModel, double globalCorr=1.0)
{
  std::vector<double> corrs;
  for (const auto& point : data) {
    const double lambda = polModel(point.ptM);
    corrs.push_back(acceptanceCorrection(point.K, lambda) * globalCorr);
  }

  return corrs;
}

/**
 * point-by-point cross section corrections for chic measurements
 */
template<typename ChiPolModel>
std::vector<double> chicCorrections(const CrossSectionMeasurement& data, CSModel psiCSModel,
                                    CSModel chiCSModel, ChiPolModel polModel,
                                    double brPsiChi, double mChi, double nuBrPsiChi,
                                    double globalCorr=1.0)
{
  std::vector<double> corrs;
  for (const auto& point : data) {
    const double dirCSPsi = sampleFunc(point.low, point.high + 1e-5, 0.5 / mChi, psiCSModel);
    const double dirCSChi = sampleFunc(point.low, point.high + 1e-5, 0.5 / mChi, chiCSModel);

    const double csChi = dirCSChi + dirCSPsi * brPsiChi / nuBrPsiChi;
    const double lambda = polModel(point.ptM, dirCSChi / csChi);
    corrs.push_back(acceptanceCorrection(point.K, lambda) * globalCorr);
  }

  return corrs;
}

template<typename JpsiPolModel>
std::vector<double> jpsiCorrections(const CrossSectionMeasurement& data, CSModel psiCSModel,
                                    CSModel chi1CSModel, CSModel chi2CSModel, CSModel jpsiCSModel,
                                    JpsiPolModel polModel, double brPsiChi1, double brPsiChi2,
                                    double brChi1Jpsi, double brChi2Jpsi, double brPsiJpsi,
                                    double globalCorr=1.0)
{
  std::vector<double> corrs;
  for (const auto& point : data) {
    const double dirCSPsi = sampleFunc(point.low, point.high + 1e-5, 0.5 / M_JPSI, psiCSModel);
    const double dirCSChi1 = sampleFunc(point.low, point.high + 1e-5, 0.5 / M_JPSI, chi1CSModel);
    const double dirCSChi2 = sampleFunc(point.low, point.high + 1e-5, 0.5 / M_JPSI, chi2CSModel);
    const double dirCSJpsi = sampleFunc(point.low, point.high + 1e-5, 0.5 / M_JPSI, jpsiCSModel);

    const double csChi1 = dirCSChi1 + dirCSPsi * brPsiChi1;
    const double fdChi1 = csChi1 * brChi1Jpsi;
    const double csChi2 = dirCSChi2 * dirCSPsi * brPsiChi2;
    const double fdChi2 = csChi2 * brChi2Jpsi;

    const double csJpsi = dirCSJpsi + fdChi1 + fdChi2 + dirCSJpsi * brPsiJpsi;
    const double lambda = polModel(point.ptM, dirCSChi1 / csChi1, dirCSChi2 / csChi2,
                                   csJpsi / dirCSJpsi, fdChi1 / csJpsi, fdChi2 / csJpsi);

    corrs.push_back(acceptanceCorrection(point.K, lambda) * globalCorr);
  }

  return corrs;
}

TGraphAsymmErrors correctCSData(const CrossSectionMeasurement& rawData, const char* name,
                                const std::vector<double>& corrFact)
{
  auto graph = asTGraph(rawData);
  for (int i = 0; i < graph.GetN(); ++i) {
    double x, y;
    graph.GetPoint(i, x, y);
    graph.SetPoint(i, x, y * corrFact[i]);
    graph.SetPointEYhigh(i, graph.GetErrorYhigh(i) * corrFact[i]);
    graph.SetPointEYlow(i, graph.GetErrorYlow(i) * corrFact[i]);
  }

  graph.SetName(name);
  return graph;
}


#endif
