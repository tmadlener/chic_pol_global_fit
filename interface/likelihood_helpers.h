#ifndef H_LIKELIHOOD_HELPERS__
#define H_LIKELIHOOD_HELPERS__

#include "data_structures.h"
#include "constants.h"

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







#endif
