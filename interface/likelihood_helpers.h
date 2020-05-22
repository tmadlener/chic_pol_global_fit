#ifndef H_LIKELIHOOD_HELPERS__
#define H_LIKELIHOOD_HELPERS__

#include "data_structures.h"
#include "constants.h"

#include "TF1.h"
#include "TGraphAsymmErrors.h"

#include <cmath>
#include <functional>
#include <utility>
#include <numeric>

using CSModel = std::function<double(double)>;
using PolModel = std::function<double(double)>;
using PolFeedDownTrafo = std::function<double(double)>;
using CosthRatioModel = std::function<double(double)>;

double power_law(double x, double beta, double gamma)
{
  return x * std::pow(1 + x * x / (beta - 2) / gamma, -beta);
}

double norm_plaw(double x, double beta, double gamma, double xnorm=PTMNORM)
{
  return power_law(x, beta, gamma) / power_law(xnorm, beta, gamma);
}

double sig_dir(double ptm, double sig, double beta, double gamma)
{
  return sig * norm_plaw(ptm, beta, gamma);
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
 * lambda theta assuming a fully transverse and a fully longitudinal component
 * given the fraction of the longitudinal component
 */
double lambdath(const double fLong)
{
  return (1 - 3 * fLong) / (1 + fLong);
}

double fracLong(const double lth)
{
  return (1 - lth) / (3 + lth);
}

/**
 * ptm dependent lambda theta from the shape parameters of the two total and
 * longitudional cross sections and the longitudinal fraction at a given
 * reference point
 */
double lambdaTheta(double ptm, double fLong, double betaLong, double betaTotal, double gamma)
{
  const double total = norm_plaw(ptm, betaTotal, gamma);
  const double longCont = fLong * norm_plaw(ptm, betaLong, gamma);

  return lambdath(longCont / total);
}

/**
 * ptm dependent lambda theta from the shape parameters and lambda theta at a
 * reference point
 */
double lambdaThetaAtPtm(double ptm, double lth, double betaLong, double betaTotal, double gamma)
{
  return lambdaTheta(ptm, fracLong(lth), betaLong, betaTotal, gamma);
}

/**
 * clip the passed value into the passed range
 */
double clipRange(const std::pair<double, double>& range, const double val)
{
  return val < range.first ? range.first : (val > range.second ? range.second : val);
}

/**
 * Calculate the polarization dependent correction for the measured cross section
 */
double acceptanceCorrection(const double kFactor, const double lambda,
                            const bool clip=false, const std::pair<double, double>& range={0, 0})
{
  const double lam = clip ? clipRange(range, lambda) : lambda;

  return (1 + 1. / 3. * kFactor) / (1 + (1 - lam) / (3 + lam) * kFactor);
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
 * Overload for evaluating the cross sections at only one ptM point instead of
 * integrating them over a whole bin
 */
std::pair<double, double> crossSecAndLambda(const double ptM,
                                            const std::vector<CSModel>& csModels,
                                            const std::vector<PolModel>& polModels,
                                            const std::vector<PolFeedDownTrafo>& feedDownTrafos,
                                            const std::vector<double>& brFracs)
{
  // constructing a "helper" BinInfo here in such a way that it will only be
  // evaluated once at its minimum value during the call to sampleFunc in
  // partialCrossSections
  return crossSecAndLambda(BinInfo{ptM, ptM, ptM}, 1e9, csModels, polModels, feedDownTrafos, brFracs);
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
                           const double mass,
                           const bool clipCorrections=false,
                           const std::pair<double, double>& lambdaRange={-1.0, 1.0})
{
  double loglike = 0;
  for (const auto& point : data) {
    double cs, lambda;
    std::tie(cs, lambda) = crossSecAndLambda(point, 0.5 / mass,
                                             csModels, polModels, feedDownTrafos, brFracs);

    const double accCorr = acceptanceCorrection(point.K, lambda, clipCorrections, lambdaRange);
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
                                const double globNuiss, const double mass,
                                const bool clipCorrections=false,
                                const std::pair<double, double>& lambdaRangeN={-0.6, 1.0},
                                const std::pair<double, double>& lambdaRangeD={-1./3., 1.0})
{
  double loglike = 0;
  for (const auto& point : data) {
    double csN, lambdaN, csD, lambdaD;
    std::tie(csN, lambdaN) = crossSecAndLambda(point, 0.5 / mass,
                                               csModelsN, polModelsN, fdTrafosN, brFracsN);
    std::tie(csD, lambdaD) = crossSecAndLambda(point, 0.5 / mass,
                                               csModelsD, polModelsD, fdTrafosD, brFracsD);

    const double accCorr = acceptanceCorrection(point.K, lambdaN, clipCorrections, lambdaRangeN) / acceptanceCorrection(point.K, lambdaD, clipCorrections, lambdaRangeD);
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



/**
 * Analytic costh ratio function
 */
double costhRatio(const double costh, const double lamN, const double lamD, const double norm)
{
  return norm * (1 + lamN * costh * costh) / (1 + lamD * costh * costh);
}

double loglikeCosthRatio(const CosthRatioMeasurement& data,
                         const CosthRatioModel& model)
{
  double loglike = 0;

  for (const auto& point : data) {
    const double pred = model(point.costh);
    // use the correct uncertainties depending on where the prediction lies
    // relative to the data. prediction above data, means upper uncertainties.
    const double uncer = (pred > point.ratio) ? point.u_high : point.u_low;
    const double relDiff = (pred - point.ratio) / uncer;
    loglike += -0.5 * relDiff * relDiff;
  }

  return loglike;
}

const std::vector<double> getCorrectionFactors(const CrossSectionMeasurement& data,
                                               const std::vector<CSModel>& csModels,
                                               const std::vector<PolModel>& polModels,
                                               const std::vector<PolFeedDownTrafo>& fdTrafos,
                                               const std::vector<double> brFracs,
                                               const double globNuiss,
                                               const double mass,
                                               const bool clipCorrs=false,
                                               const std::pair<double, double>& lambdaRange={-1.0, 1.0})
{
  std::vector<double> corrections;
  for (const auto& point : data) {
    double lambda;
    std::tie(std::ignore, lambda) = crossSecAndLambda(point, 0.5 / mass,
                                                      csModels, polModels, fdTrafos, brFracs);

    corrections.push_back(acceptanceCorrection(point.K, lambda, clipCorrs, lambdaRange) * globNuiss);
  }

  return corrections;
}


TGraphAsymmErrors correctGraph(const CrossSectionMeasurement& data, const char* name,
                               const std::vector<double>& corrFactors)
{
  auto graph = asTGraph(data);

  for (int i = 0; i < graph.GetN(); ++i) {
    double x, y;
    graph.GetPoint(i, x, y);
    graph.SetPoint(i, x, y * corrFactors[i]);
    graph.SetPointEYhigh(i, graph.GetErrorYhigh(i) * corrFactors[i]);
    graph.SetPointEYlow(i, graph.GetErrorYlow(i) * corrFactors[i]);
  }

  graph.SetName(name);
  return graph;
}


TGraphAsymmErrors correctedCSGraph(const CrossSectionMeasurement& data,
                                   const std::vector<CSModel>& csModels,
                                   const std::vector<PolModel>& polModels,
                                   const std::vector<PolFeedDownTrafo>& fdTrafos,
                                   const std::vector<double> brFracs,
                                   const double globNuiss,
                                   const double mass,
                                   const char* name,
                                   const bool clipCorrs=false,
                                   const std::pair<double, double>& lambdaRange={-1.0, 1.0}){
  const auto corrections = getCorrectionFactors(data, csModels, polModels, fdTrafos,
                                                brFracs, globNuiss, mass, clipCorrs, lambdaRange);

  return correctGraph(data, name, corrections);
}

/**
 * wrap the model in a lambda and evaluate it at x[0]
 *
 * NOTE: taking the model by const value to make sure everything that it is
 * closed over is still present when it is actually evaluated
 */
template<typename Model>
TF1 modelAsTF1(const Model model, const char* name, const double min, const double max)
{
  return TF1(name,
             [model](double* x, double*) { return model(x[0]); },
             min, max, 0);
}

/**
 * Get the total model including feed-down contributions
 *
 * NOTE: Taking the CSModels and PolModels by const value to make sure
 * everything that they close over is still present when it is actually
 * evaluated
 */
std::pair<CSModel, PolModel> combineModels(const std::vector<CSModel> csModels,
                                           const std::vector<PolModel> polModels,
                                           const std::vector<PolFeedDownTrafo>& fdTrafos,
                                           const std::vector<double> brFracs)
{
  const CSModel csModel = [=](double ptm) -> double {
    const auto csLambda = crossSecAndLambda(ptm, csModels, polModels, fdTrafos, brFracs);
    return csLambda.first;
  };

  const PolModel polModel = [=](double ptm) -> double {
    const auto csLambda = crossSecAndLambda(ptm, csModels, polModels, fdTrafos, brFracs);
    return csLambda.second;
  };

  return {csModel, polModel};
}

#endif
