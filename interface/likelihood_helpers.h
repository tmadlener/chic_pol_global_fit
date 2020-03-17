#ifndef H_LIKELIHOOD_HELPERS__
#define H_LIKELIHOOD_HELPERS__

#include "data_structures.h"
#include "constants.h"

#include <cmath>
#include <functional>

using CSModel = std::function<double(double)>;
using PsiPolModel = std::function<double(double)>;
using ChiPolModel = std::function<double(double, double)>;
using JpsiPolModel = std::function<double(double, double, double, double, double, double)>;


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

double psi2SCrossSection(const CrossSectionMeasurement& data, CSModel csModel, PsiPolModel polModel,
                         double nuissLumi, double nuissBrFrac)
{
  double loglike = 0;

  for (const auto& point : data) {
    const double dirCS = sampleFunc(point.low, point.high + 1e-5, 0.5 / M_PSI2S, csModel);
    const double lambda = polModel(point.ptM);
    const double accCorr = acceptanceCorrection(point.K, lambda);

    const double crossSec = point.xSec * accCorr * nuissLumi * nuissBrFrac;
    const double uncer = point.uncer * accCorr * nuissLumi * nuissBrFrac;

    const double relDiff = (dirCS - crossSec) / uncer;
    loglike += -0.5 * relDiff * relDiff;
  }

  return loglike;
}

double psi2SPolarization(const PolarizationMeasurement& data, PsiPolModel polModel)
{
  double loglike = 0;

  for (const auto& point : data) {
    const double lambda = polModel(point.ptM);

    const double relDiff = (lambda - point.lth) / point.uncer;
    loglike += -0.5 * relDiff * relDiff;
  }

  return loglike;
}

double chicCrossSection(const CrossSectionMeasurement& data, CSModel psiCSModel, CSModel chiCSModel,
                        ChiPolModel chiPolModel,
                        const double nuLumi, const double nuBrPsiChi, const double nuBrChiDet,
                        const double brPsiChi, const double mChi)
{
  double loglike = 0;
  for (const auto& point : data) {
    const double dirCSPsi = sampleFunc(point.low, point.high + 1e-5, 0.5 / mChi, psiCSModel);
    const double dirCSChi = sampleFunc(point.low, point.high + 1e-5, 0.5 / mChi, chiCSModel);

    const double csChi = dirCSChi + dirCSPsi * brPsiChi / nuBrPsiChi;
    const double lambda = chiPolModel(point.ptM, dirCSChi / csChi);
    const double accCorr = acceptanceCorrection(point.K, lambda);

    const double crossSec = point.xSec * accCorr * nuLumi * nuBrChiDet;
    const double uncer = point.uncer * accCorr * nuLumi * nuBrChiDet;

    const double relDiff = (csChi - crossSec) / uncer;
    loglike += -0.5 * relDiff * relDiff;
  }

  return loglike;
}

double chicCrossSectionRatio(const CrossSectionMeasurement& data, CSModel psiCSModel,
                             CSModel chi1CSModel, CSModel chi2CSModel,
                             ChiPolModel chi1PolModel, ChiPolModel chi2PolModel,
                             const double nuBrPsiChi1, const double nuBrPsiChi2)
{
  double loglike = 0;
  for (const auto& point : data) {
    const double dirCSPsi = sampleFunc(point.low, point.high + 1e-5, 0.5 / M_JPSI, psiCSModel);
    const double dirCSChi1 = sampleFunc(point.low, point.high + 1e-5, 0.5 / M_JPSI, chi1CSModel);
    const double dirCSChi2 = sampleFunc(point.low, point.high + 1e-5, 0.5 / M_JPSI, chi2CSModel);

    const double csChi1 = dirCSChi1 + dirCSPsi * B_PSIP_CHIC1[0] / nuBrPsiChi1;
    const double csChi2 = dirCSChi2 + dirCSPsi * B_PSIP_CHIC2[0] / nuBrPsiChi2;

    const double lambda1 = chi1PolModel(point.ptM, dirCSChi1 / csChi1);
    const double lambda2 = chi2PolModel(point.ptM, dirCSChi2 / csChi2);

    const double accCorr = acceptanceCorrection(point.K, lambda2) / acceptanceCorrection(point.K, lambda1);

    const double crossSecRatio = point.xSec * accCorr * nuBrPsiChi2 / nuBrPsiChi1;
    const double uncer = point.uncer * accCorr * nuBrPsiChi2 / nuBrPsiChi1;

    const double relDiff = (csChi2 / csChi1 - crossSecRatio) / uncer;
    loglike += -0.5 * relDiff * relDiff;
  }

  return loglike;
}

double jpsiCrossSection(const CrossSectionMeasurement& data, CSModel psiCSModel,
                        CSModel chi1CSModel, CSModel chi2CSModel, CSModel jpsiCSModel,
                        JpsiPolModel polModel,
                        const double nuBrPsiChi1, const double nuBrPsiChi2, const double nuBrPsiJpsi,
                        const double nuBrChi1Jpsi, const double nuBrChi2Jpsi,
                        const double nuBrJpsiMm, const double nuLumi)
{
  double loglike = 0;

  for (const auto& point : data) {
    const double dirCSPsi = sampleFunc(point.low, point.high + 1e-5, 0.5 / M_JPSI, psiCSModel);
    const double dirCSChi1 = sampleFunc(point.low, point.high + 1e-5, 0.5 / M_JPSI, chi1CSModel);
    const double dirCSChi2 = sampleFunc(point.low, point.high + 1e-5, 0.5 / M_JPSI, chi2CSModel);
    const double dirCSJpsi = sampleFunc(point.low, point.high + 1e-5, 0.5 / M_JPSI, jpsiCSModel);

    const double csChi1 = dirCSChi1 + dirCSPsi * B_PSIP_CHIC1[0] / nuBrPsiChi1;
    const double csChi2 = dirCSChi2 + dirCSPsi * B_PSIP_CHIC2[0] / nuBrPsiChi2;

    const double fdChi1 = csChi1 * B_CHIC1_JPSI[0] / nuBrChi1Jpsi;
    const double fdChi2 = csChi2 * B_CHIC2_JPSI[0] / nuBrChi2Jpsi;

    const double csJpsi = dirCSJpsi + dirCSPsi * B_PSIP_JPSI[0] / nuBrPsiJpsi + fdChi1 + fdChi2;

    const double lambda = polModel(point.ptM, dirCSChi1 / csChi1, dirCSChi2 / csChi2,
                                   csJpsi / dirCSJpsi, fdChi1 / csJpsi, fdChi2 / csJpsi);
    const double accCorr = acceptanceCorrection(point.K, lambda);

    const double crossSec = point.xSec * accCorr * nuLumi * nuBrJpsiMm;
    const double uncer = point.uncer * accCorr * nuLumi * nuBrJpsiMm;

    const double relDiff = (csJpsi - crossSec) / uncer;
    loglike += -0.5 * relDiff * relDiff;
  }

  return loglike;
}

double jpsiPolarization(const PolarizationMeasurement& data, JpsiPolModel polModel,
                        CSModel psiCSModel, CSModel chi1CSModel, CSModel chi2CSModel, CSModel jpsiCSModel,
                        const double nuBrPsiChi1, const double nuBrPsiChi2, const double nuBrPsiJpsi,
                        const double nuBrChi1Jpsi, const double nuBrChi2Jpsi)
{
  double loglike = 0;
  for (const auto& point : data) {
    const double dirCSPsi = sampleFunc(point.low, point.high + 1e-5, 0.5 / M_JPSI, psiCSModel);
    const double dirCSChi1 = sampleFunc(point.low, point.high + 1e-5, 0.5 / M_JPSI, chi1CSModel);
    const double dirCSChi2 = sampleFunc(point.low, point.high + 1e-5, 0.5 / M_JPSI, chi2CSModel);
    const double dirCSJpsi = sampleFunc(point.low, point.high + 1e-5, 0.5 / M_JPSI, jpsiCSModel);

    const double csChi1 = dirCSChi1 + dirCSPsi * B_PSIP_CHIC1[0] / nuBrPsiChi1;
    const double csChi2 = dirCSChi2 + dirCSPsi * B_PSIP_CHIC2[0] / nuBrPsiChi2;

    const double fdChi1 = csChi1 * B_CHIC1_JPSI[0] / nuBrChi1Jpsi;
    const double fdChi2 = csChi2 * B_CHIC2_JPSI[0] / nuBrChi2Jpsi;

    const double csJpsi = dirCSJpsi + dirCSPsi * B_PSIP_JPSI[0] / nuBrPsiJpsi + fdChi1 + fdChi2;

    const double lambda = polModel(point.ptM, dirCSChi1 / csChi1, dirCSChi2 / csChi2,
                                   csJpsi / dirCSJpsi, fdChi1 / csJpsi, fdChi2 / csJpsi);

    const double relDiff = (lambda - point.lth) / point.uncer;
    loglike += -0.5 * relDiff * relDiff;
  }

  return loglike;
}

#endif
