#ifndef MATRIX_HELPERS_H__
#define MATRIX_HELPERS_H__

#include <vector>
#include <utility>
#include <cmath>
#include <iostream>

/**
 * Get the number of parameters from the covariance matrix
 */
size_t n_parameters(const std::vector<double>& cov)
{
  const size_t nPars = std::sqrt(cov.size());
  if (nPars * nPars != cov.size()) {
    std::cerr << "covariance matrix is not square\n";
    return {};
  }
  return nPars;
}

/**
 * Get the variances (i.e. the diagonal) from the covariance matrix
 */
std::vector<double> variances(const std::vector<double>& cov)
{
  const size_t nPars = n_parameters(cov);
  std::vector<double> vars(nPars);
  for (size_t i = 0; i < nPars; ++i) {
    vars[i] = cov[i * nPars + i];
  }
  return vars;
}

/**
 * Get the correlation matrix from the covariance matrix, with pre-computed
 * variances
 */
std::vector<double> corr_matrix(const std::vector<double>& cov,
                                const std::vector<double>& vars)
{
  const size_t nPars = n_parameters(cov);

  std::vector<double> corr(nPars * nPars);

  for (size_t i = 0; i < nPars; ++i) {
    for (size_t j = 0; j < nPars; ++j) {
      const double pwCov = vars[i] * vars[j];
      corr[i * nPars + j] = cov[i * nPars + j] / std::sqrt(pwCov);
    }
  }

  return corr;
}

/**
 * Get the correlation matrix from the covariance matrix
 */
std::vector<double> corr_matrix(const std::vector<double>& cov)
{
  return corr_matrix(cov, variances(cov));
}

/**
 * Get a covariance matrix with the passed variances and correlations
 */
std::vector<double> cov_matrix(const std::vector<double>& corr,
                               const std::vector<double>& vars)
{
  const size_t nPars = vars.size();
  if (nPars * nPars != corr.size()) {
    std::cerr << "correlation matrix and variances do not match in their number of variables\n";
    return {};
  }

  std::vector<double> cov(nPars * nPars);
  for (size_t i = 0; i < nPars; ++i) {
    for (size_t j = 0; j < nPars; ++j) {
      const double pwCov = vars[i] * vars[j];
      cov[i * nPars + j] = corr[i * nPars + j] * std::sqrt(pwCov);
    }
  }

  return cov;
}

/**
 * Change the variances of the passed covariance matrix, while keeping the
 * correlations intact.
 *
 * The factors are the indices and multiplicative factors by which the variances
 * should be changed for the individual parameters
 */
std::vector<double> change_variance(const std::vector<double>& cov,
                                    const std::vector<std::pair<int, double>>& factors)
{
  std::vector<double> vars = variances(cov);
  std::vector<double> corr = corr_matrix(cov, vars);

  for (const auto& fact : factors) {
    vars[fact.first] *= fact.second;
  }

  return cov_matrix(corr, vars);
}

/**
 * Reduce the variance of all variables by the passed factor
 */
std::vector<double> reduce_variance(const std::vector<double>& cov, const double factor)
{
  std::vector<double> vars = variances(cov);
  for (auto& var : vars) {
    var *= factor;
  }
  return cov_matrix(corr_matrix(cov), vars);
}


#endif
