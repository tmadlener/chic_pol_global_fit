#ifndef H_MULTIVARIATE_NORMAL_DISTRIBUTON__
#define H_MULTIVARIATE_NORMAL_DISTRIBUTON__

#include "Eigen/Dense"
#include <random>
#include <vector>
#include <iostream>


/**
 * Multivariate normal distribution generator
 *
 * see: https://stackoverflow.com/a/40245513
 */
class MultivariateNormalDistribution {
public:
  MultivariateNormalDistribution(std::vector<double>& means,
                                 std::vector<double>& cov)
  {
    m_mean = Eigen::Map<Eigen::Matrix<double, 1, Eigen::Dynamic>>(means.data(), means.size());

    Eigen::MatrixXd covMatrix = Eigen::Map<Eigen::MatrixXd>(cov.data(), means.size(), means.size());

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(covMatrix);
    if ((solver.eigenvalues().array() < 0).any()) {
      std::cerr << "WARNING: covariance matrix is not positive definite!\n";
    }
    m_transform = solver.eigenvectors() * solver.eigenvalues().cwiseSqrt().asDiagonal();
  }

  std::vector<double> operator()() const
  {
    static std::mt19937 gen{std::random_device{}()};
    static std::normal_distribution<> dist;

    const Eigen::VectorXd vals = m_mean +
      m_transform * Eigen::VectorXd{m_mean.size()}.unaryExpr([&](const double) { return dist(gen); });

    return std::vector<double>(&vals[0], vals.data() + vals.size());
  }

private:
  Eigen::VectorXd m_mean;
  Eigen::MatrixXd m_transform;
};


#endif
