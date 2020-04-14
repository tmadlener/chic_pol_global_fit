#ifndef H_MULTIVARIATE_NORMAL_DISTRIBUTON__
#define H_MULTIVARIATE_NORMAL_DISTRIBUTON__

#include "Eigen/Dense"

#include <random>
#include <vector>
#include <iostream>
#include <algorithm>


/**
 * Multivariate normal distribution generator
 *
 * see: https://stackoverflow.com/a/40245513
 */
template<typename RNG=std::mt19937>
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


    // properly seed the mt19937 (see: https://codereview.stackexchange.com/a/109266)
    std::array<unsigned, RNG::state_size * RNG::word_size> random_data;
    std::random_device source{};
    std::generate(random_data.begin(), random_data.end(), std::ref(source));
    std::seed_seq seeds(random_data.begin(), random_data.end());
    m_gen.seed(seeds);
  }

  std::vector<double> operator()() const
  {
    const Eigen::VectorXd vals = m_mean +
      m_transform * Eigen::VectorXd{m_mean.size()}.unaryExpr([&](const double) { return m_dist(m_gen); });

    return std::vector<double>(&vals[0], vals.data() + vals.size());
  }

private:
  Eigen::VectorXd m_mean;
  Eigen::MatrixXd m_transform;

  mutable RNG m_gen{};
  mutable std::normal_distribution<> m_dist;
};


#endif
