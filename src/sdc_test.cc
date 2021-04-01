#include "sdc.h"

#include <iostream>
#include <string>
#include <utility>

constexpr static double ccbarMassSDCs = 3.0;

/**
 * Read the total and longitudinal SDCs of the psi states. I.e build a MultiSDC
 * that represents all the intermediate states that contribute to direct psi
 * production and that can be evaluated using corresponding LDMEs.
 */
std::pair<sdc::MultiSDC, sdc::MultiSDC> readPsiSDCs(std::string dataDir) {
  // Read all the relevant sdcs from file
  const auto sdc_3S1_1_unpol = sdc::read_from_file(dataDir + "/3s1singletunpol.txt");
  const auto sdc_3S1_1_trans = sdc::read_from_file(dataDir + "/3s1singlettransverse.txt");
  const auto sdc_3S1_8_unpol = sdc::read_from_file(dataDir + "/3s1octetunpol.txt");
  const auto sdc_3S1_8_long = sdc::read_from_file(dataDir + "/3s1octetlong.txt");
  const auto sdc_3PJ_8_unpol = sdc::read_from_file(dataDir + "/3pjoctetunpol.txt");
  const auto sdc_3PJ_8_long = sdc::read_from_file(dataDir + "/3pjoctetlong.txt");
  const auto sdc_1S0_8 = sdc::read_from_file(dataDir + "/1s0octetunpol.txt");

  // the singlet is unpolarized, hence the three components 0, +1, -1 are equal
  // -> long = 1/3 * total
  const auto sdc_1S0_8_long = 1. / 3 * sdc_1S0_8;

  // TODO: the actually correct combinations of these to have total and
  // longitudinal SDCs
  auto longitudinal = sdc::MultiSDC({sdc_3S1_1_trans, sdc_3S1_8_long, sdc_3PJ_8_long, sdc_1S0_8_long});
  longitudinal.scale_support(1. / ccbarMassSDCs); // we want them to be in pT/M
  auto total = sdc::MultiSDC({sdc_3S1_1_unpol, sdc_3S1_8_unpol, sdc_3PJ_8_unpol, sdc_1S0_8});
  total.scale_support(1. / ccbarMassSDCs); // we want them to be in pT/M

  return {total, longitudinal};
}

/**
 * Read the total and longitudinal SDCs for the chic1. I.e. build a MultiSDC
 * that represents all the intermediate states that contribute to direct chi
 * production as observable after their decay to a J/psi since that is what the
 * fitting code uses (and also what the chi polarization measurements provide).
 */
std::pair<sdc::MultiSDC, sdc::MultiSDC> readChic1SDCs(std::string dataDir) {
  const auto sdc_3P1_1_unpol = sdc::read_from_file(dataDir + "/3p1singletunpol.txt");
  const auto sdc_3P1_1_h1 = sdc::read_from_file(dataDir + "/3p1singletpol1.txt");
  const auto sdc_3S1_8_unpol = sdc::read_from_file(dataDir + "/3s1octetunpol.txt");
  const auto sdc_3S1_8_long = sdc::read_from_file(dataDir + "/3s1octetlong.txt");

  // TODO: check with Pietro what the correct combinations are to arrive at the
  // total and longitudinal SDCs
  auto longitudinal = sdc::MultiSDC({sdc_3P1_1_h1, sdc_3S1_8_long});
  longitudinal.scale_support(1. / ccbarMassSDCs);

  auto total = sdc::MultiSDC({sdc_3P1_1_unpol, sdc_3S1_8_unpol});
  longitudinal.scale_support(1. / ccbarMassSDCs);

  return {total, longitudinal};
}

/**
 * Read the total and longitudinal SDCs for the chic2. I.e. build a MultiSDC
 * that represents all the intermediate states that contribute to direct chi
 * production as observable after their decay to a J/psi since that is what the
 * fitting code uses (and also what the chi polarization measurements provide).
 */
std::pair<sdc::MultiSDC, sdc::MultiSDC> readChic2SDCs(std::string dataDir) {
  const auto sdc_3P2_1_unpol = sdc::read_from_file(dataDir + "/3p2singletunpol.txt");
  const auto sdc_3P2_1_h1 = sdc::read_from_file(dataDir + "/3p2singletpol1.txt");
  const auto sdc_3P2_1_h2 = sdc::read_from_file(dataDir + "/3p2singletpol2.txt");

  const auto sdc_3S1_8_unpol = sdc::read_from_file(dataDir + "/3s1octetunpol.txt");
  const auto sdc_3S1_8_long = sdc::read_from_file(dataDir + "/3s1octetlong.txt");

  // TODO: check with Pietro what the correct combinations are to arrive at the
  // total and longitudinal SDCs
  auto longitudinal = sdc::MultiSDC({sdc_3P2_1_unpol, sdc_3P2_1_h1, sdc_3P2_1_h2});
  longitudinal.scale_support(1. / ccbarMassSDCs);

  auto total = sdc::MultiSDC({sdc_3S1_8_unpol, sdc_3S1_8_long});
  total.scale_support(1. / ccbarMassSDCs);

  return {total, longitudinal};
}

int main() {

  const auto [psiSDCTotal, psiSDCLong] = readPsiSDCs("../SDCs_Chung");

  std::cout << psiSDCTotal << std::endl;
  std::cout << "==========================" << std::endl;
  std::cout << psiSDCLong << std::endl;

  // test the multi sdc evaluation
  //
  // NOTE: LDME 3 has to be small because SDC is negative for that. This should
  // automatically be good for the fit
  std::cout << psiSDCTotal(10.5, {1.3, 3.4, 0.1, 7.8}, sdc::util::log_interpolate) << " "
            << psiSDCTotal(10.5, {1.3, 3.4, 0.1, 7.8}) << std::endl;

  std::cout << "--------------------------" << std::endl;
  const auto [chic1SDCTotal, chic1SDCLong] = readChic1SDCs("../SDCs_Chung");

  std::cout << chic1SDCTotal << std::endl;
  std::cout << "==========================" << std::endl;
  std::cout << chic1SDCLong << std::endl;

  std::cout << "--------------------------" << std::endl;
  const auto [chic2SDCTotal, chic2SDCLong] = readChic2SDCs("../SDCs_Chung");

  std::cout << chic2SDCTotal << std::endl;
  std::cout << "==========================" << std::endl;
  std::cout << chic2SDCLong << std::endl;

  return 0;
}
