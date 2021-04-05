#ifndef GLOBALFIT_READSDCS_H
#define GLOBALFIT_READSDCS_H

#include "sdc.h"

#include <string>
#include <utility>

constexpr static double ccbarMassSDCs = 3.0;

/**
 * Read the total and longitudinal SDCs of the psi states. I.e build a MultiSDC
 * that represents all the intermediate states that contribute to direct psi
 * production and that can be evaluated using corresponding LDMEs.
 *
 * The order of contributions is (see also nrqcd_helpers.h):
 * - 3S1_1
 * - 3S1_8
 * - 3PJ_8
 * - 1S0_8
 */
sdc::StateSDCs readPsiSDCs(std::string dataDir, sdc::SDCType sdcType=sdc::SDCType::LP_NLO) {
  // For further processing using the naming convention:
  // * S_U -> unpolarized SDC (== S_total)
  // * S_T -> transverse SDC (NOTE: only referring to either J_z = +1 or -1, not the combination)
  // * S_L -> longitudinal SDC

  // 3S1 singlet
  const auto sdc_3S1_1_unpol = sdc::read_from_file(dataDir + "/3s1singletunpol.txt", sdcType);
  const auto sdc_3S1_1_trans = sdc::read_from_file(dataDir + "/3s1singlettransverse.txt", sdcType);
  // It holds: S_U = S_L + 2 * S_T
  // --> S_L = S_U - 2 * S_T
  const auto sdc_3S1_1_long = sdc_3S1_1_unpol - 2 * sdc_3S1_1_trans;

  // 3S1 octet
  const auto sdc_3S1_8_unpol = sdc::read_from_file(dataDir + "/3s1octetunpol.txt", sdcType);
  const auto sdc_3S1_8_long = sdc::read_from_file(dataDir + "/3s1octetlong.txt", sdcType);
  // Following the same logic as above:
  // --> S_T = 0.5 * (S_U - S_L)

  // 3PJ octet
  const auto sdc_3PJ_8_unpol = sdc::read_from_file(dataDir + "/3pjoctetunpol.txt", sdcType);
  const auto sdc_3PJ_8_long = sdc::read_from_file(dataDir + "/3pjoctetlong.txt", sdcType);

  // 1S0 octet
  const auto sdc_1S0_8 = sdc::read_from_file(dataDir + "/1s0octetunpol.txt", sdcType);
  // the octet is unpolarized, hence the three components 0, +1, -1 are equal
  // -> long = 1/3 * total
  const auto sdc_1S0_8_long = 1. / 3 * sdc_1S0_8;

  // Put everything into the combination MultiSDCs
  auto longitudinal = sdc::MultiSDC({sdc_3S1_1_long, sdc_3S1_8_long, sdc_3PJ_8_long, sdc_1S0_8_long});
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
 *
 * The order of the contributions is (see als nrqcd_helpers.h):
 * - 3P1_1
 * - 3S1_8
 */
sdc::StateSDCs readChic1SDCs(std::string dataDir, sdc::SDCType sdcType=sdc::SDCType::LP_NLO) {
  // For further processing using the naming convention:
  // * S_U -> unpolarized SDC (== total)
  // * S_0 -> SDC of J_z = 0
  // * S_1 -> SDC of J_z = +/- 1
  //
  // To calculate the longitudinal SDC from the chi states we want them to yield
  // the observable J/psi polariztion (in the radiative decay) not the direct
  // chic polarization. To achieve this we first calculate the the lambda theta
  // for each SDC to get to the longitudinal fraction and then use that to
  // determine the longitudinal SDC from the total one. The formulas used can be
  // found in PRD 83, 096001 (2011)

  // 3P1 singlet
  const auto sdc_3P1_1_unpol = sdc::read_from_file(dataDir + "/3p1singletunpol.txt", sdcType);
  const auto sdc_3P1_1_h1 = sdc::read_from_file(dataDir + "/3p1singletpol1.txt", sdcType);
  // We have S_total = S_U = S_0 + 2 * S_1
  // --> S_0 = S_U - 2 * S_1
  const auto sdc_3P1_1_h0 = sdc_3P1_1_unpol - 2 * sdc_3P1_1_h1;
  // lth = (3P1_1_h0 - 3P1_1_h1) / (3P1_1_h0 + 3 * 3P1_1_h1)
  const auto lth_3P1 = (sdc_3P1_1_h0 - sdc_3P1_1_h1) / (sdc_3P1_1_h0 + 3 * sdc_3P1_1_h1);
  const auto fL_3P1 = (1 - lth_3P1) / (3 + lth_3P1);
  const auto sdc_3P1_1_long = fL_3P1 * sdc_3P1_1_unpol;

  // 3S1 octet
  const auto sdc_3S1_8_unpol = sdc::read_from_file(dataDir + "/3s1octetunpol.txt", sdcType);
  const auto sdc_3S1_8_long_psi = sdc::read_from_file(dataDir + "/3s1octetlong.txt", sdcType);
  // For calculating the longitudinal fraction of the 3S1 contributions, we
  // first assume PSI production to arrive at lth_3S1_psi. Then we calculate the
  // transmission from any 3S1 state to a chic1 and use the resulting lambda to
  // calculate the longitudinal fraction of the chi state (as observable after
  // its decay to a J/psi)
  const auto fL_3S1_psi = sdc_3S1_8_long_psi / sdc_3S1_8_unpol;
  const auto lth_3S1_psi = (1 - 3 * fL_3S1_psi) / (1 + fL_3S1_psi);
  // lth_chi = lth_psi / (3 + lth_psi)
  const auto lth_3S1_chi = (lth_3S1_psi) / (3 + lth_3S1_psi);
  const auto fL_3S1 = (1 - lth_3S1_chi) / (3 + lth_3S1_chi);
  const auto sdc_3S1_8_long = fL_3S1 * sdc_3S1_8_unpol;

  // Put everything into the combination MultiSDCs
  auto longitudinal = sdc::MultiSDC({sdc_3P1_1_long, sdc_3S1_8_long});
  longitudinal.scale_support(1. / ccbarMassSDCs);

  auto total = sdc::MultiSDC({sdc_3P1_1_unpol, sdc_3S1_8_unpol});
  total.scale_support(1. / ccbarMassSDCs);

  return {total, longitudinal};
}

/**
 * Read the total and longitudinal SDCs for the chic2. I.e. build a MultiSDC
 * that represents all the intermediate states that contribute to direct chi
 * production as observable after their decay to a J/psi since that is what the
 * fitting code uses (and also what the chi polarization measurements provide).
 *
 * The order of the contributions is (see als nrqcd_helpers.h):
 * - 3P2_1
 * - 3S1_8
 */
sdc::StateSDCs readChic2SDCs(std::string dataDir, sdc::SDCType sdcType=sdc::SDCType::LP_NLO) {
  // For further processing using the naming convention:
  // * S_U -> unpolarized SDC (== total)
  // * S_0 -> SDC of J_z = 0
  // * S_1 -> SDC of J_z = +/- 1
  // * S_2 -> SDC of J_z = +/- 2
  //
  // The calculations for the longitudinal SDCs follows the same calculations
  // (with different formulas) as the chic1 and the general approach is detailed
  // there. Here only the different formulas are stated where necessary

  // 3P2 octet
  const auto sdc_3P2_1_unpol = sdc::read_from_file(dataDir + "/3p2singletunpol.txt", sdcType);
  const auto sdc_3P2_1_h1 = sdc::read_from_file(dataDir + "/3p2singletpol1.txt", sdcType);
  const auto sdc_3P2_1_h2 = sdc::read_from_file(dataDir + "/3p2singletpol2.txt", sdcType);
  // We have S_total = S_U = S_0 + 2 * S_1 + 2 * S_2
  // --> S0 = S_U - 2 * (S_1 + S_2)
  const auto sdc_3P2_1_h0 = sdc_3P2_1_unpol - 2 * (sdc_3P2_1_h1 + sdc_3P2_1_h2);
  // lth = (-3 * S_0 - 3 * S_1 + 6 * S_2) / (5 * S_0 + 9 * S_1 + 6 * S_2)
  const auto lth_3P2 = (-3 * sdc_3P2_1_h0 - 3 * sdc_3P2_1_h1 + 6 * sdc_3P2_1_h2) /
      (5 * sdc_3P2_1_h0 + 9 * sdc_3P2_1_h1 + 6 * sdc_3P2_1_h2);
  const auto fL_3P2 = (1 - lth_3P2) / (3 + lth_3P2);
  const auto sdc_3P2_1_long = fL_3P2 * sdc_3P2_1_unpol;

  // 3S1 octet
  const auto sdc_3S1_8_unpol = sdc::read_from_file(dataDir + "/3s1octetunpol.txt", sdcType);
  const auto sdc_3S1_8_long_psi = sdc::read_from_file(dataDir + "/3s1octetlong.txt", sdcType);
  // See above for reasoning
  const auto fL_3S1_psi = sdc_3S1_8_long_psi / sdc_3S1_8_unpol;
  const auto lth_3S1_psi = (1 - 3 * fL_3S1_psi) / (1 + fL_3S1_psi);
  // lth_chi = 21 * lth_psi / (60 + 13 * lth_psi)
  const auto lth_3S1_chi = 21 * lth_3S1_psi / (60 + 13 * lth_3S1_psi);
  const auto fL_3S1 = (1 - lth_3S1_chi) / (3 + lth_3S1_chi);
  const auto sdc_3S1_8_long = fL_3S1 * sdc_3S1_8_unpol;

  // Put everything into the combination MultiSDCs
  auto longitudinal = sdc::MultiSDC({sdc_3P2_1_long, sdc_3S1_8_long});
  longitudinal.scale_support(1. / ccbarMassSDCs);

  auto total = sdc::MultiSDC({sdc_3P2_1_unpol, sdc_3S1_8_unpol});
  total.scale_support(1. / ccbarMassSDCs);

  return {total, longitudinal};
}

#endif
