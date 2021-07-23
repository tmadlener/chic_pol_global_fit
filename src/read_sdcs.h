#ifndef GLOBALFIT_READSDCS_H
#define GLOBALFIT_READSDCS_H

#include "sdc.h"

#include <string>
#include <utility>

constexpr static double ccbarMassSDCs = 3.0; // GeV
constexpr static double cQuarkMass = 1.5; // GeV
constexpr static double cQuarkMass2 = cQuarkMass * cQuarkMass; // GeV^2

sdc::SDCType asOrderEnum(const std::string& sdcOrder) {
  const auto it = std::find(sdc::SDCTypeNames.begin(), sdc::SDCTypeNames.end(), sdcOrder);
  if (it == sdc::SDCTypeNames.end()) {
    std::cerr << "Argument to --order is not valid. Must be one of \'LO\', \'NLO\' or \'LP+NLO\'" << std::endl;
    std::exit(1);
  }

  return util::to_enum<sdc::SDCType>(std::distance(sdc::SDCTypeNames.begin(), it));
}

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

  // 3PJ octet (multiplied by m_c^2 to have LDMEs in same units as S-wave states)
  const auto sdc_3PJ_8_unpol = cQuarkMass2 * sdc::read_from_file(dataDir + "/3pjoctetunpol.txt", sdcType);
  const auto sdc_3PJ_8_long = cQuarkMass2 * sdc::read_from_file(dataDir + "/3pjoctetlong.txt", sdcType);

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

  const auto inputfile_3P1 = [dataDir, sdcType]() {
    if (sdcType == sdc::SDCType::NLO) {
      return dataDir + "/3p1singlet_NLO.txt";
    }
    else if (sdcType == sdc::SDCType::LP_NLO) {
      return dataDir + "/3p1singlet_NLOLP.txt";
    }
    return std::string("");
  }();

  // 3P1 singlet (re-use the sdc types to get the right column)
  const auto sdc_3P1_1_h1 = cQuarkMass2 * sdc::read_from_file(inputfile_3P1, util::to_enum<sdc::SDCType>(1));
  const auto sdc_3P1_1_h0 = cQuarkMass2 * sdc::read_from_file(inputfile_3P1, util::to_enum<sdc::SDCType>(0));

  const auto sdc_3P1_1_unpol = sdc_3P1_1_h0 + 2 * sdc_3P1_1_h1;

  // // 3P1 singlet (multiplied by m_c^2 to have LDMEs in same units as S-wave states)
  // const auto sdc_3P1_1_unpol = cQuarkMass2 * sdc::read_from_file(dataDir + "/3p1singletunpol.txt", sdcType);
  // const auto sdc_3P1_1_h1 = cQuarkMass2 * sdc::read_from_file(dataDir + "/3P1_1_h1_Carlos_interpolation.txt", sdcType);
  // // We have S_total = S_U = S_0 + 2 * S_1
  // // --> S_0 = S_U - 2 * S_1
  // const auto sdc_3P1_1_h0 = sdc_3P1_1_unpol - 2 * sdc_3P1_1_h1;

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
  // lth_chi = lth_psi / (4 + lth_psi)
  const auto lth_3S1_chi = (lth_3S1_psi) / (4 + lth_3S1_psi);
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
   const auto inputfile_3P2 = [dataDir, sdcType]() {
    if (sdcType == sdc::SDCType::NLO) {
      return dataDir + "/3p2singlet_NLO.txt";
    }
    else if (sdcType == sdc::SDCType::LP_NLO) {
      return dataDir + "/3p2singlet_NLOLP.txt";
    }
    return std::string("");
  }();

  // The calculations for the longitudinal SDCs follows the same calculations
  // (with different formulas) as the chic1 and the general approach is detailed
  // there. Here only the different formulas are stated where necessary

 // 3P1 singlet (re-use the sdc types to get the right column)
  const auto sdc_3P2_1_h2 = cQuarkMass2 * sdc::read_from_file(inputfile_3P2, util::to_enum<sdc::SDCType>(2));
  const auto sdc_3P2_1_h1 = cQuarkMass2 * sdc::read_from_file(inputfile_3P2, util::to_enum<sdc::SDCType>(1));
  const auto sdc_3P2_1_h0 = cQuarkMass2 * sdc::read_from_file(inputfile_3P2, util::to_enum<sdc::SDCType>(0));

  const auto sdc_3P2_1_unpol = sdc_3P2_1_h0 + 2 * sdc_3P2_1_h1 + 2 * sdc_3P2_1_h2;

  // 3P2 singlet (multiplied by m_c^2 to have LDMEs in same units as S-waves states)
  // const auto sdc_3P2_1_unpol = cQuarkMass2 * sdc::read_from_file(dataDir + "/3p2singletunpol.txt", sdcType);
  // // const auto sdc_3P2_1_h1 = cQuarkMass2 * sdc::read_from_file(dataDir + "/3P2_1_h1_Carlos_interpolation.txt", sdcType);
  // const auto sdc_3P2_1_h0 = cQuarkMass2 * sdc::read_from_file(dataDir + "/3P2_1_h0_Carlos_interpolation.txt", sdcType);
  // const auto sdc_3P2_1_h2 = cQuarkMass2 * sdc::read_from_file(dataDir + "/3P2_1_h2_Carlos_interpolation.txt", sdcType);
  // // We have S_total = S_U = S_0 + 2 * S_1 + 2 * S_2
  // // --> S0 = S_U - 2 * (S_1 + S_2)
  // // const auto sdc_3P2_1_h0 = sdc_3P2_1_unpol - 2 * (sdc_3P2_1_h1 + sdc_3P2_1_h2);
  // const auto sdc_3P2_1_h1 = 0.5 * (sdc_3P2_1_unpol - sdc_3P2_1_h0 - 2 * sdc_3P2_1_h2);

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


sdc::SDC getCombinedSDC(std::string const& dir, std::string const& sdc, double lpFraction=0) {
  const std::string nloDir = dir + "/NLO/";
  const std::string lpDir = dir + "/LP/";

  // all files only have one column
  const auto nlo = sdc::read_from_file(nloDir + sdc + ".txt", 1);
  const auto lp = sdc::read_from_file(lpDir + sdc + ".txt", 1);

  return nlo + lpFraction * lp;
}

/**
 * Read the total and longitudinal SDCs of the psi states using the NLO as
 * baseline and adding a fraction of the LP corrections on top of them.
 *
 * Returns a MultiSDC where the contributions have the following order (see also nrqcd_helpers.h)
 * - 3S1_1
 * - 3S1_8
 * - 3PJ_8
 * - 1S0_8
 */
sdc::StateSDCs readPsiSDCs(std::string const& sdcDir, double lpFraction=0) {
  // 3S1 singlet
  const auto s_3S1_1_total = getCombinedSDC(sdcDir, "3S1_1_total", lpFraction);
  const auto s_3S1_1_trans = getCombinedSDC(sdcDir, "3S1_1_trans", lpFraction);
  // S_L = S_total - 2 * S_T
  const auto s_3S1_1_long = s_3S1_1_total - 2 * s_3S1_1_trans;

  // 3S1 octet
  const auto s_3S1_8_total = getCombinedSDC(sdcDir, "3S1_8_total", lpFraction);
  const auto s_3S1_8_long = getCombinedSDC(sdcDir, "3S1_8_long", lpFraction);

  // 3PJ octet (multiplied by m_c^2 to have LDMEs in same units as S-wave states)
  const auto s_3PJ_8_total = cQuarkMass2 * getCombinedSDC(sdcDir, "3PJ_8_total", lpFraction);
  const auto s_3PJ_8_long = cQuarkMass2 * getCombinedSDC(sdcDir, "3PJ_8_long", lpFraction);

  // 1S0 octet
  const auto s_1S0_8 = getCombinedSDC(sdcDir, "1S0_8_total", lpFraction);
  // the octet is unpolarized, hence the three components 0, +1, -1 are equal
  // -> long = 1/3 * total
  const auto s_1S0_8_long = 1. / 3 * s_1S0_8;

  // Put everything into the combination MultiSDCs
  auto longitudinal = sdc::MultiSDC({s_3S1_1_long, s_3S1_8_long, s_3PJ_8_long, s_1S0_8_long});
  longitudinal.scale_support(1. / ccbarMassSDCs); // we want them to be in pT/M

  auto total = sdc::MultiSDC({s_3S1_1_total, s_3S1_8_total, s_3PJ_8_total, s_1S0_8});
  total.scale_support(1. / ccbarMassSDCs); // we want them to be in pT/M

  return {total, longitudinal};
}


/**
 * Read the total and longitudinal SDCs for the chic1. I.e. build a MultiSDC
 * that represents all the intermediate states that contribute to direct chi
 * production as observable after their decay to a J/psi since that is what the
 * fitting code uses (and also what the chi polarization measurements provide).
 * To calculate the SDCs use the NLO as baseline and add a fraction of the LP
 * corrections on top of that
 *
 * The order of the contributions is (see als nrqcd_helpers.h):
 * - 3P1_1
 * - 3S1_8
 */
sdc::StateSDCs readChic1SDCs(std::string const& sdcDir, double lpFraction=0) {
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
  const auto s_3P1_1_h1 = cQuarkMass2 * getCombinedSDC(sdcDir, "3P1_1_h1", lpFraction);
  const auto s_3P1_1_h0 = cQuarkMass2 * getCombinedSDC(sdcDir, "3P1_1_h0", lpFraction);

  const auto s_3P1_1_total = s_3P1_1_h0 + 2 * s_3P1_1_h1;
  // lth = (3P1_1_h0 - 3P1_1_h1) / (3P1_1_h0 + 3 * 3P1_1_h1)
  const auto lth_3P1 = (s_3P1_1_h0 - s_3P1_1_h1) / (s_3P1_1_h0 + 3 * s_3P1_1_h1);
  const auto fL_3P1 = (1 - lth_3P1) / (3 + lth_3P1);
  const auto s_3P1_1_long = fL_3P1 * s_3P1_1_total;

  // 3S1 octet
  const auto s_3S1_8_total = getCombinedSDC(sdcDir, "3S1_8_total", lpFraction);
  const auto s_3S1_8_long_psi = getCombinedSDC(sdcDir, "3S1_8_long", lpFraction);

  // For calculating the longitudinal fraction of the 3S1 contributions, we
  // first assume PSI production to arrive at lth_3S1_psi. Then we calculate the
  // transmission from any 3S1 state to a chic1 and use the resulting lambda to
  // calculate the longitudinal fraction of the chi state (as observable after
  // its decay to a J/psi)
  const auto fL_3S1_psi = s_3S1_8_long_psi / s_3S1_8_total;
  const auto lth_3S1_psi = (1 - 3 * fL_3S1_psi) / (1 + fL_3S1_psi);
  // lth_chi = lth_psi / (4 + lth_psi)
  const auto lth_3S1_chi = (lth_3S1_psi) / (4 + lth_3S1_psi);
  const auto fL_3S1 = (1 - lth_3S1_chi) / (3 + lth_3S1_chi);
  const auto s_3S1_8_long = fL_3S1 * s_3S1_8_total;

  // Put everything into the combination MultiSDCs
  auto longitudinal = sdc::MultiSDC({s_3P1_1_long, s_3S1_8_long});
  longitudinal.scale_support(1. / ccbarMassSDCs);

  auto total = sdc::MultiSDC({s_3P1_1_total, s_3S1_8_total});
  total.scale_support(1. / ccbarMassSDCs);

  return {total, longitudinal};
}

/**
 * Read the total and longitudinal SDCs for the chic2. I.e. build a MultiSDC
 * that represents all the intermediate states that contribute to direct chi
 * production as observable after their decay to a J/psi since that is what the
 * fitting code uses (and also what the chi polarization measurements provide).
 * To calculate the SDCs use the NLO as baseline and add a fraction of the LP
 * corrections on top of that
 *
 * The order of the contributions is (see als nrqcd_helpers.h):
 * - 3P2_1
 * - 3S1_8
 */
sdc::StateSDCs readChic2SDCs(std::string sdcDir, double lpFraction=0) {
  // For further processing using the naming convention:
  // * S_U -> unpolarized SDC (== total)
  // * S_0 -> SDC of J_z = 0
  // * S_1 -> SDC of J_z = +/- 1
  // * S_2 -> SDC of J_z = +/- 2

  // The calculations for the longitudinal SDCs follows the same calculations
  // (with different formulas) as the chic1 and the general approach is detailed
  // there. Here only the different formulas are stated where necessary

  // 3P2 singlet
  const auto s_3P2_1_h0 = cQuarkMass2 * getCombinedSDC(sdcDir, "3P2_1_h0", lpFraction);
  const auto s_3P2_1_h1 = cQuarkMass2 * getCombinedSDC(sdcDir, "3P2_1_h1", lpFraction);
  const auto s_3P2_1_h2 = cQuarkMass2 * getCombinedSDC(sdcDir, "3P2_1_h2", lpFraction);

  const auto s_3P2_1_total = s_3P2_1_h0 + 2 * s_3P2_1_h1 + 2 * s_3P2_1_h2;
  // lth = (-3 * S_0 - 3 * S_1 + 6 * S_2) / (5 * S_0 + 9 * S_1 + 6 * S_2)
  const auto lth_3P2 = (-3 * s_3P2_1_h0 - 3 * s_3P2_1_h1 + 6 * s_3P2_1_h2) /
      (5 * s_3P2_1_h0 + 9 * s_3P2_1_h1 + 6 * s_3P2_1_h2);
  const auto fL_3P2 = (1 - lth_3P2) / (3 + lth_3P2);
  const auto s_3P2_1_long = fL_3P2 * s_3P2_1_total;

  // 3S1 octet
  const auto s_3S1_8_total = getCombinedSDC(sdcDir, "3S1_8_total", lpFraction);
  const auto s_3S1_8_long_psi = getCombinedSDC(sdcDir, "3S1_8_long", lpFraction);

  // See above for reasoning
  const auto fL_3S1_psi = s_3S1_8_long_psi / s_3S1_8_total;
  const auto lth_3S1_psi = (1 - 3 * fL_3S1_psi) / (1 + fL_3S1_psi);
  // lth_chi = 21 * lth_psi / (60 + 13 * lth_psi)
  const auto lth_3S1_chi = 21 * lth_3S1_psi / (60 + 13 * lth_3S1_psi);
  const auto fL_3S1 = (1 - lth_3S1_chi) / (3 + lth_3S1_chi);
  const auto s_3S1_8_long = fL_3S1 * s_3S1_8_total;

  // Put everything into the combination MultiSDCs
  auto longitudinal = sdc::MultiSDC({s_3P2_1_long, s_3S1_8_long});
  longitudinal.scale_support(1. / ccbarMassSDCs);

  auto total = sdc::MultiSDC({s_3P2_1_total, s_3S1_8_total});
  total.scale_support(1. / ccbarMassSDCs);

  return {total, longitudinal};
}

#endif
