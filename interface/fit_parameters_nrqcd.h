#ifndef FIT_PARAMETER_NRQCD_H
#define FIT_PARAMETER_NRQCD_H

#include "simple_compile_time_map.h"

#include "TTree.h"

constexpr std::array<ParameterIndex, 23> NRQCD_PARAMETERS = {{
    // place the ones that do not change on top for easier maintenance and less
    // repetition below
    // Branching fraction nuisance parameters
    {"br_psip_dp", 0},
    {"br_psip_mm", 1},
    {"br_psip_c2", 2},
    {"br_psip_c1", 3},
    {"br_psip_jpsi", 4},
    {"br_c2_jpsi", 5},
    {"br_c1_jpsi", 6},
    {"br_jpsi_mm", 7},

    // Luminosity nuisance parameters
    {"L_CMS", 8},
    {"L_ATLAS", 9},

    // (irrelevant) Costh ratio normalizations
    {"norm_costh_1", 10},
    {"norm_costh_2", 11},
    {"norm_costh_3", 12},

    // Actually interesting LDMEs
    // chic0 LDMEs are free parameters, chic1 and chic2 follow from them
    {"l_3S1_8_c0", 13},
    {"l_3P0_1_c0", 14},

    // psi and jpsi LDMEs
    // Singlets: constrained by theory, independent for the two states
    {"l_3S1_1_jpsi", 15},
    {"l_3S1_1_psip", 16},

    // For the octets, we want to potentially constrain parameters between the
    // two states, resp. we want to at least constrain ratios of parameters to
    // be the same between the two states. Hence, they are defined as follows:
    // - 1S0_8: free parameters for both states
    {"l_1S0_8_jpsi", 17},
    {"l_1S0_8_psip", 18},
    // The other two octets are related to this one for the jpsi
    // - 3PJ_8 = r_3PJ_8_1S0_8 * 1S0_8 and r_3PJ_8_1S0_8 (the ratio between the
    //   two) is the free parameter (for both states independently)
    // - 3S1_8 = r_3S1_8_1S0_8 * 1S0_8, where again r_3S1_8_1S0_8 is free
    {"l_r_3PJ_8_1S0_8_jpsi", 19},
    {"l_r_3S1_8_1S0_8_jpsi", 20},
    // And finally these same octets for the psip are related to this ratios via
    // double ratios:
    //
    // - 3PJ_8 = 1S0_8_psip * r_3PJ_8_1S0_8_jpsi * l_rr_3PJ_8_1S0_8_psip_jpsi,
    //   where the latter = r_3PJ_8_1S0_8_psip / r_3PJ_8_1S0_8_jpsi. I.e. the
    //   double ratio of the above defined ratios for the jpsi, dividing the psip
    //   state by the jpsi state
    // - 3S1_8 = 1S0_8_psip * r_3S1_8_1S0_8_jpsi * l_rr_3S1_8_1S0_8_psip_jpsi, in
    //   the same fashion as the 3PJ
    {"l_rr_3PJ_8_1S0_8_psip_jpsi", 21},
    {"l_rr_3S1_8_1S0_8_psip_jpsi", 22},
}};

constexpr static int IPAR(const char* name) { return getParIdx(NRQCD_PARAMETERS, name); }
constexpr static int NPARS(int max = 0, int index = 0) {
  return (index == NRQCD_PARAMETERS.size()) ? max + 1 : // +1 due to 0-indexing
      NRQCD_PARAMETERS[index].second > max ? NPARS(NRQCD_PARAMETERS[index].second, index + 1) : NPARS(max, index + 1);
}

/**
 * Get the first name corresponding to the passed parameter index
 */
static constexpr const char* PARNAME(const int index)
{
  return getParName(NRQCD_PARAMETERS, index);
}

template <size_t N = NPARS()>
void parametersIndicesAsTTree(TTree* tree, const std::array<ParameterIndex, N>& parameters = NRQCD_PARAMETERS) {
  std::array<int, N> indices;
  int i = 0;
  for (const auto& par : parameters) {
    indices[i] = par.second;
    tree->Branch(par.first, &indices[i]);
    i++;
  }
  tree->Fill();
}

#endif
