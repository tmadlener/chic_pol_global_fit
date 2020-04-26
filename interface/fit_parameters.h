#ifndef FIT_PARAMETER_H__
#define FIT_PARAMETER_H__

#include "simple_compile_time_map.h"

// make it possible to switch the parametrization of the cross section
// functions:
// 0 - Usual fit with 6 independent betas. The fit with the most freedom.
// 1 - Fit with only 1 beta. This implies constant polarizations and feed-down
// fractions
// 2 - Fit with same betas for psi, one transverse beta and two independent
// longitudinal betas for the chi. Implies constant psi polarization, but allows
// for smooth changes of the chi polarizations and feed-downs

#ifndef FIT_OPTION
#define FIT_OPTION 0
#endif

#if FIT_OPTION > 2
#error "FIT_OPTION can only be 0, 1 or 2"
#endif

constexpr std::array<ParameterIndex, 27> PARAMETERS = {{
    {"sigma_psip",  0},
    {"sigma_chic2", 1},
    {"sigma_chic1", 2},
    {"sigma_jpsi", 3},

    {"f_long_psi", 4},
    {"f_long_c1", 5},
    {"f_long_c2", 6},

    {"gamma", 7},
    {"beta_long_psi", 8},
#if FIT_OPTION == 0
    {"beta_trans_psi", 9},
    {"beta_long_c1", 10},
    {"beta_trans_c1", 11},
    {"beta_long_c2", 12},
    {"beta_trans_c2", 13},

    {"br_psip_dp", 14},
    {"br_psip_mm", 15},
    {"br_psip_c2", 16},
    {"br_psip_c1", 17},
    {"br_psip_jpsi", 18},
    {"br_c2_jpsi", 19},
    {"br_c1_jpsi", 20},
    {"br_jpsi_mm", 21},

    {"L_CMS", 22},
    {"L_ATLAS", 23},

    {"norm_costh_1", 24},
    {"norm_costh_2", 25},
    {"norm_costh_3", 26},
#elif FIT_OPTION == 1
    {"beta_trans_psi", 8},
    {"beta_long_c1", 8},
    {"beta_trans_c1", 8},
    {"beta_long_c2", 8},
    {"beta_trans_c2", 8},

    {"br_psip_dp", 9},
    {"br_psip_mm", 10},
    {"br_psip_c2", 11},
    {"br_psip_c1", 12},
    {"br_psip_jpsi", 13},
    {"br_c2_jpsi", 14},
    {"br_c1_jpsi", 15},
    {"br_jpsi_mm", 16},

    {"L_CMS", 17},
    {"L_ATLAS", 18},

    {"norm_costh_1", 19},
    {"norm_costh_2", 20},
    {"norm_costh_3", 21},
#else // now only 2 is still available
    {"beta_trans_psi", 8},
    {"beta_long_c1", 9},
    {"beta_trans_c1", 10},
    {"beta_long_c2", 11},
    {"beta_trans_c2", 10},

    {"br_psip_dp", 12},
    {"br_psip_mm", 13},
    {"br_psip_c2", 14},
    {"br_psip_c1", 14},
    {"br_psip_jpsi", 15},
    {"br_c2_jpsi", 16},
    {"br_c1_jpsi", 17},
    {"br_jpsi_mm", 18},

    {"L_CMS", 19},
    {"L_ATLAS", 20},

    {"norm_costh_1", 21},
    {"norm_costh_2", 22},
    {"norm_costh_3", 23},
#endif
  }};

/**
 * Convenience overload
 */
static constexpr int IPAR(const char* name)
{
  return getParIdx(PARAMETERS, name);
}

/**
 * Get the highest index of all defined parameters
 */
static constexpr int NPARS(int max=0, int index=0)
{
  return (index == PARAMETERS.size()) ? max + 1: // +1 due to 0-indexing
    PARAMETERS[index].second > max ? NPARS(PARAMETERS[index].second, index + 1) :
    NPARS(max, index + 1);
}

/**
 * Get the first name corresponding to the passed parameter index
 */
static constexpr const char* PARNAME(const int index)
{
  return getParName(PARAMETERS, index);
}


#endif
