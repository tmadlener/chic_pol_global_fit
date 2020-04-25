#ifndef FIT_PARAMETER_H__
#define FIT_PARAMETER_H__

#include "simple_compile_time_map.h"

#ifndef USE_SAME_BETAS
#define USE_SAME_BETAS 0
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
#if USE_SAME_BETAS == 0
    {"beta_long_psi", 8},
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
#else
    {"beta_long_psi", 8},
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
 * Get the number of actually defined parameters
 */
static constexpr int NPARS()
{
  return PARAMETERS.back().second + 1;
}

/**
 * Get the first name corresponding to the passed parameter index
 */
static constexpr const char* PARNAME(const int index)
{
  return getParName(PARAMETERS, index);
}

#endif
