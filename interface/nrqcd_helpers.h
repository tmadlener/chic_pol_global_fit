#ifndef GLOBALFIT_NRQCD_HELPERS_H
#define GLOBALFIT_NRQCD_HELPERS_H

#include <array>

/**
 * The SDCs that contribute to the psi states and their indices in the MultiSDC
 */
enum class PsiSDCs {
  s3S1_1 = 0, ///< 3S1 singlet
  s3S1_8 = 1, ///< 3S1 octet
  s3PJ_8 = 2, ///< 3PJ octet
  s1S0_8 = 3, ///< 1S0 octet
};

constexpr static std::array<const char*, 4> PsiSDCsNames = {"3S1_1", "3S1_8", "3PJ_8", "1S0_8"};
constexpr static std::array<PsiSDCs, 4> allPsiSDCs = {PsiSDCs::s3S1_1, PsiSDCs::s3S1_8, PsiSDCs::s3PJ_8,
                                                      PsiSDCs::s1S0_8};

/**
 * The SDCs that contribute to the chic1 states and their indices in the Multi SDCs
 */
enum class Chic1SDCs {
  s3P1_1 = 0, ///< 3P1 singlet
  s3S1_8 = 1, ///< 3S1 octet
};

constexpr static std::array<const char*, 2> Chic1SDCsNames = {"3P1_1", "3S1_8"};
constexpr static std::array<Chic1SDCs, 2> allChic1SDCs = {Chic1SDCs::s3P1_1, Chic1SDCs::s3S1_8};

/**
 * The SDCs that contribute to the chic2 states and their indices in the Multi SDCs
 */
enum class Chic2SDCs {
  s3P2_1 = 0, ///< 3P2 singlet
  s3S1_8 = 1, ///< 3S1 octet
};

constexpr static std::array<const char*, 2> Chic2SDCsNames = {"3P2_1", "3S1_8"};
constexpr static std::array<Chic2SDCs, 2> allChic2SDCs = {Chic2SDCs::s3P2_1, Chic2SDCs::s3S1_8};

#endif
