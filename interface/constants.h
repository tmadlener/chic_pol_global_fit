#ifndef H_CONSTANTS__
#define H_CONSTANTS__

#include <array>

// Some PDG branching fractions and their uncertainties from the 2020 PDG
constexpr std::array<double, 2> B_PSIP_MM = {8.0e-3, 0.6e-3};
constexpr std::array<double, 2> B_PSIP_CHIC2 = {9.52e-2, 0.20e-2};
constexpr std::array<double, 2> B_PSIP_CHIC1 = {9.75e-2, 0.24e-2};
constexpr std::array<double, 2> B_PSIP_JPSI = {61.4e-2, 0.6e-2};
constexpr std::array<double, 2> B_PSIP_PIPI = {34.68e-2, 0.3e-2};
constexpr std::array<double, 2> B_CHIC2_JPSI = {19.0e-2, 0.5e-2};
constexpr std::array<double, 2> B_CHIC1_JPSI = {34.3e-2, 1.0e-2};
constexpr std::array<double, 2> B_JPSI_MM = {5.961e-2, 0.033e-2};

// masses in GeV
const double M_JPSI = 3.097;
const double M_PSI2S = 3.686;
const double M_CHIC2 = 3.556;
const double M_CHIC1 = 3.511;

const double PTMNORM = 5.;

#endif
