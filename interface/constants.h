#ifndef H_CONSTANTS__
#define H_CONSTANTS__

#include <array>

// Some PDG branching fractions and their uncertainties from the 2015 PDG
constexpr std::array<double, 2> B_PSIP_MM = {7.9e-3, 0.9e-3};
constexpr std::array<double, 2> B_PSIP_CHIC2 = {9.11e-2, 0.31e-2};
constexpr std::array<double, 2> B_PSIP_CHIC1 = {9.55e-2, 0.31e-2};
constexpr std::array<double, 2> B_PSIP_JPSI = {60.9e-2, 0.6e-2};
constexpr std::array<double, 2> B_PSIP_PIPI = {34.46e-2, 0.3e-2};
constexpr std::array<double, 2> B_CHIC2_JPSI = {19.2e-2, 0.7e-2};
constexpr std::array<double, 2> B_CHIC1_JPSI = {33.9e-2, 1.2e-2};
constexpr std::array<double, 2> B_JPSI_MM = {5.961e-2, 0.033e-2};

// masses in GeV
const double M_JPSI = 3.097;
const double M_PSI2S = 3.686;
const double M_CHIC2 = 3.556;
const double M_CHIC1 = 3.511;

const double MIN_PTM = 2;
const double PTMNORM = 2.;

#endif
