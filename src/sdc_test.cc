#include "sdc.h"

#include <iostream>

int main() {
  const auto sdc_1S0_8 = sdc::read_from_file("../SDCs_Chung/1s0octetunpol.txt");

  const auto sdc = sdc::MultiSDC({sdc_1S0_8, sdc_1S0_8});

  std::cout << sdc_1S0_8 << std::endl;

  std::cout << sdc_1S0_8(10.5) << " " << sdc_1S0_8(10.5, sdc::util::log_interpolate) << std::endl;

  std::cout << sdc << std::endl;

  std::cout << sdc(95, {1.2, 3.4}, sdc::util::log_interpolate) <<std::endl;

  return 0;
}
