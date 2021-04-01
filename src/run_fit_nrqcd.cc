#include "ArgParser.h"
#include "global_fit_nrqcd.C"

#include <iostream>

int main(int argc, char* argv[]) {
  const auto parser = ArgParser(argc, argv);

  global_fit_nrqcd();

  return 0;
}
