#ifndef SRC_UTILS_HPP
#define SRC_UTILS_HPP

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include "Math/Minimizer.h"
#include "./configuration.hpp"
#include "./theoretical_moments.hpp"
#include "./weights.hpp"

using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;

void writeOutput(
  const std::string outputFilePath,
  const ROOT::Math::Minimizer* min,
  const Configuration config
);
void writeOutput(const string text, const string outputFilePath);

namespace utils {
  int momCount(std::vector<Input> inputs);
}

#endif
