#ifndef UTILS_HPP
#define UTILS_HPP

#include <fstream>
#include <iostream>
#include <iomanip>
#include "Math/Minimizer.h"
#include "./configuration.hpp"

using std::string;
using ROOT::Math::Minimizer;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;

void writeOutput(const string outputFilePath, Minimizer* min, const Configuration config);
void writeOutput(const string text, const string outputFilePath);

#endif
