#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
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

void writeOutput(const string configFilePath, const string outputFilePath, Minimizer* min);

#endif
