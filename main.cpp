#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <stdexcept>

using std::vector;
using std::cout;
using std::endl;
using std::complex;
using std::string;
using std::runtime_error;

#include "./src/constants.hpp"
// Theoretical Moments
#include "./src/adler_function.hpp"
#include "./src/alpha_s.hpp"
#include "./src/numerics.hpp"

// Chisquared
#include "./src/chisquared.hpp"
#include "./src/s0_sets.hpp"
#include "./src/weights.hpp"

// MINUIT
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"

// Config
#include "json.hpp"
#include <fstream>
using json = nlohmann::json;

using ROOT::Math::Minimizer;
using ROOT::Math::Factory;
using ROOT::Math::Functor;

// UBLAS
#include <boost/numeric/ublas/matrix.hpp>
using ublas::matrix;

using std::ifstream;

vector<double> readExpMomentsFromFile() {
  vector<double> expMoms;
  std::ifstream expMomFile;
  expMomFile.open("./data/expMom.dat");
  if(!expMomFile) {
    std::cerr << "Unable to open file expMom.dat";
    exit(1);
  }

  double x;
  while (expMomFile >> x) {
    expMoms.push_back(x);
  }
  expMomFile.close();
  return expMoms;
}

matrix<double> readCovMatrixFromFile(const int &size) {
  matrix<double> covMat(size, size);
  ifstream covMatFile;
  covMatFile.open("./data/covMat.dat");
  if(!covMatFile) {
    std::cerr << "Unable to open file expMom.dat";
    exit(1);
  }

  double row = 0;
  double col = 0;
  double x;
  while (covMatFile >> x) {
    covMat(row, col) = x;
    if (col == size-1) {
      col = 0;
      row++;
    } else {
      col++;
    }
  }
  covMatFile.close();
  return covMat;
}

bool areSame(double a, double b) {
  return fabs(a - b) < 1.0e-14;
}

bool compareVectors(const vector<double> &vec1, const vector<double> &vec2) {
  for(std::vector<double>::size_type i = 0; i != vec1.size(); i++) {
    if(!areSame(vec1[i], vec2[i])) {
      cout << vec1[i] << " = " << vec2[i] << endl;
      return false;
    }
  }
  return true;
}

bool compareSquareMatrix(const matrix<double> &m1, const matrix<double> &m2) {
  for (auto i = 0; i < m1.size1(); i++) {
    cout << endl;
    for (auto j = 0; j < m1.size2(); j++) {
      if (!areSame(m1(i, j), m2(i, j))) {
        cout << "at m(" << i << ", " << j << ") \t" << m1(i, j) << " != " << m2(i, j) << endl;
        return false;
      }
    }
  }
  return true;
}

int main () {
  cout.precision(17);

  std::ifstream configFile("./configuration.json");
  json config;
  configFile >> config;

  const Constants constants(config);
  Chisquared chisquared(config, constants);

  cout << "Are same: " << compareVectors(readExpMomentsFromFile(), chisquared.expMom_.getExpPlusPionPoleMoments()) << endl;
  cout << "Are CovMat same: " << compareSquareMatrix(readCovMatrixFromFile(9), chisquared.expMom_.covarianceMatrix) << endl;


  chisquared.log(0.32307175541329564, 2.1e-2, -0.309486307083497, -3.0869411117237483e-2);

  // MINUIT
  // Minimizer* min = Factory::CreateMinimizer("Minuit2", "Migrad");

  // // set tolerances
  // min->SetMaxFunctionCalls(10000000); // for Minuit2
  // min->SetMaxIterations(10000000); // for GSL
  // min->SetTolerance(1e-15);
  // min->SetPrintLevel(1); // activate logging

  // // function wrapper
  // Functor chi2(chisquared, 4);

  // min->SetFunction(chi2);

  // // set free variables to be minimized
  // if (config["variables"]["astau"]["fixed"]) {
  //   min->SetFixedVariable(0, "astau", config["variables"]["astau"]["value"]);
  // } else {
  //   min->SetVariable(0, "astau", config["variables"]["astau"]["value"], config["variables"]["astau"]["stepSize"]);
  // }
  // if (config["variables"]["aGGInv"]["fixed"]) {
  //   min->SetFixedVariable(1, "aGGInv", config["variables"]["aGGInv"]["value"]);
  // } else {
  //   min->SetVariable(1, "aGGInv", config["variables"]["aGGInv"]["value"], config["variables"]["aGGInv"]["stepSize"]);
  // }
  // if (config["variables"]["rhoVpA"]["fixed"]) {
  //   min->SetFixedVariable(2, "rhoVpA", config["variables"]["rhoVpA"]["value"]);
  // } else {
  //   min->SetVariable(2, "rhoVpA", config["variables"]["rhoVpA"]["value"], config["variables"]["rhoVpA"]["stepSize"]);
  // }
  // if (config["variables"]["c8VpA"]["fixed"]) {
  //   min->SetFixedVariable(3, "c8VpA", config["variables"]["c8VpA"]["value"]);
  // } else {
  //   min->SetVariable(3, "c8VpA", config["variables"]["c8VpA"]["value"], config["variables"]["c8VpA"]["stepSize"]);
  // }

  // // minimize!
  // min->Minimize();

  // min->PrintResults();
  return 0;
}
