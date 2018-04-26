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
using json = nlohmann::json;

using ROOT::Math::Minimizer;
using ROOT::Math::Factory;
using ROOT::Math::Functor;

// UBLAS
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
using ublas::matrix;
using ublas::prod;

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

matrix<double> readMatrixFromFile(const int &size, const string filePath) {
  matrix<double> m(size, size);
  ifstream file;
  file.open(filePath);
  if(!file) {
    std::cerr << "Unable to open file expMom.dat";
    exit(1);
  }

  double row = 0;
  double col = 0;
  double x;
  while (file >> x) {
    m(row, col) = x;
    if (col == size-1) {
      col = 0;
      row++;
    } else {
      col++;
    }
  }
  file.close();
  return m;
}

// following https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
bool almostEqualRelative(const double &a, const double &b) {
  // Machine epsilon https://en.wikipedia.org/wiki/Machine_epsilon
  double maxRelDiff = 1.0e-14;//2.22e-16;
  double diff = abs(a - b);
  double largest = (b > a) ? b : a;
  if ( diff <= largest * maxRelDiff )
    return true;
  return false;
}

bool compareVectors(const vector<double> &vec1, const vector<double> &vec2) {
  for(std::vector<double>::size_type i = 0; i != vec1.size(); i++) {
    if(!almostEqualRelative(vec1[i], vec2[i])) {
      cout << vec1[i] << " = " << vec2[i] << endl;
      return false;
    }
  }
  return true;
}

bool compareSquareMatrix(const matrix<double> &m1, const matrix<double> &m2) {
  for (auto i = 0; i < m1.size1(); i++) {
    for (auto j = 0; j < m1.size2(); j++) {
      if (!almostEqualRelative(m1(i, j), m2(i, j))) {
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

  // cout << "Are ExpMom same: " << compareVectors(readExpMomentsFromFile(), chisquared.expMom_.getExpPlusPionPoleMoments()) << endl;
  // cout << "Is CovMat same: " << compareSquareMatrix(readMatrixFromFile(9, "./data/covMat.dat"), chisquared.expMom_.covarianceMatrix) << endl;
  // Numerics num(constants);

  // cout << "invInvCovMat: " << prod(chisquared.expMom_.covarianceMatrix, chisquared.expMom_.inverseCovarianceMatrix) << endl;

  // beta 4th order
  // cout << "Matthias CHI2: " << chisquared(0.32307175541329564, 2.1e-2, -0.309486307083497, -3.0869411117237483e-2) << endl;
  // cout << "MY CHI2: " << chisquared(0.32306509503077424, 2.1e-2, -0.30906464837998099, -0.029781360856290042)  << endl;

  // beta 5th order
  cout << "Matthias CHI2: " << chisquared(0.32326096168471358, 2.1e-2, -0.31488720134123538, -2.6524803026353995e-2) << endl;
  cout << "MY CHI2: " << chisquared(0.31541087265066009, 2.1e-2, 0.81120226780379001, 5.2853783566688879)  << endl;


  // MINUIT
  Minimizer* min = Factory::CreateMinimizer("Minuit2", "Migrad");

  // set tolerances
  min->SetMaxFunctionCalls(10000000); // for Minuit2
  min->SetMaxIterations(10000000); // for GSL
  min->SetTolerance(1e-15);
  min->SetPrintLevel(1); // activate logging
  min->SetStrategy( 2 );

  // function wrapper
  Functor chi2(chisquared, 4);

  min->SetFunction(chi2);

  // set free variables to be minimized
  if (config["variables"]["astau"]["fixed"]) {
    min->SetFixedVariable(0, "astau", config["variables"]["astau"]["value"]);
  } else {
    min->SetVariable(0, "astau", config["variables"]["astau"]["value"], config["variables"]["astau"]["stepSize"]);
  }
  if (config["variables"]["aGGInv"]["fixed"]) {
    min->SetFixedVariable(1, "aGGInv", config["variables"]["aGGInv"]["value"]);
  } else {
    min->SetVariable(1, "aGGInv", config["variables"]["aGGInv"]["value"], config["variables"]["aGGInv"]["stepSize"]);
  }
  if (config["variables"]["rhoVpA"]["fixed"]) {
    min->SetFixedVariable(2, "rhoVpA", config["variables"]["rhoVpA"]["value"]);
  } else {
    min->SetVariable(2, "rhoVpA", config["variables"]["rhoVpA"]["value"], config["variables"]["rhoVpA"]["stepSize"]);
  }
  if (config["variables"]["c8VpA"]["fixed"]) {
    min->SetFixedVariable(3, "c8VpA", config["variables"]["c8VpA"]["value"]);
  } else {
    min->SetVariable(3, "c8VpA", config["variables"]["c8VpA"]["value"], config["variables"]["c8VpA"]["stepSize"]);
  }

  // minimize!
  min->Minimize();

  // min->PrintResults();
  return 0;
}
