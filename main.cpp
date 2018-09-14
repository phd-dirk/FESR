#define BOOST_UBLAS_NDEBUG 1 // numerical issues with ublas matrix inverse

#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <string>
#include <stdexcept>
#include "./src/types.hpp"

using std::vector;
using std::cout;
using std::endl;
using std::complex;
using std::string;
using std::runtime_error;

#include "./src/configuration.hpp"
// Theoretical Moments
#include "./src/adler_function.hpp"
#include "./src/alpha_s.hpp"
#include "./src/numerics.hpp"
#include "./src/duality_violations.hpp"

// Chisquared
#include "./src/chisquared.hpp"
#include "./src/weights.hpp"

// MINUIT
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"

using ROOT::Math::Minimizer;
using ROOT::Math::Factory;
using ROOT::Math::Functor;

// UBLAS
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
using ublas::matrix;
using ublas::prod;

using std::ifstream;
using std::ofstream;

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


int dof(const json config) {
  int numS0s = config["parameters"]["s0Set"].size();
  int numVar = 0;
  if (!config["variables"]["astau"]["fixed"]) {
    numVar++;
  }
  if (!config["variables"]["aGGInv"]["fixed"]) {
    numVar++;
  }
  if (!config["variables"]["rhoVpA"]["fixed"]) {
    numVar++;
  }
  if (!config["variables"]["c8VpA"]["fixed"]) {
    numVar++;
  }
  return numS0s-numVar;
}

int dof(const Configuration config) {
  int dof = 1;
  // int dof = config.s0Set.size();
  // if(!config.astau.isFixed)
  //   dof--;
  // if(!config.aGGInv.isFixed)
  //   dof--;
  // if(!config.rhoVpA.isFixed)
  //   dof--;
  // if(!config.c8VpA.isFixed)
  //   dof--;
  return dof;
}

void writeOutput(const string configFilePath, const string outputFilePath, Minimizer* min) {
  ifstream configFile;
  configFile.open(configFilePath);
  ofstream outputFile;
  outputFile.open(outputFilePath, std::ios::app);
  outputFile << std::setprecision(15);

  string line;
  if (configFile.is_open())
  {
    if(outputFile.is_open())
    {
      // write config.json to output file
      while( getline (configFile, line) )
      {
        outputFile << line << "\n";
      }

      const double *xs = min->X();
      const double *errors = min->Errors();

      outputFile << "astau: \t = " << xs[0] << "\t +/- \t " << errors[0] << "\n";
      outputFile << "aGG: \t = " << xs[1] << "\t +/- \t " << errors[1] << "\n";
      outputFile << "rho: \t = " << xs[2] << "\t +/- \t " << errors[2] << "\n";
      outputFile << "c8: \t = " << xs[3] << "\t +/- \t " << errors[3] << "\n";

      outputFile << endl;
      outputFile.close();
    }
    configFile.close();
  }
}

int main (int argc, char* argv[]) {
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  cout.precision(17);
  string outputFilePath = "./output/fits.csv";
  string configFilePath = "./configuration.json";
  if (argc > 1) {
    printf("Error no outputfile argument. Try ./build/FESR ./output/folder/fits.dat");
    outputFilePath = argv[1];
  }

  ifstream configFile(configFilePath);
  json jsonConfig;
  configFile >> jsonConfig;
  const Configuration config(jsonConfig);
  const Chisquared chisquared(config);
  const TheoreticalMoments th(config);
  // const DualityViolations dv;
  // cout << "lol" << endl;
  // cout << dv.moment(1., Weight(1), 3.0, 0.0018963380072094527, -2.2, 3.9) << endl;

  // MINUIT
  Minimizer* min = Factory::CreateMinimizer("Minuit2", "Migrad");

  // set tolerances
  // min->SetMaxFunctionCalls(10000000); // for Minuit2
  // min->SetMaxIterations(10000000); // for GSL
  min->SetTolerance(1e-15);
  min->SetStrategy(2);
  min->SetPrintLevel(3); // activate logging

  // function wrapper
  Functor chi2(chisquared, 4);

  min->SetFunction(chi2);

  // set free variables to be minimized
  if (config.astau.isFixed) {
    min->SetFixedVariable(0, "astau", config.astau.value);
  } else {
    min->SetVariable(0, "astau", config.astau.value, config.astau.stepSize);
  }
  if (config.aGGInv.isFixed) {
    min->SetFixedVariable(1, "aGGInv", config.aGGInv.value);
  } else {
    min->SetVariable(1, "aGGInv", config.aGGInv.value, config.aGGInv.stepSize);
  }
  if (config.rhoVpA.isFixed) {
    min->SetFixedVariable(2, "rhoVpA", config.rhoVpA.value);
  } else {
    min->SetVariable(2, "rhoVpA", config.rhoVpA.value, config.rhoVpA.stepSize);
  }
  if (config.c8VpA.isFixed) {
    min->SetFixedVariable(3, "c8VpA", config.c8VpA.value);
  } else {
    min->SetVariable(3, "c8VpA", config.c8VpA.value, config.c8VpA.stepSize);
  }
  // min->SetVariable(4, "vKappa", 3.0, 0.1);
  // min->SetVariableLowerLimit(4, 0);
  // min->SetVariable(5, "vGamma", 1.2, 0.1);
  // min->SetFixedVariable(6, "vAlpha", -2.2);
  // min->SetFixedVariable(7, "vBeta", 3.9);
  // min->SetVariable(8, "aKappa", 3.0, 0.1);
  // min->SetVariableLowerLimit(8, 0);
  // min->SetFixedVariable(9, "aGamma", 1.3);
  // min->SetFixedVariable(10, "aAlpha", 4.7);
  // min->SetFixedVariable(11, "aBeta", 1.8);

  // minimize!
  min->Minimize();
  // std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
  // std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()*1e-6  <<std::endl;
  // const double *xs = min->X();
  // const double *errors = min->Errors();
  // // chisquared.log(xs[0], xs[1], xs[2], xs[3]);
  // // const double chi2AtMin = chisquared(config.s0Set, xs[0], xs[1], xs[2], xs[3]);



  writeOutput("./configuration.json", outputFilePath, min);
  return 0;
}
