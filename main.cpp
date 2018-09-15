#define BOOST_UBLAS_NDEBUG 1 // numerical issues with ublas matrix inverse

#include <chrono>
#include <iostream>
#include <vector>
#include <complex>
#include <stdexcept>
#include "./src/types.hpp"

using std::vector;
using std::cout;
using std::endl;
using std::complex;
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

#include "./src/utils.hpp"

int main (int argc, char* argv[]) {
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  cout.precision(17);
  string outputFilePath = "./output/fit.txt";
  string configFilePath = "./configuration.json";
  if (argc > 1) {
    printf("Error no outputfile argument. Try ./build/FESR ./output/folder/fits.dat");
    outputFilePath = argv[1];
  }

  const Configuration config(configFilePath);
  const Chisquared chisquared(config);

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

  writeOutput("./configuration.json", outputFilePath, min, config);
  return 0;
}
