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
  string configFilePath = "./configuration.json";
  if (argc == 2) {
    configFilePath = argv[1];
  }

  try {
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
    min->SetTolerance(1e-10);
    min->SetStrategy(2);
    min->SetPrintLevel(3); // activate logging

    // function wrapper
    Functor chi2(chisquared, 12);

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
    min->SetVariable(4, "deltaV", 3.9, 0.1);
    min->SetVariable(5, "gammaV", 0.29, 0.1);
    min->SetVariable(8, "alphaV", -1.02, 0.1);
    min->SetVariable(9, "betaV", 3.64, 0.1);
    min->SetVariable(10, "deltaA", 1.82, 0.1);
    min->SetVariable(11, "gammaA", 1.46, 0.1);
    min->SetVariable(11, "alphaA", -2.5, 0.1);
    min->SetVariable(11, "betaA", 2.91, 0.1);

    // minimize!
    min->Minimize();

    std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()*1e-6  <<std::endl;

    writeOutput(configFilePath, min, config);
    return 0;
  }
  catch (const std::exception& e) {
    cout << e.what() << endl;
    // writeOutput(e.what(), outputFilePath);
    return 1;
  }
}
