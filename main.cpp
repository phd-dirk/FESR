
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
#include "./src/ope.hpp"
#include "./src/alpha_s.hpp"
#include "./src/numerics.hpp"
#include "./src/duality_violations.hpp"
#include "./src/condensates.hpp"

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
#include <boost/spirit/include/karma.hpp>
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
    const ThMoms thMoms(
      3, 3, { 2.8e-3, 5.0e-3, 97.0e-3 },
      Condensates(
        0.3156,
        { -pow(0.272, 3), -pow(0.272, 3), 0.8*-pow(0.272, 3) },
        { 2.8e-3, 5.0e-3, 97.0e-3 }
      ),
      pow(1.77686, 2),
      0.13957018,
      92.21e-3,
      2.2e-3, 1.3, 0.4,
      0.19e-3, 1.8, 0.21,
      {
        {
          Weight(1),
          {
            3.1572314596000002, 3.0, 2.800, 2.600, 2.400, 2.300,
            2.200, 2.100, 2.000
          }
        }
      },
      { "FO", true, true, true, false, true },
      0.97420,
      1.0198
    );

    std::cout << "thmoms \t" << thMoms(
      pow(1.77686, 2), Weight(1), 0.3179, 0.021, -0.15, 0.24, 5, pow(1.77686, 2),
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.13957018, 92.21e-3,
     2.2e-3, 1.3, 0.4, 0.19e-3, 1.8, 0.21
    ) << std::endl;
    return 0;



    // MINUIT
    Minimizer* min = Factory::CreateMinimizer("Minuit2", "Migrad");

    // set tolerances
    // min->SetMaxFunctionCalls(10000000); // for Minuit2
    // min->SetMaxIterations(10000000); // for GSL
    // min->SetTolerance(1e-8);
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
    if (config.deltaV.isFixed) {
      min->SetFixedVariable(4, "deltaV", config.deltaV.value);
    } else {
      min->SetVariable(4, "deltaV", config.deltaV.value, config.deltaV.stepSize);
    }
    if (config.gammaV.isFixed) {
      min->SetFixedVariable(5, "gammaV", config.gammaV.value);
    } else {
      min->SetVariable(5, "gammaV", config.gammaV.value, config.gammaV.stepSize);
    }
    if (config.alphaV.isFixed) {
      min->SetFixedVariable(6, "alphaV", config.alphaV.value);
    } else {
      min->SetVariable(6, "alphaV", config.alphaV.value, config.alphaV.stepSize);
    }
    if (config.betaV.isFixed) {
      min->SetFixedVariable(7, "betaV", config.betaV.value);
    } else {
      min->SetVariable(7, "betaV", config.betaV.value, config.betaV.stepSize);
    }
    if (config.deltaA.isFixed) {
      min->SetFixedVariable(8, "deltaA", config.deltaA.value);
    } else {
      min->SetVariable(8, "deltaA", config.deltaA.value, config.deltaA.stepSize);
    }
    if (config.gammaA.isFixed) {
      min->SetFixedVariable(9, "gammaA", config.gammaA.value);
    } else {
      min->SetVariable(9, "gammaA", config.gammaA.value, config.gammaA.stepSize);
    }
    if (config.alphaA.isFixed) {
      min->SetFixedVariable(10, "alphaA", config.alphaA.value);
    } else {
      min->SetVariable(10, "alphaA", config.alphaA.value, config.alphaA.stepSize);
    }
    if (config.betaA.isFixed) {
      min->SetFixedVariable(11, "betaA", config.betaA.value);
    } else {
      min->SetVariable(11, "betaA", config.betaA.value, config.betaA.stepSize);
    }

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
