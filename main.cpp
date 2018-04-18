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

void log(string name, complex<double> x) {
  cout << name << ": \t" << x << endl;
}

int main () {
  cout.precision(17);

  std::ifstream configFile("./configuration.json");
  json config;
  configFile >> config;

  const int nc = config["parameters"]["nc"];
  const int nf = config["parameters"]["nf"];
  const int order = config["parameters"]["order"];

  const vector<double> s0s = config["parameters"]["s0Set"];

  const Constants constants(config);
  const uint weightId = config["parameters"]["weight"];
  const Weight weight(weightId);
  Chisquared chisquared(order, s0s, weight, config, constants);

  chisquared.log(0.32307175541329564);


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
