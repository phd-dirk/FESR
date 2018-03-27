#include <iostream>
#include <vector>
#include <complex>
#include <string>

using std::vector;
using std::cout;
using std::endl;
using std::complex;
using std::string;

#include "./src/constants.hpp"
// Theoretical Moments
#include "./src/theoretical_moments.hpp"

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

using ROOT::Math::Minimizer;
using ROOT::Math::Factory;
using ROOT::Math::Functor;

void log(string name, complex<double> x) {
  cout << name << ": \t" << x << endl;
}

int main () {
  cout.precision(17);

  int nc = 3;
  int nf = 3;
  int order = 5;

  vector<double> s0s = s0Set;

  Constants constants(nc, nf);
  Chisquared chisquared(order, s0s, wD00, constants);

  // MINUIT
  Minimizer* min = Factory::CreateMinimizer("Minuit2", "Migrad");

  // set tolerances
  min->SetMaxFunctionCalls(10000000); // for Minuit2
  min->SetMaxIterations(10000000); // for GSL
  min->SetTolerance(1e-14);
  min->SetPrintLevel(1); // activate logging

  // function wrapper
  Functor chi2(chisquared, 4);

  min->SetFunction(chi2);

  // set free variables to be minimized
  // min->SetVariable(0, "astau", 0.31927, 0.2e-2);
  min->SetFixedVariable(0, "astau", 0.31927);
  min->SetVariable(1, "aGGInv", 0.21e-1, 1.e-2);
  // min->SetFixedVariable(1, "aGGInv", 0.21e-1);
  // min->SetVariable(2, "rhoVpA",  -0.1894, 0.1);
  min->SetFixedVariable(2, "rhoVpA",  -0.1894);
  // min->SetVariable(3, "c8VpA",  0.16315, 0.3);
  min->SetFixedVariable(3, "c8VpA",  0.16315);

  // minimize!
  min->Minimize();

  const double *xs = min->X();

  cout << "Minimum chi2("
       << xs[0] << ", " << xs[1] << ", "
       << xs[2] << ", " << xs[3] << "): "
       << min->MinValue() << endl;


  min->PrintResults();
  cout << "f() :\t" << chi2(xs) << endl;

  // TheoreticalMoments thMom(order, s0s, wD00, constants);
  // log("delVpA0FO", thMom.cIntVpAD0FO(s0s[0], wD00, 0.31927, 5));
  // log("deltaP", thMom.deltaP(constants.kSTau, wR00));
  double xs2[4] = { 0.31921, 0.21e-1, -0.18939792247957590, 0.16314594513667133 };
  cout << "fotro() \t" << chi2(xs2) << endl;
  return 0;
}
