#include <iostream>
#include <vector>
#include <complex>
#include <string>

using std::vector;
using std::cout;
using std::endl;
using std::complex;
using std::string;

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

  TheoreticalMoments thMom(nc, nf, order, s0s, wD00);
  log("D68", thMom.D68CInt(s0s[0], wD00, -0.1894, 0.16315));
  Chisquared chisquared(nc, nf, order, s0s, wD00);

  // // // MINUIT
  // Minimizer* min = Factory::CreateMinimizer("Minuit2", "Migrad");

  // // set tolerances
  // min->SetMaxFunctionCalls(10000000); // for Minuit2
  // min->SetMaxIterations(10000000); // for GSL
  // min->SetTolerance(1e-17);
  // // min->SetPrintLevel(1); // activate logging

  // // function wrapper
  // Functor chi2(chisquared, 2);

  // double step[2] = { 0.2e-2, 1e-2 };

  // // starting point
  // double variable[2] = { 0.31927, 2.1e-2 };

  // min->SetFunction(chi2);

  // // set free variables to be minimized
  // min->SetVariable(0, "astau", variable[0], step[0]);
  // min->SetVariable(1, "aGGInv", variable[1], step[1]);

  // // minimize!
  // min->Minimize();

  // const double *xs = min->X();

  // cout << "Minimum chi2(" << xs[0] << "," << xs[1] << "): "
  //      << min->MinValue() << endl;

  return 0;
}
