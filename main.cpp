#include <iostream>
#include <vector>

using std::vector;
using std::cout;
using std::endl;

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

double RosenBrock(const double *xx) {
  const Double_t x = xx[0];
  const Double_t y = xx[1];
  const Double_t tmp1 = y-x*x;
  const Double_t tmp2 = 1-x;
  return 100*tmp1*tmp1*tmp2*tmp2;
}

int main () {
  cout.precision(17);

  int nc = 3;
  int nf = 3;
  int order = 5;

  vector<double> s0s = s0Set;

  Chisquared chisquared(nc, nf, order, s0s, wD00);
  cout << chisquared() << endl;


  // MINUIT
  // Minimizer* min = Factory::CreateMinimizer("Minuit2", "Migrad");

  // // set tolerances
  // min->SetMaxFunctionCalls(10000000); // for Minuit2
  // min->SetMaxIterations(10000000); // for GSL
  // min->SetTolerance(1e-17);
  // // min->SetPrintLevel(1); // activate logging

  // // function wrapper
  // Functor f(&RosenBrock, 2);
  // double step[2] = {0.01, 0.01};
  // // starting point
  // double variable[2] = {-1., 1.2};

  // min->SetFunction(f);

  // // set free variables to be minimized
  // min->SetVariable(0, "x", variable[0], step[0]);
  // min->SetVariable(1, "y", variable[1], step[1]);

  // // minimize!
  // min->Minimize();

  // const double *xs = min->X();
  // cout << "Minimum f(" << xs[0] << "," << xs[1] << "): "
  //      << min->MinValue() << endl;

  // cout << "min " << RosenBrock(xs) << endl;

  return 0;
}
