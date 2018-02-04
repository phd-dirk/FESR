#include "./src/data.hpp"
#include "./src/s0_sets.hpp"
#include "./src/state.hpp"
#include "./src/weights.hpp"
#include "./src/theoretical_moments.hpp"
#include "./src/numerics.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <gsl/gsl_integration.h>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <string>
#include <complex>


#include "./lib/cRunDec/CRunDec.h"
using std::cout;
using std::endl;
using std::complex;
using std::log;
using std::sqrt;
using std::printf;

std::string projectRoot = "/Users/knowledge/Developer/PhD/FESR";

inline double test(double x) {
  return x*x;
}


int main () {
  cout.precision(17);
  ExperimentalMoments expMom(projectRoot+"/aleph.json", 0.99363, s0Set, wD00, wD00);
  State state(s0Set, projectRoot+"/aleph.json", wR00);
  renormalizeState(state, 0.99363);
  cout << state.weight(complex<double>(1., 2.)) << endl;

  AdlerFunction adler(3, 3, 4);

  cout << adler.D0(-1., 2.) << endl;


  // gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

  // double result, error;
  // double expected = -4.0;
  // double alpha = 1.0;

  // gsl_function F;
  // F.function = &f;
  // F.params = &alpha;

  // gsl_integration_qag(&F, 0, 1, 1e-13, 0, 1000, 6, w, &result, &error);

  // printf("result = % .18f\n", result);

  Numerics numerics(1e-7, 1e-7);
  cout << numerics.integrate(Numerics::test) << endl;

  return 0;
}
