#include "./src/data.hpp"
#include "./src/s0_sets.hpp"
#include "./src/state.hpp"
#include "./src/weights.hpp"
#include "./src/theoretical_moments.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <complex>

#include "./lib/cRunDec/CRunDec.h"
using std::cout;
using std::endl;
using std::complex;

std::string projectRoot = "/Users/knowledge/Developer/PhD/FESR";

int main () {
  cout.precision(17);
  ExperimentalMoments expMom(projectRoot+"/aleph.json", 0.99363, s0Set, wD00, wD00);
  State state(s0Set, projectRoot+"/aleph.json", wR00);
  renormalizeState(state, 0.99363);
  cout << state.weight(complex<double>(1., 2.)) << endl;

  AdlerFunction adler(3, 3, 4);
  cout << adler.c_[5][5] << endl;
  cout << adler.D0(-1., 2.) << endl;

  return 0;
}
