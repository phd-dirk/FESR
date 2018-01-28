#include "./src/data.hpp"
#include "./src/s0_sets.hpp"
#include "./src/state.hpp"
#include "./src/weights.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <complex>

using std::cout;
using std::endl;
using std::complex;

std::string projectRoot = "/Users/knowledge/Developer/PhD/FESR";


int main () {
  cout.precision(17);
  State state(s0Set, 80, projectRoot+"/aleph.json", wR00);
  renormalizeState(state, 0.99363);

  cout << state.weight(complex<double>(1., 2.)) << endl;

  return 0;
}
