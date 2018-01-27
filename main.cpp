#include <iostream>
#include <fstream>
#include <string>
#include "./src/data.hpp"
#include "./src/s0_sets.hpp"
#include "./src/state.hpp"
#include <boost/numeric/ublas/matrix.hpp>

using std::cout;
using std::endl;

std::string projectRoot = "/Users/knowledge/Developer/PhD/FESR";

int main ()
{
  State state(s0Set, 80, projectRoot+"/aleph.json");
  cout << state.data.sfm2[0] << endl;
  renormalizeState(state, 0.99363);
  cout << state.data.sfm2[0] << endl;

  std::vector<double> vec(80);
  vec[0] = 1.2;
  const Data data(80, projectRoot+"/aleph.json");
  std::cout << data.corerr(1,1) << std::endl;
  return 0;
}
