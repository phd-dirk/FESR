#include <iostream>
#include <fstream>
#include <string>
#include "./src/data.hpp"
#include <boost/numeric/ublas/matrix.hpp>

std::string projectRoot = "/Users/knowledge/Developer/PhD/FESR";

int main ()
{
  std::vector<double> vec(80);
  vec[0] = 1.2;
  const Data data = readData(80, projectRoot+"/aleph.json");
  std::cout << data.corerr(0,0) << std::endl;
  return 0;
}
