#include <iostream>
#include <fstream>
#include "./src/data.hpp"


int main ()
{
  const sfm2 sfm2Array = readData("./../aleph.json");
  std::cout << sfm2Array[22] << std::endl;
  return 0;
}
