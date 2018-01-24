#include "json.hpp"
#include <vector>
#include <string>
#include <fstream>

typedef std::vector<double> sfm2;

using json = nlohmann::json;

sfm2 getSfm2(const json& json)
{
  return json["data"]["sbin"];
}

sfm2 readData(const std::string filename) {
  std::ifstream file(filename);
  json json;
  file >> json;
  return getSfm2(json);
}

