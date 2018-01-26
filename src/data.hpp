#include "json.hpp"
#include <vector>
#include <string>
#include <fstream>
#include <boost/numeric/ublas/matrix.hpp>

using json = nlohmann::json;
using namespace boost::numeric::ublas;

const int n = 80;

// Is the central data storage
// Example:
//    Data data(80);
//    cout << data.sfm2[0] << endl;
//    cout << data.corerr(0, 0) << endl;
struct Data {
public:
  Data(int size) : sbin(80) {}
  std::vector<double> sbin;
  std::vector<double> dsbin;
  std::vector<double> sfm2;
  std::vector<double> derr;
  matrix<double> corerr;
};

// Returns the corerr matrix from the json object
//
// This is necassary, because the Aleph data (given in Fortran) cannot be
// exported as a clean json file.
// e.g.
// {
//   "data": {
//     "corerr01": [...],
//     "corerr02": [...]
//   }
// }
matrix<double> getCorErr(const json& json) {
  matrix<double> m(80, 80);
  for(int i = 0; i < n; i++) {
    std::string corerrRowName = "corerr";
    if (i < 9) {
      corerrRowName += "0" + std::to_string(i+1);
    } else {
      corerrRowName += std::to_string(i+1);
    }
    std::vector<double> corerrRow = json["data"][corerrRowName];
    for(int j = 0; j < n; j++) {
      m(i, j) = corerrRow[j];
    }
  }
  return m;
}

// Initalizes the Data struct
// Takes the data.json file and reads it into memory.
// The json file needs to be of the following form:
//
// {
//   "data": {
//     "sbin": [...],
//     "dsbin": [...],
//     "sfm2": [...],
//     "derr": [...],
//     "corerr01": [...],
//     "corerr02": [...]
//     ...
//   }
// }
Data readData(const int size, const std::string filename) {
  std::ifstream file(filename);
  json json;
  file >> json;
  Data data(size);
  data.sbin = json["data"]["sbin"].get<std::vector<double>>();
  // getCorErr(json);
  // return getSfm2(json);
  data.dsbin = json["data"]["dsbin"].get<std::vector<double>>();
  data.sfm2 = json["data"]["sfm2"].get<std::vector<double>>();
  data.derr = json["data"]["derr"].get<std::vector<double>>();
  data.corerr= getCorErr(json);
  return data;
}

