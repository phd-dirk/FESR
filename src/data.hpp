#include "json.hpp"
#include <vector>
#include <string>
#include <fstream>
#include <boost/numeric/ublas/matrix.hpp>

using json = nlohmann::json;
using boost::numeric::ublas::matrix;
using std::vector;
using std::string;

const int n = 80;

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
    vector<double> corerrRow = json["data"][corerrRowName];
    for(int j = 0; j < n; j++) {
      m(i, j) = corerrRow[j];
    }
  }
  return m;
}

// Is the central data storage
// Example:
//    Data data(80);
//    cout << data.sfm2[0] << endl;
//    cout << data.corerr(0, 0) << endl;
struct Data {
public:
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
  Data(const int size, const string filename) {
    std::ifstream file(filename);
    json json;
    file >> json;
    this->sbin = json["data"]["sbin"].get<vector<double>>();
    this->dsbin = json["data"]["dsbin"].get<vector<double>>();
    this->sfm2 = json["data"]["sfm2"].get<vector<double>>();
    this->derr = json["data"]["derr"].get<vector<double>>();
    this->corerr= getCorErr(json);
  }
  vector<double> sbin;
  vector<double> dsbin;
  vector<double> sfm2;
  vector<double> derr;
  matrix<double> corerr;
};
