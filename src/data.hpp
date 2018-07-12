#ifndef SRC_DATA_H
#define SRC_DATA_H

#include "./types.hpp"
#include "./weights.hpp"
#include "json.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <complex>
#include <functional>
#include <algorithm>

using std::string;
using std::abs;
using std::complex;
using std::function;
using std::transform;
using std::sqrt;

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
inline mat getCorErr(const json &json, const int &binCount) {
  mat m(80, 80);
  for(int i = 0; i < binCount; i++) {
    std::string corerrRowName = "corerr";
    if (i < 9) {
      corerrRowName += "0" + std::to_string(i+1);
    } else {
      corerrRowName += std::to_string(i+1);
    }
    vec corerrRow = json["data"][corerrRowName];
    for(int j = 0; j < binCount; j++) {
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
  Data(const string &filename, const double &normalizationFactor) {
    std::ifstream file(filename);
    json json;
    file >> json;
    this->sbins = json["data"]["sbin"].get<vec>();
    this->binCount = sbins.size();
    this->dsbins = json["data"]["dsbin"].get<vec>();
    this->sfm2s = json["data"]["sfm2"].get<vec>();
    this->derrs = json["data"]["derr"].get<vec>();
    this->corerrs = getCorErr(json, binCount);

    normalizeData(normalizationFactor);
  }

  // normalize sfm2s and derr (multiply data vector by 'factor' scalar)
  void normalizeData(const double &factor) {
    for(auto &value: sfm2s ) {
      value *= factor;
    }
    for(auto &value: derrs ) {
      value *= factor;
    }
  }

  int binCount;
  vec sbins;
  vec dsbins;
  vec sfm2s;
  vec derrs;
  mat corerrs;
};



// Return an vector with every entry multiplied by a factor
//
// Is used in state.hpp to mutate the state
inline vec renormalize(double renormalizationFactor, vec v) {
  vec renormalizedVector;
  std::cout << "max_size: " << renormalizedVector.max_size() << std::endl;
  for(auto const& value: v) {
    renormalizedVector.push_back(value*renormalizationFactor);
  }
  return renormalizedVector;
}

#endif
