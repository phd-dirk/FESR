#include "./constants.hpp"
#include "json.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <complex>
#include <functional>

using json = nlohmann::json;
using boost::numeric::ublas::matrix;
using std::vector;
using std::string;
using std::abs;
using std::complex;
using std::function;


// number of data points
const int n = 80;

// Selects the closest bin number from s0
// If s0 is exactly between two bins we select the smaller one
int closestBinToS0(double s0, vector<double> sbins, vector<double> dsbins) {
  int i = dsbins.size()-1;
  // -1.e-6: if exactly between two bins choose smaller one as closest
  // e.g. s0 = 2.1; sbins[71] = 2.05; sbins[72] = 2.15 => closest bin: 71
  while(abs(s0-sbins[i]-1.e-6) > dsbins[i]/2.) {
    i--;
  }
  return i;
}

//
complex<double> expSpectralMoment(
                                  double s0, vector<double> sfm2s,
                                  vector<double> sbins, vector<double> dsbins,
                                  function<complex<double>(complex<double>)> weight,
                                  function<complex<double>(complex<double>)> wTau,
                                  double sTau, double be
                                  ) {
  complex<double> mom(0., 0.);
  for(int i = 0; i <= closestBinToS0(s0, sbins, dsbins); i++) {
    double s0UpperLimit = sbins[i]+dsbins[i]/2.;
    double s0LowerLimit = sbins[i]-dsbins[i]/2.;
    complex<double> wRatio = s0/sTau*(weight(s0LowerLimit/s0) - weight(s0UpperLimit/s0)/
                              wTau(s0LowerLimit/sTau) - wTau(s0UpperLimit/sTau));
    mom += sTau/s0/be*sfm2s[i]*wRatio;
  }
  return mom;
}

// Return an vector with every entry multiplied by a factor
//
// Is used in state.hpp to mutate the state
vector<double> renormalize(double renormalizationFactor, vector<double> vec) {
  vector<double> renormalizedVector;
  std::cout << "max_size: " << renormalizedVector.max_size() << std::endl;
  for(auto const& value: vec) {
    renormalizedVector.push_back(value*renormalizationFactor);
  }
  return renormalizedVector;
}

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
    this->sbins = json["data"]["sbin"].get<vector<double>>();
    this->dsbins = json["data"]["dsbin"].get<vector<double>>();
    this->sfm2s = json["data"]["sfm2"].get<vector<double>>();
    this->derrs = json["data"]["derr"].get<vector<double>>();
    this->corerrs = getCorErr(json);
  }
  vector<double> sbins;
  vector<double> dsbins;
  vector<double> sfm2s;
  vector<double> derrs;
  matrix<double> corerrs;
};
