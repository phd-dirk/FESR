#ifndef SRC_DATA_H
#define SRC_DATA_H

#include "./constants.hpp"
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

using json = nlohmann::json;
using boost::numeric::ublas::matrix;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::abs;
using std::complex;
using std::function;
using std::transform;

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
inline matrix<double> getCorErr(const json &json, const int &binCount) {
  matrix<double> m(80, 80);
  for(int i = 0; i < binCount; i++) {
    std::string corerrRowName = "corerr";
    if (i < 9) {
      corerrRowName += "0" + std::to_string(i+1);
    } else {
      corerrRowName += std::to_string(i+1);
    }
    vector<double> corerrRow = json["data"][corerrRowName];
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
    this->sbins = json["data"]["sbin"].get<vector<double>>();
    this->binCount = sbins.size();
    this->dsbins = json["data"]["dsbin"].get<vector<double>>();
    this->sfm2s = json["data"]["sfm2"].get<vector<double>>();
    this->derrs = json["data"]["derr"].get<vector<double>>();
    this->corerrs = getCorErr(json, binCount);

    normalizeData(normalizationFactor);
  }

  // normalize sfm2s and derr (multiply data vector by 'factor' scalar)
  void normalizeData(const double &factor) {
    transform(sfm2s.begin(), sfm2s.end(), sfm2s.begin(),
              std::bind1st(std::multiplies<double>(), factor));
    transform(derrs.begin(), derrs.end(), derrs.begin(),
              std::bind1st(std::multiplies<double>(), factor));
  }

  int binCount;
  vector<double> sbins;
  vector<double> dsbins;
  vector<double> sfm2s;
  vector<double> derrs;
  matrix<double> corerrs;
};

class ExperimentalMoments : public Constants {
 public:
  // Init. data, vector of all wRatios, exp. Spectral Moments, error Matrix
  // and Covariance matrix
  //
  // The Spectral Moments and Covariance matrix can then bet exported via the
  // public getter functions
  ExperimentalMoments(const string &filename, const double &normalizationFactor,
                      const vector<double> &s0s,
                      function<complex<double>(complex<double>)> weight,
                      function<complex<double>(complex<double>)> wTau) :
      data(Data(filename, normalizationFactor)), s0s(s0s), weight(weight),
      wTau(wTau) {

    // init weightRatios
    setWeightRatios(); // (s0s.size() x data.binCount) e.g. (9 x 80)

    // init exp. moments
    setExperimentalMoments();

    // // init error matrix
    this->errorMatrix = setErrorMatrix();

    // // init jacobian matrix
    this->jacobianMatrix = setJacobianMatrix();

    // // init covariance matrix
    this->covarianceMatrix = setCovarianceMatrix();
  }

  vector<double> getExpPlusPionPoleMoments() {
    vector<double> expPlusPionMoments(s0s.size());
    for (int i = 0; i < s0s.size(); i++) {
      expPlusPionMoments[i] = experimentalMoments[i]
          + pionPoleMoment(s0s[i]);
    }
    return expPlusPionMoments;
  }
  double getExpPlusPionMoment(int i) {
    return experimentalMoments[i] + pionPoleMoment(s0s[i]);
  }

  vector<double> getExperimentalMoments() {
    return this->experimentalMoments;
  }

  matrix<double> getWeightRatios() {
    return this->weightRatios;
  }
  double getErrorMatrix(int i, int j) {
    return this->errorMatrix(i, j);
  }
  double getJacobianMatrix(int i, int j) {
    return this->jacobianMatrix(i, j);
  }

  // Selects the closest bin number from s0
  // If s0 is exactly between two bins we select the smaller one
  int closestBinToS0(const double &s0) {
    int i = data.dsbins.size()-1;
    // -1.e-6: if exactly between two bins choose smaller one as closest
    // e.g. s0 = 2.1; sbins[71] = 2.05; sbins[72] = 2.15 => closest bin: 71
    while(abs(s0-data.sbins[i]-1.e-6) > data.dsbins[i]/2.) {
      i--;
    }
    return i;
  }

 private:
  Data data;
  vector<double> s0s, experimentalMoments;
  function<complex<double>(complex<double>)> weight, wTau;
  matrix<double> weightRatios, errorMatrix, jacobianMatrix, covarianceMatrix;


  // returns weightRatios
  void setWeightRatios() {
    matrix<double> wRatios(s0s.size(), data.binCount);

    for (int i = 0; i < s0s.size(); i++) {
      for (int j = 0; j < data.binCount; j++) {
        double s0UpperLimit = data.sbins[j]+data.dsbins[j]/2.;
        double s0LowerLimit = data.sbins[j]-data.dsbins[j]/2.;
        wRatios(i, j) = (s0s[i]/kSTauMass*(
            (weight(s0LowerLimit/s0s[i]) - weight(s0UpperLimit/s0s[i]))/
            (wTau(s0LowerLimit/kSTauMass) - wTau(s0UpperLimit/kSTauMass)))).real();
      }
    }

    this->weightRatios = wRatios;
  }


  // return the experimental spectral moment
  void setExperimentalMoments() {
    vector<double> moments(s0s.size());
    for (int i = 0; i <= s0s.size(); i++) {
      for(int j = 0; j <= closestBinToS0(s0s[i]); j++) {
        moments[i] += kSTauMass/s0s[i]/kBe*data.sfm2s[j]*weightRatios(i, j);
      }
    }
    this->experimentalMoments = moments;
  }

  // returns the error matrix
  matrix<double> setErrorMatrix() {
    matrix<double> errorMatrix(data.binCount+2, data.binCount+2);
    for (int i = 0; i < data.binCount; i++) {
      for (int j = 0; j < data.binCount; j++) {
        errorMatrix(i, j) = data.corerrs(i, j)*data.derrs[i]*data.derrs[j]/1.e2;
      }
    }
    errorMatrix(data.binCount, data.binCount) = pow(kDBe, 2);
    errorMatrix(data.binCount+1, data.binCount+1) = pow(kDRTauVex, 2);
    return errorMatrix;
  }

  // returns the Jacobian Matrix
  matrix<double> setJacobianMatrix() {
    matrix<double> jacobi(data.binCount+2, data.binCount+2);
    for (int i = 0; i < s0s.size(); i++) {
      for (int j = 0; j < data.binCount; j++) {
        jacobi(j, i) = kSTauMass/s0s[i]/kBe*weightRatios(i, j);
      }
      jacobi(data.binCount, i) = (pionPoleMoment(s0s[i]) - getExpPlusPionMoment(i))/kBe;
      jacobi(data.binCount+1, i) = pionPoleMoment(s0s[i])/kPiFac();
    }
    cout << pionPoleMoment(s0s[0]) << endl;
    return jacobi;
  }

  // returns the covariance matrix
  matrix<double> setCovarianceMatrix() {
    matrix<double> covMat(s0s.size(), s0s.size());
    for (int i = 0; i < s0s.size(); i++) {
      for (int j = 0; j < s0s.size(); j++) {
        for (int k = 0; k < data.binCount; k++) {
          for (int l = 0; l < data.binCount; l++) {
            covMat(i,j) += jacobianMatrix(k, i)*errorMatrix(k, l)*jacobianMatrix(l, j);
          }
        }
      }
    }
    return covMat;
  }

  double kPiFac() {
    return 24.*pow(Constants::kPi*Constants::kVud*Constants::kFPi, 2)*Constants::kSEW;
  }

  double pionPoleMoment(const double &s0) {
    double axialMoment = 0;
    double pseudoMoment = 0;
    axialMoment += kPiFac()/s0*wR00(pow(Constants::kPionMinusMass, 2)/s0).real();
    pseudoMoment += axialMoment*(-2.*pow(kPionMinusMass, 2)/(kSTauMass + 2.*pow(kPionMinusMass, 2)));
    return axialMoment + pseudoMoment;
  }

}; // END ExperimentalMoments


// Return an vector with every entry multiplied by a factor
//
// Is used in state.hpp to mutate the state
inline vector<double> renormalize(double renormalizationFactor, vector<double> vec) {
  vector<double> renormalizedVector;
  std::cout << "max_size: " << renormalizedVector.max_size() << std::endl;
  for(auto const& value: vec) {
    renormalizedVector.push_back(value*renormalizationFactor);
  }
  return renormalizedVector;
}




#endif
