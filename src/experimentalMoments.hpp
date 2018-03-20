#ifndef SRC_EXPERIMENTALMOMENTS_HPP
#define SRC_EXPERIMENTALMOMENTS_HPP

#include "./data.hpp"
#include "./numerics.hpp"
#include <boost/numeric/ublas/matrix.hpp>

namespace ublas = boost::numeric::ublas;

class ExperimentalMoments : public Numerics {
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
    Numerics(1e-13, 0), data(Data(filename, normalizationFactor)), s0s(s0s),
    weight(weight), wTau(wTau) {

    // init weightRatios
    setWeightRatios(); // (s0s.size() x data.binCount) e.g. (9 x 80)

    // init exp. moments
    setExperimentalMoments();

    // init error matrix

    setErrorMatrix();

    // init jacobian matrix
    setJacobianMatrix();

    // init covariance matrix
    setCovarianceMatrix();
  }

  double operator ()(int i) {
    return getExpPlusPionMoment(i);
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
  matrix<double> getCovarianceMatrix() {
    return this->covarianceMatrix;
  }
  double getCovarianceMatrix(int i, int j) {
    return this->covarianceMatrix(i, j);
  }
  ublas::matrix<double> getInverseCovarianceMatrix() {
    // Remove correlations with R_tau, V+A in Aleph fit
    ublas::matrix<double> covMat = this->covarianceMatrix;
    ublas::matrix<double> invCovMat(s0s.size(), s0s.size());
    for (int i = 1; i < 9; i++) {
      covMat(0, i) = 0.;
      covMat(i, 0) = 0.;
    }
    invertMatrix(covMat, invCovMat);
    return invCovMat;
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

  double kPiFac() {
    return 24.*pow(Constants::kPi*Constants::kVud*Constants::kFPi, 2)*Constants::kSEW;
  }
  double kDPiFac() {
    return kPiFac()*sqrt(4.*pow(kDVud/kVud, 2) + pow(kDSEW/kSEW, 2) + 4.*pow(kDFPi/kFPi, 2));
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
        wRatios(i, j) = (s0s[i]/kSTau*(
            (weight(s0LowerLimit/s0s[i]) - weight(s0UpperLimit/s0s[i]))/
            (wTau(s0LowerLimit/kSTau) - wTau(s0UpperLimit/kSTau)))).real();
      }
    }

    this->weightRatios = wRatios;
  }


  // return the experimental spectral moment
  void setExperimentalMoments() {
    vector<double> moments(s0s.size());
    for (int i = 0; i < s0s.size(); i++) {
      for(int j = 0; j <= closestBinToS0(s0s[i]); j++) {
        moments[i] += kSTau/s0s[i]/kBe*data.sfm2s[j]*weightRatios(i, j);
      }
    }
    this->experimentalMoments = moments;
  }

  // returns the error matrix
  void setErrorMatrix() {
    matrix<double> errMat(data.binCount+2, data.binCount+2);
    for (int i = 0; i < data.binCount+2; i++) {
      for (int j = 0; j < data.binCount+2; j++) {
        if (i < data.binCount && j < data.binCount) {
          errMat(i, j) = data.corerrs(i, j)*data.derrs[i]*data.derrs[j]/1.e2;
        } else {
          errMat(i, j) = 0.;
        }
      }
    }
    errMat(data.binCount, data.binCount) = pow(kDBe, 2);
    errMat(data.binCount+1, data.binCount+1) = pow(kDPiFac(), 2);
    this->errorMatrix = errMat;
  }

  // returns the Jacobian Matrix
  void setJacobianMatrix() {
    matrix<double> jacobi(data.binCount+2, data.binCount+2);
    for (int i = 0; i < s0s.size(); i++) {
      for (int j = 0; j < data.binCount+2; j++) {
        if (j <= closestBinToS0(s0s[i])) {
          jacobi(j, i) = kSTau/s0s[i]/kBe*weightRatios(i, j);
        } else {
          jacobi(j, i) = 0.;
        }
      }
      jacobi(data.binCount, i) = (pionPoleMoment(s0s[i]) - getExpPlusPionMoment(i))/kBe;
      jacobi(data.binCount+1, i) = pionPoleMoment(s0s[i])/kPiFac();
    }
    this->jacobianMatrix = jacobi;
  }

  // returns the covariance matrix
  void setCovarianceMatrix() {
    matrix<double> covMat(s0s.size(), s0s.size());
    for (int i = 0; i < s0s.size(); i++) {
      for (int j = 0; j < s0s.size(); j++) {
        covMat(i, j) = 0.;
        for (int k = 0; k < data.binCount+2; k++) {
          for (int l = 0; l < data.binCount+2; l++) {
            covMat(i,j) = covMat(i, j) + jacobianMatrix(k, i)*errorMatrix(k, l)*jacobianMatrix(l, j);
            // if (i == 0 && j == 0 && k == 0) {
            //   cout << "l=" << l << "\t" << covMat(i,j) << endl;
            // }
          }
        }
      }
    }
    // cout << std::setprecision(17);
    // cout << jacobianMatrix(0, 0)*errorMatrix(0, 0)*jacobianMatrix(0, 0) << endl;
    // cout << jacobianMatrix(0, 0)*errorMatrix(0, 1)*jacobianMatrix(1, 0) << endl;
    this->covarianceMatrix = covMat;
  }


  double pionPoleMoment(const double &s0) {
    double axialMoment = 0;
    double pseudoMoment = 0;
    axialMoment += kPiFac()/s0*wR00(pow(Constants::kPionMinusMass, 2)/s0).real();
    pseudoMoment += axialMoment*(-2.*pow(kPionMinusMass, 2)/(kSTau + 2.*pow(kPionMinusMass, 2)));
    return axialMoment + pseudoMoment;
  }

}; // END ExperimentalMoments

#endif
