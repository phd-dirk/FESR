#ifndef SRC_EXPERIMENTALMOMENTS_HPP
#define SRC_EXPERIMENTALMOMENTS_HPP

#include "./data.hpp"
#include "./numerics.hpp"
#include "./constants.hpp"
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
                      const vector<double> &s0s, const Weight &weight, Constants constants) :
    Numerics(constants), const_(constants), data(Data(filename, normalizationFactor)),
    s0s(s0s), weight_(weight) {

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

  void log() {
    cout << "ExpMom = \t" << getExpPlusPionMoment(0) << endl;
    cout << "InvCov = \t" << getInverseCovarianceMatrix() << endl;
  }

  vector<double> getExpPlusPionPoleMoments() {
    vector<double> expPlusPionMoments(s0s.size());
    for (uint i = 0; i < s0s.size(); i++) {
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
    for (uint i = 1; i < s0s.size(); i++) {
      covMat(0, i) = 0.;
      covMat(i, 0) = 0.;
    }
    // employ uncertainity of R_VA = 3.4718(72) (HFLAV 2017)
    covMat(0, 0) = pow(0.0072, 2);
    invertMatrix(covMat, invCovMat);
    return invCovMat;
  }

  // Selects the closest bin number from s0
  // If s0 is exactly between two bins we select the smaller one
  int closestBinToS0(const double &s0) {
    int i = data.dsbins.size()-1;
    if (s0 > data.sbins.back())
      return i;

    // -1.e-6: if exactly between two bins choose smaller one as closest
    // e.g. s0 = 2.1; sbins[71] = 2.05; sbins[72] = 2.15 => closest bin: 71
    while(abs(s0-data.sbins[i]-1.e-6) > data.dsbins[i]/2.)
      i--;

    return i;
  }

  double kPiFac() {
    return 24.*pow(const_.kPi*const_.kVud*const_.kFPi, 2)*const_.kSEW;
  }
  double kDPiFac() {
    return kPiFac()*sqrt(4.*pow(const_.kDVud/const_.kVud, 2)
                         + pow(const_.kDSEW/const_.kSEW, 2)
                         + 4.*pow(const_.kDFPi/const_.kFPi, 2));
  }

 private:
  Constants const_;
  Data data;
  vector<double> s0s, experimentalMoments;
  Weight weight_;
  matrix<double> weightRatios, errorMatrix, jacobianMatrix, covarianceMatrix;

  // returns weightRatios
  void setWeightRatios() {
    matrix<double> wRatios(s0s.size(), data.binCount);

    for (uint i = 0; i < s0s.size(); i++) {
      for (int j = 0; j < data.binCount; j++) {
        double s0UpperLimit = data.sbins[j]+data.dsbins[j]/2.;
        double s0LowerLimit = data.sbins[j]-data.dsbins[j]/2.;
        wRatios(i, j) = (s0s[i]/const_.kSTau*(
            (weight_.wD(s0LowerLimit/s0s[i]) - weight_.wD(s0UpperLimit/s0s[i]))/
            (weight_.wTau(s0LowerLimit/const_.kSTau) - weight_.wTau(s0UpperLimit/const_.kSTau)))).real();
      }
    }

    this->weightRatios = wRatios;
  }


  // return the experimental spectral moment
  void setExperimentalMoments() {
    vector<double> moments(s0s.size());
    for (uint i = 0; i < s0s.size(); i++) {
      for(int j = 0; j <= closestBinToS0(s0s[i]); j++) {
        moments[i] += const_.kSTau/s0s[i]/const_.kBe*data.sfm2s[j]*weightRatios(i, j);
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
    errMat(data.binCount, data.binCount) = pow(const_.kDBe, 2);
    errMat(data.binCount+1, data.binCount+1) = pow(kDPiFac(), 2);
    this->errorMatrix = errMat;
  }

  // returns the Jacobian Matrix
  void setJacobianMatrix() {
    matrix<double> jacobi(data.binCount+2, data.binCount+2);
    for (uint i = 0; i < s0s.size(); i++) {
      for (int j = 0; j < data.binCount+2; j++) {
        if (j <= closestBinToS0(s0s[i])) {
          jacobi(j, i) = const_.kSTau/s0s[i]/const_.kBe*weightRatios(i, j);
        } else {
          jacobi(j, i) = 0.;
        }
      }
      jacobi(data.binCount, i) = (pionPoleMoment(s0s[i]) - getExpPlusPionMoment(i))/const_.kBe;
      jacobi(data.binCount+1, i) = pionPoleMoment(s0s[i])/kPiFac();
    }
    this->jacobianMatrix = jacobi;
  }

  // returns the covariance matrix
  void setCovarianceMatrix() {
    matrix<double> covMat(s0s.size(), s0s.size());
    for (uint i = 0; i < s0s.size(); i++) {
      for (uint j = 0; j < s0s.size(); j++) {
        covMat(i, j) = 0.;
        for (int k = 0; k < data.binCount+2; k++) {
          for (int l = 0; l < data.binCount+2; l++) {
            covMat(i,j) = covMat(i, j) + jacobianMatrix(k, i)*errorMatrix(k, l)*jacobianMatrix(l, j);
          }
        }
      }
    }
    this->covarianceMatrix = covMat;
  }


  double pionPoleMoment(const double &s0) {
    double axialMoment = 0;
    double pseudoMoment = 0;
    axialMoment += kPiFac()/s0*weight_.wR(pow(const_.kPionMinusMass, 2)/s0).real();
    pseudoMoment += axialMoment*(-2.*pow(const_.kPionMinusMass, 2)/(const_.kSTau + 2.*pow(const_.kPionMinusMass, 2)));
    return axialMoment + pseudoMoment;
  }

}; // END ExperimentalMoments

#endif
