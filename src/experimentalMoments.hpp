#ifndef SRC_EXPERIMENTALMOMENTS_HPP
#define SRC_EXPERIMENTALMOMENTS_HPP

#include "./types.hpp"
#include "./configuration.hpp"
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
  ExperimentalMoments(const string &filename, const Configuration &config) :
    Numerics(), config_(config), data(Data(filename, config.RVANormalization)),
    s0s(config.s0Set), weight_(config.weight), covarianceMatrix(s0s.size(), s0s.size()),
    inverseCovarianceMatrix(s0s.size(), s0s.size()) {

    // init weightRatios
    initWeightRatios(); // (s0s.size() x data.binCount) e.g. (9 x 80)

    // init exp. moments
    initExperimentalMoments();

    // init error matrix
    initErrorMatrix();

    // init jacobian matrix
    initJacobianMatrix();

    // init covariance matrix
    initCovarianceMatrix();

    // Remove correlations with R_tau,V+A in Aleph fit
    for (uint i = 1; i < s0s.size(); i++) {
      covarianceMatrix(0, i) = 0.;
      covarianceMatrix(i, 0) = 0.;
    }
    // employ uncertainity of R_VA = 3.4718(72) (HFLAV 2017)
    covarianceMatrix(0, 0) = pow(0.0072, 2);

    // init inverse covariance matrix
    invertMatrix(covarianceMatrix, inverseCovarianceMatrix);
  }
  double operator ()(int i) const {
    return getExpPlusPionMoment(i);
  }

  void log() const {
    cout << "s0 \t" << s0s[0] << endl;
    cout << "ExpMom = \t" << getExpPlusPionMoment(0) << endl;
    cout << "CovMom = \t" << covarianceMatrix << endl;
    cout << "InvCov = \t" << inverseCovarianceMatrix << endl;
  }

  vector<double> getExpPlusPionPoleMoments() const {
    vector<double> expPlusPionMoments(s0s.size());
    for (uint i = 0; i < s0s.size(); i++) {
      expPlusPionMoments[i] = experimentalMoments[i]
          + pionPoleMoment(s0s[i]);
    }
    return expPlusPionMoments;
  }
  double getExpPlusPionMoment(int i) const {
    return experimentalMoments[i] + pionPoleMoment(s0s[i]);
  }

  vector<double> getExperimentalMoments() {
    return this->experimentalMoments;
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

  double kPiFac() const {
    return 24.*pow(M_PI*config_.kVud*config_.kFPi, 2)*config_.kSEW;
  }
  double kDPiFac() const {
    return kPiFac()*sqrt(4.*pow(config_.kDVud/config_.kVud, 2)
                         + pow(config_.kDSEW/config_.kSEW, 2)
                         + 4.*pow(config_.kDFPi/config_.kFPi, 2));
  }


  // returns weightRatios
  void initWeightRatios() {
    matrix<double> wRatios(s0s.size(), data.binCount);

    for (uint i = 0; i < s0s.size(); i++) {
      for (int j = 0; j < data.binCount; j++) {
        double s0UpperLimit = data.sbins[j]+data.dsbins[j]/2.;
        double s0LowerLimit = data.sbins[j]-data.dsbins[j]/2.;
        wRatios(i, j) = (s0s[i]/config_.sTau*(
            (weight_.wD(s0LowerLimit/s0s[i]) - weight_.wD(s0UpperLimit/s0s[i]))/
            (weight_.wTau(s0LowerLimit/config_.sTau) - weight_.wTau(s0UpperLimit/config_.sTau)))).real();
      }
    }

    this->weightRatios = wRatios;
  }


  // return the experimental spectral moment
  void initExperimentalMoments() {
    vector<double> moments(s0s.size());
    for (uint i = 0; i < s0s.size(); i++) {
      for(int j = 0; j <= closestBinToS0(s0s[i]); j++) {
        moments[i] += config_.sTau/s0s[i]/config_.be*data.sfm2s[j]*weightRatios(i, j);
      }
    }
    this->experimentalMoments = moments;
  }

  // returns the error matrix
  void initErrorMatrix() {
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
    errMat(data.binCount, data.binCount) = pow(config_.dBe, 2);
    errMat(data.binCount+1, data.binCount+1) = pow(kDPiFac(), 2);
    this->errorMatrix = errMat;
  }

  // returns the Jacobian Matrix
  void initJacobianMatrix() {
    matrix<double> jacobi(data.binCount+2, s0s.size());
    for (uint i = 0; i < s0s.size(); i++) {
      for (int j = 0; j < data.binCount+2; j++) {
        if (j <= closestBinToS0(s0s[i])) {
          jacobi(j, i) = config_.sTau/s0s[i]/config_.be*weightRatios(i, j);
        } else {
          jacobi(j, i) = 0.;
        }
      }
      jacobi(data.binCount, i) = (pionPoleMoment(s0s[i]) - getExpPlusPionMoment(i))/config_.be;
      jacobi(data.binCount+1, i) = pionPoleMoment(s0s[i])/kPiFac();
    }
    this->jacobianMatrix = jacobi;
  }

  // returns the covariance matrix
  void initCovarianceMatrix() {
    for (uint i = 0; i < s0s.size(); i++) {
      for (uint j = 0; j < s0s.size(); j++) {
        covarianceMatrix(i, j) = 0.;
        for (int k = 0; k < data.binCount+2; k++) {
          for (int l = 0; l < data.binCount+2; l++) {
            covarianceMatrix(i,j) += jacobianMatrix(k, i)*errorMatrix(k, l)*jacobianMatrix(l, j);
          }
        }
      }
    }
  }

  double pionPoleMoment(const double &s0) const {
    double axialMoment = 0;
    double pseudoMoment = 0;
    axialMoment += kPiFac()/s0*weight_.wR(pow(config_.kPionMinusMass, 2)/s0).real();
    pseudoMoment += axialMoment*(-2.*pow(config_.kPionMinusMass, 2)/(config_.sTau + 2.*pow(config_.kPionMinusMass, 2)));
    return axialMoment + pseudoMoment;
  }

  Configuration config_;
  Data data;
  vector<double> s0s, experimentalMoments;
  Weight weight_;
  matrix<double> weightRatios, errorMatrix, jacobianMatrix, covarianceMatrix, inverseCovarianceMatrix;
}; // END ExperimentalMoments

#endif
