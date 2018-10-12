#ifndef SRC_EXPERIMENTALMOMENTS_HPP
#define SRC_EXPERIMENTALMOMENTS_HPP

#include "./types.hpp"
#include "./configuration.hpp"
#include "./data.hpp"

class ExperimentalMoments {
 public:
  // Init. data, vector of all wRatios, exp. Spectral Moments, error Matrix
  // and Covariance matrix
  //
  // The Spectral Moments and Covariance matrix can then bet exported via the
  // public getter functions
  ExperimentalMoments(const string &filename, const Configuration &config);

  vec operator ()() const;

  mat getCovMat () const {
    return covMat;
  }

  void initExpMoms();

  // return the experimental spectral moment
  double expMom(const double &s0, const Weight &w) const;
  double pionPoleMoment(const double &s0, const Weight &w) const;
  double expPlusPionMom(const double &s0, const Weight &w) const {
    return expMom(s0, w) + pionPoleMoment(s0, w);
  }

  // Selects the closest bin number from s0
  // If s0 is exactly between two bins we select the smaller one
  int closestBinToS0(const double &s0) const;

  // returns weight ratio
  double wRatio(
    const double &s0, const Weight &w, const double &sbin,
    const double &dsbin
  ) const;


  void log() const {
    // cout << "s0 \t" << s0s[0] << endl;
    // cout << "ExpMom = \t" << getExpPlusPionMoment(0) << endl;
    // cout << "CovMom = \t" << covarianceMatrix << endl;
    // cout << "InvCov = \t" << inverseCovarianceMatrix << endl;
  }

  double kPiFac() const {
    return 24.*pow(M_PI*config_.vud_*config_.kFPi, 2)*config_.kSEW;
  }
  double kDPiFac() const {
    return kPiFac()*sqrt(
      4.*pow(config_.dVud_/config_.vud_, 2)
      + pow(config_.kDSEW/config_.kSEW, 2)
      + 4.*pow(config_.kDFPi/config_.kFPi, 2)
    );
  }

  // returns the error matrix
  mat errMat() const {
    mat errMat(data_.binCount+2, data_.binCount+2);
    for (int i = 0; i < data_.binCount+2; i++) {
      for (int j = 0; j < data_.binCount+2; j++) {
        if (i < data_.binCount && j < data_.binCount) {
          errMat(i, j) = data_.corerrs(i, j)*data_.derrs[i]*data_.derrs[j]/1.e2;
        } else {
          errMat(i, j) = 0.;
        }
      }
    }
    errMat(data_.binCount, data_.binCount) = pow(config_.dBe_, 2);
    errMat(data_.binCount+1, data_.binCount+1) = pow(kDPiFac(), 2);
    return errMat;
  }

  // returns the Jacobian Matrix
  mat jacMat() const {
    mat jac(data_.binCount+2, config_.momCount_);
    int xMom = 0;
    for(auto const &input: config_.inputs_) {
      vec s0s = input.s0s;
      Weight w = input.weight;
      for(auto const &s0: s0s) {
        for (int j = 0; j < data_.binCount+2; j++) {
          if (j <= closestBinToS0(s0)) {
            jac(j, xMom) = config_.sTau_/s0/config_.be_*wRatio(s0, w, data_.sbins[j], data_.dsbins[j]);
          } else {
            jac(j, xMom) = 0.;
          }
          jac(data_.binCount, xMom) = (
            pionPoleMoment(s0, w) - expPlusPionMom(s0, w)
          )/config_.be_;
          jac(data_.binCount+1, xMom) = pionPoleMoment(s0, w)/kPiFac();
        }
        xMom++;
      }
    }
    return jac;
  }

  // returns the covariance matrix
  void initCovMat() {
    mat jac = jacMat();
    mat err = errMat();
    for(uint i = 0; i < config_.momCount_; i++) {
      for(uint j = 0; j < config_.momCount_; j++) {
        covMat(i, j) = 0.0;
        for (int k = 0; k < data_.binCount+2; k++) {
          for (int l = 0; l < data_.binCount+2; l++) {
            covMat(i, j) += jac(k, i)*err(k, l)*jac(l, j);
          }
        }
      }
    }
  }

  Configuration config_;
  const Data data_;
  vec expMoms;
  mat covMat;
}; // END ExperimentalMoments

#endif
