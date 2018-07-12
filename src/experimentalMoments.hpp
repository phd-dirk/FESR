#ifndef SRC_EXPERIMENTALMOMENTS_HPP
#define SRC_EXPERIMENTALMOMENTS_HPP

#include "./types.hpp"
#include "./configuration.hpp"
#include "./data.hpp"
#include "./numerics.hpp"

class ExperimentalMoments : public Numerics {
 public:
  // Init. data, vector of all wRatios, exp. Spectral Moments, error Matrix
  // and Covariance matrix
  //
  // The Spectral Moments and Covariance matrix can then bet exported via the
  // public getter functions
  ExperimentalMoments(const string &filename, const Configuration &config) :
    Numerics(), config_(config), data_(Data(filename, config.RVANormalization)),
    inputs_(config_.inputs)
  {
    // moment count
    for(auto const &input: inputs_) {
      vec s0s = input.s0s;
      momCount_ += s0s.size();
    }

    // cache experimental moments
    initExpMoms();
  }

  vec operator ()() const {
    return expMoms;
  }

  void initExpMoms() {
    for(auto const &input: inputs_) {
      vec s0s = input.s0s;
      Weight w = input.weight;
      for(auto const &s0: s0s) {
        expMoms.push_back(expPlusPionMom(s0, w));
      }
    }
  }

  // Selects the closest bin number from s0
  // If s0 is exactly between two bins we select the smaller one
  int closestBinToS0(const double &s0) const {
    int i = data_.dsbins.size()-1;
    if (s0 > data_.sbins.back())
      return i;

    // -1.e-6: if exactly between two bins choose smaller one as closest
    // e.g. s0 = 2.1; sbins[71] = 2.05; sbins[72] = 2.15 => closest bin: 71
    while(abs(s0-data_.sbins[i]-1.e-6) > data_.dsbins[i]/2.)
      i--;

    return i;
  }

  // returns weight ratio
  double wRatio(const double &s0, const Weight &w, const double &sbin,
                const double &dsbin) const {
    double binRight = sbin+dsbin/2.;
    double binLeft = sbin-dsbin/2.;
    return (s0/config_.sTau*
            (
             (w.wD(binLeft/s0) - w.wD(binRight/s0))/
             (w.wTau(binLeft/config_.sTau) - w.wTau(binRight/config_.sTau))
             )
            ).real();
  }

  // return the experimental spectral moment
  double expMom(const double &s0, const Weight &w) const {
    double mom = 0;
    for(int j = 0; j <= closestBinToS0(s0); j++) {
      mom += config_.sTau/s0/config_.be*data_.sfm2s[j]
        *wRatio(s0, w, data_.sbins[j], data_.dsbins[j]);
    }
    return mom;
  }

  double pionPoleMoment(const double &s0, const Weight &w) const {
    double axialMoment = 0;
    double pseudoMoment = 0;
    axialMoment += kPiFac()/s0*w.wR(pow(config_.kPionMinusMass, 2)/s0).real();
    pseudoMoment += axialMoment*(-2.*pow(config_.kPionMinusMass, 2)/(config_.sTau + 2.*pow(config_.kPionMinusMass, 2)));
    return axialMoment + pseudoMoment;
  }

  double expPlusPionMom(const double &s0, const Weight &w) const {
    return expMom(s0, w) + pionPoleMoment(s0, w);
  }

  void log() const {
    // cout << "s0 \t" << s0s[0] << endl;
    // cout << "ExpMom = \t" << getExpPlusPionMoment(0) << endl;
    // cout << "CovMom = \t" << covarianceMatrix << endl;
    // cout << "InvCov = \t" << inverseCovarianceMatrix << endl;
  }

  double kPiFac() const {
    return 24.*pow(M_PI*config_.kVud*config_.kFPi, 2)*config_.kSEW;
  }
  double kDPiFac() const {
    return kPiFac()*sqrt(4.*pow(config_.kDVud/config_.kVud, 2)
                         + pow(config_.kDSEW/config_.kSEW, 2)
                         + 4.*pow(config_.kDFPi/config_.kFPi, 2));
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
    errMat(data_.binCount, data_.binCount) = pow(config_.dBe, 2);
    errMat(data_.binCount+1, data_.binCount+1) = pow(kDPiFac(), 2);
    return errMat;
  }

  // returns the Jacobian Matrix
  mat jacMat() const {
    mat jac(data_.binCount+2, momCount_);
    int xMom = 0;
    for(auto const &input: inputs_) {
      vec s0s = input.s0s;
      Weight w = input.weight;
      for(auto const &s0: s0s) {
        for (int j = 0; j < data_.binCount+2; j++) {
          if (j <= closestBinToS0(s0)) {
            jac(j, xMom) = config_.sTau/s0/config_.be*wRatio(s0, w, data_.sbins[j], data_.dsbins[j]);
          } else {
            jac(j, xMom) = 0.;
          }
          jac(data_.binCount, xMom) = (pionPoleMoment(s0, w)
                                          - expPlusPionMom(s0, w))/config_.be;
          jac(data_.binCount+1, xMom) = pionPoleMoment(s0, w)/kPiFac();
        }
        xMom++;
      }
    }
    return jac;
  }

  // returns the covariance matrix
  mat covMom() const {
    mat jac = jacMat();
    mat err = errMat();
    mat cov(momCount_, momCount_);
    for(int i = 0; i < momCount_; i++) {
      for(int j = 0; j < momCount_; j++) {
        for (int k = 0; k < data_.binCount+2; k++) {
          for (int l = 0; l < data_.binCount+2; l++) {
            cov(i,j) += jac(k, i)*err(k, l)*jac(l, j);
          }
        }
      }
    }
    return cov;
  }

  Configuration config_;
  const Data data_;
  std::vector<Input> inputs_;
  vec expMoms;
  int momCount_ = 0;
}; // END ExperimentalMoments

#endif
