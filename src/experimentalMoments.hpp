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
    Numerics(), config_(config), data_(Data(filename, config.RVANormalization)),
    inputs_(config_.inputs)
  {
    // init jacobian matrix
    // initJacobianMatrix();

    // init covariance matrix
    // initCovarianceMatrix();

    // Remove correlations with R_tau,V+A in Aleph fit
    // for (uint i = 1; i < s0s.size(); i++) {
    //   covarianceMatrix(0, i) = 0.;
    //   covarianceMatrix(i, 0) = 0.;
    // }
    // employ uncertainity of R_VA = 3.4718(72) (HFLAV 2017)
    // covarianceMatrix(0, 0) = pow(0.0072, 2);

    // init inverse covariance matrix
    // invertMatrix(covarianceMatrix, inverseCovarianceMatrix);
  }
  vec operator ()() const {
    vec expMoms;
    for(auto const &input: inputs_) {
      vec s0s = input.s0s;
      Weight w = input.weight;
      for(auto const &s0: s0s) {
        expMoms.push_back(expMom(s0, w) + pionPoleMoment(s0, w));
      }
    }
    return expMoms;
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

  void log() const {
    // cout << "s0 \t" << s0s[0] << endl;
    // cout << "ExpMom = \t" << getExpPlusPionMoment(0) << endl;
    cout << "CovMom = \t" << covarianceMatrix << endl;
    cout << "InvCov = \t" << inverseCovarianceMatrix << endl;
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
  mat errMatrix() {
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
  // void initJacobianMatrix() {
  //   mat<double> jacobi(data.binCount+2, s0s.size());
  //   for (uint i = 0; i < s0s.size(); i++) {
  //     for (int j = 0; j < data.binCount+2; j++) {
  //       if (j <= closestBinToS0(s0s[i])) {
  //         jacobi(j, i) = config_.sTau/s0s[i]/config_.be*weightRatios(i, j);
  //       } else {
  //         jacobi(j, i) = 0.;
  //       }
  //     }
  //     jacobi(data.binCount, i) = (pionPoleMoment(s0s[i]) - getExpPlusPionMoment(i))/config_.be;
  //     jacobi(data.binCount+1, i) = pionPoleMoment(s0s[i])/kPiFac();
  //   }
  //   this->jacobianMatrix = jacobi;
  // }

  // returns the covariance matrix
  // void initCovarianceMatrix() {
  //   for (uint i = 0; i < s0s.size(); i++) {
  //     for (uint j = 0; j < s0s.size(); j++) {
  //       covarianceMatrix(i, j) = 0.;
  //       for (int k = 0; k < data.binCount+2; k++) {
  //         for (int l = 0; l < data.binCount+2; l++) {
  //           covarianceMatrix(i,j) += jacobianMatrix(k, i)*errorMatrix(k, l)*jacobianMatrix(l, j);
  //         }
  //       }
  //     }
  //   }
  // }


  Configuration config_;
  const Data data_;
  std::vector<Input> inputs_;
  vec experimentalMoments;
  mat jacobianMatrix, covarianceMatrix,
    inverseCovarianceMatrix;
}; // END ExperimentalMoments

#endif
