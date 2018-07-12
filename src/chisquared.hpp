#ifndef SRC_CHISQUARED_H
#define SRC_CHISQUARED_H

#include "./types.hpp"
#include "./configuration.hpp"
#include "./experimentalMoments.hpp"
#include "./theoretical_moments.hpp"
#include "./numerics.hpp"
#include "./weights.hpp"
#include "json.hpp"
#include <functional>
#include <iostream>


// test invCovMat
#include <boost/numeric/ublas/io.hpp>
using ublas::prod;

using std::ifstream;

namespace ublas = boost::numeric::ublas;

mat readMatrixFromFile(const int &size, const string filePath) {
  mat m(size, size);
  ifstream file;
  file.open(filePath);
  if(!file) {
    std::cerr << "Unable to open file expMom.dat";
    exit(1);
  }

  double row = 0;
  double col = 0;
  double x;
  while (file >> x) {
    m(row, col) = x;
    if (col == size-1) {
      col = 0;
      row++;
    } else {
      col++;
    }
  }
  file.close();
  return m;
}

using std::function;
using std::cout;
using std::endl;

class Chisquared: Numerics {
 public:
  Chisquared(Configuration config) :
    inputs_(config.inputs), order_(config.order), momCount_(config.momCount),
    invCovMat_(momCount_, momCount_),
    expMom_(ExperimentalMoments("/Users/knowledge/Developer/PhD/FESR/aleph.json", config)),
    thMom_(TheoreticalMoments(config)) {
    // cache inverse covariance matrix
    initInvCovMat();
  }

  double operator ()( const double *xx) const {
    // init fit parameters
    double astau = xx[0];
    double aGGinv = xx[1];
    double rhoD6VpA = xx[2];
    double c8D8VpA = xx[3];

    return chi2(astau, aGGinv, rhoD6VpA, c8D8VpA);
  }

  double operator ()(const double &astau, const double &aGGinv, const double &rho, const double &c8) const {
    return chi2(astau, aGGinv, rho, c8);
  }

  double operator ()(const json &config) const {
    double xx[4];
    xx[0] = config["variables"]["astau"]["value"];
    xx[1] = config["variables"]["aGGInv"]["value"];
    xx[2] = config["variables"]["rhoVpA"]["value"];
    xx[3] = config["variables"]["c8VpA"]["value"];
    return operator()(xx);
  }

  double chi2(const double &astau, const double &aGGinv,
              const double &rho, const double &c8) const {
    double chi = 0;

    vec momDiff(momCount_);
    for(uint i = 0; i < momCount_; i++) {
      momDiff[i] = expMom_()[i] - thMom_(astau, aGGinv, rho, c8, order_)[i];
    }

    for(uint k = 0; k < momCount_; k++) {
      for(uint l = 0; l < momCount_; l++) {
        chi += momDiff[k] * invCovMat_(k, l) * momDiff[l];
      }
    }
    return chi;
  }

  void initInvCovMat() {
    mat covMat = expMom_.getCovMat();

    // Remove correlations with R_tau,V+A in Aleph fit
    for (uint i = 1; i < momCount_; i++) {
      covMat(0, i) = 0.;
      covMat(i, 0) = 0.;
    }
    // employ uncertainity of R_VA = 3.4718(72) (HFLAV 2017)
    covMat(0, 0) = pow(0.0072, 2);
    Numerics::invertMatrix(covMat, invCovMat_);
  }

  void log(const double &astau, const double &aGGinv, const double &rhoVpa, const double &c8Vpa) const {
    cout << "Theoretical Moments:" << endl;
    thMom_.log(astau, aGGinv, rhoVpa, c8Vpa, order_);
    cout << endl;
    cout << "Experimental Moments:" << endl;
    expMom_.log();
    cout << endl;
  }

  const std::vector<Input> inputs_;
  const int order_;
  const uint momCount_;
  mat invCovMat_;
  const ExperimentalMoments expMom_;
  const TheoreticalMoments thMom_;
 private:
};

#endif
