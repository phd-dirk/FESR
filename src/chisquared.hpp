#ifndef SRC_CHISQUARED_H
#define SRC_CHISQUARED_H

#include "./types.hpp"
#include "./configuration.hpp"
#include "./experimentalMoments.hpp"
#include "./theoretical_moments.hpp"
#include "./numerics.hpp"
#include "./weights.hpp"
#include <nlohmann/json.hpp>
#include <functional>
#include <iostream>
#include <thread>
#include <future>
#include <chrono>


// test invCovMat
#include <boost/numeric/ublas/io.hpp>
using ublas::prod;

using std::ifstream;

namespace ublas = boost::numeric::ublas;

struct Test {
  void func(int x) {
    cout << x << endl;
  }
};

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
    config_(config), inputs_(config.inputs), order_(config.order),
    momCount_(config.momCount), invCovMat_(momCount_, momCount_),
    expMom_(ExperimentalMoments("/Users/knowledge/Developer/PhD/FESR/aleph.json", config)) {
    // cache inverse covariance matrix
    initInvCovMat();
  }

  double operator ()( const double *xx) const {
    // init fit parameters
    double astau = xx[0];
    double aGGinv = xx[1];
    double rhoD6VpA = xx[2];
    double c8D8VpA = xx[3];
    double deV = xx[4];
    double gaV = xx[5];
    double alV = xx[6];
    double beV = xx[7];
    double deA = xx[8];
    double gaA = xx[9];
    double alA = xx[10];
    double beA = xx[11];

    return chi2(astau, aGGinv, rhoD6VpA, c8D8VpA,
                deV, gaV, alV, beV,
                deA, gaA, alA, alA);
  }

  double operator ()(cDbl &astau, cDbl &aGGinv, cDbl &rho, cDbl &c8,
                     cDbl &vKappa, cDbl &vGamma, cDbl &vAlpha, cDbl &vBeta,
                     cDbl &aKappa, cDbl &aGamma, cDbl &aAlpha, cDbl &aBeta) const
  {
    return chi2(astau, aGGinv, rho, c8,
                vKappa, vGamma, vAlpha, vBeta,
                aKappa, aGamma, aAlpha, aBeta);
  }

  double chi2(cDbl &astau, cDbl &aGGinv, cDbl &rho, cDbl &c8,
              cDbl &deV, cDbl &gaV, cDbl &alV, cDbl &beV,
              cDbl &deA, cDbl &gaA, cDbl &alA, cDbl &beA) const
  {
    double chi = 0;

    vec momDiff(momCount_);
    vec thMoms = calcThMoms(astau, aGGinv, rho, c8, order_, deV, gaV, alV,
                            beV, deA, gaA, alA, beA);
    for(uint i = 0; i < momCount_; i++) {
      // parallelized
      momDiff[i] = expMom_()[i] - thMoms[i];
      // without parallization
      // TheoreticalMoments th(config_);
      // momDiff[i] = expMom_()[i] - th.thMom(astau, aGGinv, rho, c8, order_, deV, gaV, alV,
      //                                      beV, deA, gaA, alA, beA)[i];
    }

    for(uint k = 0; k < momCount_; k++) {
      for(uint l = 0; l < momCount_; l++) {
        chi += momDiff[k] * invCovMat_(k, l) * momDiff[l];
      }
    }
    return chi;
  }

  vec calcThMoms(cDbl &astau, cDbl &aGGinv, cDbl &rhoVpA, cDbl &c8VpA, cDbl &order,
                 cDbl &deV, cDbl &gaV, cDbl &alV, cDbl &beV,
                 cDbl &deA, cDbl &gaA, cDbl &alA, cDbl &beA) const
  {
    vec thMoms(config_.momCount);
    std::vector<std::future<double>> ftrs(config_.momCount);
    TheoreticalMoments th(config_);
    int i = 0;
    for(auto const& input: inputs_) {
      vec s0s = input.s0s;
      Weight w = input.weight;
      for(auto const& s0: s0s) {
        ftrs[i] = std::async(&TheoreticalMoments::thMom, &th, s0, w, astau, aGGinv,
                             rhoVpA, c8VpA, order,
                             deV, gaV, alV, beV,
                             deA, gaA, alA, beA);
        i++;
      }
    }

    for(uint i=0; i<config_.momCount; i++) {
      thMoms[i] = ftrs[i].get();
    }

    return thMoms;
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
    // thMom_.log(astau, aGGinv, rhoVpa, c8Vpa, order_);
    cout << endl;
    cout << "Experimental Moments:" << endl;
    expMom_.log();
    cout << endl;
  }

  const Configuration config_;
  const std::vector<Input> inputs_;
  const int order_;
  const uint momCount_;
  mat invCovMat_;
  const ExperimentalMoments expMom_;
 private:
};

#endif
