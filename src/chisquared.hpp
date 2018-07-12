#ifndef SRC_CHISQUARED_H
#define SRC_CHISQUARED_H

#include "./configuration.hpp"
#include "./experimentalMoments.hpp"
#include "./theoretical_moments.hpp"
#include "./weights.hpp"
#include "json.hpp"
#include <vector>
#include <functional>
#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>


// test invCovMat
#include <boost/numeric/ublas/io.hpp>
using ublas::matrix;
using ublas::prod;

using std::ifstream;

namespace ublas = boost::numeric::ublas;

matrix<double> readMatrixFromFile(const int &size, const string filePath) {
  matrix<double> m(size, size);
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

using json = nlohmann::json;
using std::vector;
using std::function;
using std::cout;
using std::endl;

class Chisquared {
 public:
  Chisquared(Configuration config) :
    s0Set_(config.s0Set), order_(config.order),
    expMom_(ExperimentalMoments("/Users/knowledge/Developer/PhD/FESR/aleph.json", config)),
    thMom_(TheoreticalMoments(config)) {}

  double operator ()(const double *xx) const {
    // init fit parameters
    double astau = xx[0];
    double aGGinv = xx[1];
    double rhoD6VpA = xx[2];
    double c8D8VpA = xx[3];

    return chi2(s0Set_, astau, aGGinv, rhoD6VpA, c8D8VpA);
  }

  double operator ()(const vec s0Set, const double &astau, const double &aGGinv, const double &rho, const double &c8) const {
    return chi2(s0Set, astau, aGGinv, rho, c8);
  }

  double operator ()(const json &config) const {
    double xx[4];
    xx[0] = config["variables"]["astau"]["value"];
    xx[1] = config["variables"]["aGGInv"]["value"];
    xx[2] = config["variables"]["rhoVpA"]["value"];
    xx[3] = config["variables"]["c8VpA"]["value"];
    return operator()(xx);
  }

  double chi2(const vec s0Set, const double &astau, const double &aGGinv,
              const double &rho, const double &c8) const {
    double chi = 0;

    ublas::matrix<double> covMat = expMom_.covarianceMatrix;
    ublas::matrix<double> invCovMat = expMom_.inverseCovarianceMatrix; //readMatrixFromFile(9, "./data/invCovMat.dat");
    vec momDiff(s0Set.size());
    for(uint i = 0; i < s0Set.size(); i++) {
      const double s0 = s0Set[i];
      // momDiff[i] = expMom_(i) - thMom_(s0, astau, aGGinv, rho, c8, order_);
    }

    for(uint k = 0; k < s0Set.size(); k++) {
      for(uint l = 0; l < s0Set.size(); l++) {
        chi += momDiff[k] * invCovMat(k, l) * momDiff[l];
      }
    }

    return chi;
  }

  void log(const double &astau, const double &aGGinv, const double &rhoVpa, const double &c8Vpa) const {
    cout << "Theoretical Moments:" << endl;
    thMom_.log(astau, aGGinv, rhoVpa, c8Vpa, order_);
    cout << endl;
    cout << "Experimental Moments:" << endl;
    expMom_.log();
    cout << endl;
  }

  const vec s0Set_;
  const int order_;
  const ExperimentalMoments expMom_;
  const TheoreticalMoments thMom_;
 private:
};

#endif
