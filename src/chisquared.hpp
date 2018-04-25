#ifndef SRC_CHISQUARED_H
#define SRC_CHISQUARED_H

#include "./experimentalMoments.hpp"
#include "./theoretical_moments.hpp"
#include "./weights.hpp"
#include "json.hpp"
#include <vector>
#include <functional>
#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>

namespace ublas = boost::numeric::ublas;

using json = nlohmann::json;
using std::vector;
using std::function;
using std::cout;
using std::endl;

class Chisquared {
 public:
  Chisquared(const json &config, const Constants &constants) :
    order_(config["parameters"]["order"]), s0s_(config["parameters"]["s0Set"].get<vector<double>>()),
    expMom_(ExperimentalMoments("/Users/knowledge/Developer/PhD/FESR/aleph.json",
                               config["parameters"]["RVANormalization"],
                                s0s_, Weight(config["parameters"]["weight"].get<int>()), constants)),
    thMom_(TheoreticalMoments(order_, s0s_, Weight(config["parameters"]["weight"].get<int>()), config, constants)) {}

  double operator ()(const double *xx) {
    // init fit parameters
    double astau = xx[0];
    double aGGinv = xx[1];
    double rhoD6VpA = xx[2];
    double c8D8VpA = xx[3];

    double chi = 0;

    ublas::matrix<double> covMat = expMom_.getCovarianceMatrix();
    ublas::matrix<double> invCovMat = expMom_.getInverseCovarianceMatrix();

    vector<double> momDiff(s0s_.size());
    for(uint i = 0; i < s0s_.size(); i++) {
      momDiff[i] = expMom_(i) - thMom_(i, astau, aGGinv, rhoD6VpA, c8D8VpA, order_);
    }

    for(uint k = 0; k < s0s_.size(); k++) {
      for(uint l = 0; l < s0s_.size(); l++) {
        chi += momDiff[k] * invCovMat(k, l) * momDiff[l];
      }
    }

    return chi;
  }

  void log(const double &astau, const double &aGGinv, const double &rhoVpa, const double &c8Vpa) {
    cout << "Theoretical Moments:" << endl;
    thMom_.log(astau, aGGinv, rhoVpa, c8Vpa, order_);
    cout << endl;
    cout << "Experimental Moments:" << endl;
    expMom_.log();
    cout << endl;
  }

  const int order_;
  const vector<double> s0s_;
  ExperimentalMoments expMom_;
  TheoreticalMoments thMom_;
};

#endif
