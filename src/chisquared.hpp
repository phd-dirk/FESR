#ifndef SRC_CHISQUARED_H
#define SRC_CHISQUARED_H

#include "./experimentalMoments.hpp"
#include "./theoretical_moments.hpp"
#include "./weights.hpp"
#include <vector>
#include <functional>
#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>

namespace ublas = boost::numeric::ublas;

using std::vector;
using std::function;
using std::cout;
using std::endl;

class Chisquared {
 public:
  Chisquared(const int &order, const vector<double> &s0s,
             function<complex<double>(complex<double>)> weight, Constants constants) :
      s0s(s0s),
      expMom(ExperimentalMoments("/Users/knowledge/Developer/PhD/FESR/aleph.json",
                                 //0.99363,
                                 1,
                                 s0s, weight, wD00, constants)),
      thMom(TheoreticalMoments(order, s0s, weight, constants)) {
    order_ = order;
  }


  double operator ()(const double *xx) {
    // init fit parameters
    double astau = xx[0];
    double aGGinv = xx[1];
    double rhoD6VpA = xx[2];
    double c8D8VpA = xx[3];

    double chi = 0;

    ublas::matrix<double> covMat = expMom.getCovarianceMatrix();
    ublas::matrix<double> invCovMat = expMom.getInverseCovarianceMatrix();

    vector<double> momDiff(s0s.size());
    for(uint i = 0; i < s0s.size(); i++) {
      momDiff[i] = expMom(i) - thMom(i, astau, aGGinv, rhoD6VpA, c8D8VpA, order_);
    }

    for(uint k = 0; k < s0s.size(); k++) {
      for(uint l = 0; l < s0s.size(); l++) {
        chi += momDiff[k] * invCovMat(k, l) * momDiff[l];
      }
    }

    return chi;
  }

 private:
  int order_;
  vector<double> s0s;
  ExperimentalMoments expMom;
  TheoreticalMoments thMom;
};

#endif
