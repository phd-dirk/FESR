#ifndef SRC_CHISQUARED_H
#define SRC_CHISQUARED_H

#include "./experimentalMoments.hpp"
#include "./numerics.hpp"
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

class Chisquared : public Numerics {
 public:
  Chisquared(const int &nc, const int &nf, const int &order, const vector<double> &s0s,
             function<complex<double>(complex<double>)> weight) : Numerics(1e-13, 0),
      s0s(s0s),
      expMom(ExperimentalMoments("/Users/knowledge/Developer/PhD/FESR/aleph.json",
                                 0.99363, s0s, weight, wD00)),
      thMom(TheoreticalMoments(nc, nf, order, s0s, weight)) {
  }

  double operator ()() {
    double chi = 0;

    vector<double> momDiff(s0s.size());
    for(int i = 0; i < s0s.size(); i++) {
      momDiff[i] = expMom.getExpPlusPionMoment(i) - thMom(i);
    }

    ublas::matrix<double> covMat = expMom.getCovarianceMatrix();
    ublas::matrix<double> invCovarianceMatrix(9, 9);

    // Remove correlations with R_tau, V+A in Aleph fit
    for (int i = 1; i < 9; i++) {
      covMat(0, i) = 0.;
      covMat(i, 0) = 0.;
    }

    invertMatrix(covMat, invCovarianceMatrix);

    cout << "covMat : \t" << covMat << endl;
    cout << "invCov : \t" << invCovarianceMatrix << endl;

    for(int k = 0; k < s0s.size(); k++) {
      for(int l = 0; l < s0s.size(); l++) {
        chi += momDiff[k] * invCovarianceMatrix(k, l) * momDiff[l];
      }
    }

    return chi;
  }

 private:
  vector<double> s0s;
  ExperimentalMoments expMom;
  TheoreticalMoments thMom;
};

#endif
