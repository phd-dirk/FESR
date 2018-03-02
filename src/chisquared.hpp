#ifndef SRC_CHISQUARED_H
#define SRC_CHISQUARED_H

#include "./experimentalMoments.hpp"
#include "./theoretical_moments.hpp"
#include "./weights.hpp"
#include <vector>
#include <functional>

using std::vector;
using std::function;

class Chisquared {
 public:
  Chisquared(const int &nc, const int &nf, const int &order, const vector<double> &s0s,
             function<complex<double>(complex<double>)> weight) :
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

    for(int k = 0; k < s0s.size(); k++) {
      for(int l = 0; l < s0s.size(); l++) {
        chi += momDiff[k] * expMom.getCovarianceMatrix(k, l) * momDiff[l];
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
