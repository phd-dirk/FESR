#ifndef SRC_CHISQUARED_H
#define SRC_CHISQUARED_H

#include "./experimentalMoments.hpp"
#include "./theoretical_moments.hpp"
#include "./weights.hpp"
#include <vector>
#include <functional>

using std::vector;
using std::function;

class Chisquared : public ExperimentalMoments {
 public:
  Chisquared(const vector<double> &s0s, function<complex<double>(complex<double>)> weight) :
      ExperimentalMoments("/Users/knowledge/Developer/PhD/FESR/aleph.json",
                          0.99363, s0s, weight, wD00) {
  }

  double operator ()() {

  }
};

#endif
