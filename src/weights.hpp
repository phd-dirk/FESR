#ifndef SRC_WEIGHTS_H
#define SRC_WEIGHTS_H

#include <cmath>
#include <complex>

using std::pow;
using std::complex;

inline complex<double> wR00(complex<double> x) {
  return pow((1. - x), 2)*(1. + 2.*x);
}

inline complex<double> wD00(complex<double> x) {
  return pow((1. - x), 3)*(1. + x);
}

#endif
