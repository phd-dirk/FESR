#ifndef SRC_WEIGHTS_H
#define SRC_WEIGHTS_H

#include <cmath>
#include <complex>
#include <functional>

using std::pow;
using std::complex;
using std::function;

typedef function<complex<double>(complex<double>)> complexFunction;
typedef complex<double> dComplex;

class Weight {
 public:
  Weight(uint i) : weightId_(i) {}

  dComplex wD(dComplex x) const {
    switch(weightId_) {
    case 1:
      return wD00(x);
    }
    return wD00(x);
  }

  dComplex wR(dComplex x) const {
    switch(weightId_) {
    case 1:
      return wR00(x);
    }
    return wR00(x);
  }

  static dComplex wTau(dComplex x) {
    return wD00(x);
  }
  static dComplex wR00(dComplex x) {
    return pow((1. - x), 2)*(1. + 2.*x);
  }
  static dComplex wD00(dComplex x) {
    return pow((1. - x), 3)*(1. + x);
  }

 private:
  uint weightId_;
};

// inline complex<double> wR00(complex<double> x) {
//   return pow((1. - x), 2)*(1. + 2.*x);
// }

// inline complex<double> wD00(complex<double> x) {
//   return pow((1. - x), 3)*(1. + x);
// }

#endif
