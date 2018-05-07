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
  Weight(int i) : weightId_(i) {}

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
  // Polunomial weights
  static dComplex wR00(dComplex x) {
    return pow((1. - x), 2)*(1. + 2.*x);
  }
  static dComplex wD00(dComplex x) {
    return pow((1. - x), 3)*(1. + x);
  }
  static dComplex wR01(dComplex x) {
    return pow(1.0 - x, 2)*(1.0 + 2.0*x)*x;
  }
  static dComplex wR02(dComplex x) {
    return pow(1.0 - x, 2)*(1.0 + 2.0*x)*pow(x, 2);
  }
  static dComplex wR03(dComplex x) {
    return pow(1.0 - x, 2)*(1.0 + 2.0*x)*pow(x, 3);
  }
  static dComplex wR10(dComplex x) {
    return pow(1.0 - x, 3)*(1.0 + 2.0*x);
  }
  static dComplex wR11(dComplex x) {
    return pow(1.0 - x, 3)*(1.0 + 2.0*x)*x;
  }
  static dComplex wR12(dComplex x) {
    return pow(1.0 - x, 3)*(1.0 + 2.0*x)*pow(x, 2);
  }
  static dComplex wR13(dComplex x) {
    return pow(1.0 - x, 3)*(1.0 + 2.0*x)*pow(x, 3);
  }
  static dComplex wR20(dComplex x) {
    return pow(1.0 - x, 4)*(1.0 + 2.0*x);
  }
  static dComplex wR21(dComplex x) {
    return pow(1.0 - x, 4)*(1.0 + 2.0*x)*x;
  }
  static dComplex wR22(dComplex x) {
    return pow(1.0 - x, 4)*(1.0 + 2.0*x)*pow(x, 2);
  }
  static dComplex wR23(dComplex x) {
    return pow(1.0 - x, 4)*(1.0 + 2.0*x)*pow(x, 3);
  }

 private:
  int weightId_;
};

// inline complex<double> wR00(complex<double> x) {
//   return pow((1. - x), 2)*(1. + 2.*x);
// }

// inline complex<double> wD00(complex<double> x) {
//   return pow((1. - x), 3)*(1. + x);
// }

#endif
