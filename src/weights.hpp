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
    case 1: return wD00(x);
    case 2: return wD01(x);
    case 3: return wD02(x);
    case 4: return wD03(x);
    case 5: return wD10(x);
    case 6: return wD11(x);
    case 7: return wD12(x);
    case 8: return wD13(x);
    case 9: return wD20(x);
    case 10: return wD21(x);
    case 11: return wD22(x);
    case 12: return wD23(x);
    }
    return wD00(x);
  }

  dComplex wR(dComplex x) const {
    switch(weightId_) {
    case 1: return wR00(x);
    case 2: return wR01(x);
    case 3: return wR02(x);
    case 4: return wR03(x);
    case 5: return wR10(x);
    case 6: return wR11(x);
    case 7: return wR12(x);
    case 8: return wR13(x);
    case 9: return wR20(x);
    case 10: return wR21(x);
    case 11: return wR22(x);
    case 12: return wR23(x);
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
  static dComplex wD01(dComplex x) {
    return pow(1.0 - x, 3)*(3.0 + 9.0*x + 8.0*pow(x, 2))/10.0;
  }
  static dComplex wR02(dComplex x) {
    return pow(1.0 - x, 2)*(1.0 + 2.0*x)*pow(x, 2);
  }
  static dComplex wD02(dComplex x) {
    return 2.0*pow(1.0 - x, 3)*(1.0 + 3.0*x + 6.0*pow(x, 2) + 5.0*pow(x, 3))/15.0;
  }
  static dComplex wR03(dComplex x) {
    return pow(1.0 - x, 2)*(1.0 + 2.0*x)*pow(x, 3);
  }
  static dComplex wD03(dComplex x) {
    return pow(1.0 - x, 3)*(1.0 + 3.0*x + 6.0*pow(x, 2) + 10.0*pow(x, 3) + 8.0*pow(x, 4))/14.0;
  }
  static dComplex wR10(dComplex x) {
    return pow(1.0 - x, 3)*(1.0 + 2.0*x);
  }
  static dComplex wD10(dComplex x) {
    return pow(1.0 - x, 4)*(7.0 + 8.0*x)/10.0;
  }
  static dComplex wR11(dComplex x) {
    return pow(1.0 - x, 3)*(1.0 + 2.0*x)*x;
  }
  static dComplex wD11(dComplex x) {
    return pow(1.0 - x, 4)*pow(1.0 + 2.0*x, 2)/6.0;
  }
  static dComplex wR12(dComplex x) {
    return pow(1.0 - x, 3)*(1.0 + 2.0*x)*pow(x, 2);
  }
  static dComplex wD12(dComplex x) {
    return pow(1.0 - x, 4)*(13.0 + 52.0*x + 130.0*pow(x, 2) + 120.0*pow(x, 3))/210.0;
  }
  static dComplex wR13(dComplex x) {
    return pow(1.0 - x, 3)*(1.0 + 2.0*x)*pow(x, 3);
  }
  static dComplex wD13(dComplex x) {
    return pow(1.0 - x, 4)*(2.0 + 8.0*x + 20.0*pow(x, 2) + 40.0*pow(x, 3) + 35.0*pow(x, 4))/70.0;
  }
  static dComplex wR20(dComplex x) {
    return pow(1.0 - x, 4)*(1.0 + 2.0*x);
  }
  static dComplex wD20(dComplex x) {
    return 2.0*pow(1.0 - x, 5)*(4.0 + 5.0*x)/15.0;
  }
  static dComplex wR21(dComplex x) {
    return pow(1.0 - x, 4)*(1.0 + 2.0*x)*x;
  }
  static dComplex wD21(dComplex x) {
    return pow(1.0 - x, 5)*(11.0 + 55.0*x + 60.0*pow(x, 2))/105.0;
  }
  static dComplex wR22(dComplex x) {
    return pow(1.0 - x, 4)*(1.0 + 2.0*x)*pow(x, 2);
  }
  static dComplex wD22(dComplex x) {
    return pow(1.0 - x, 5)*(1.0 + 5.0*x + 15.0*pow(x, 2) + 15.0*pow(x, 3))/30.0;
  }
  static dComplex wR23(dComplex x) {
    return pow(1.0 - x, 4)*(1.0 + 2.0*x)*pow(x, 3);
  }
  static dComplex wD23(dComplex x) {
    return pow(1.0 - x, 5)*(17.0 + 85.0*x + 255.0*pow(x, 2) + 595.0*pow(x, 3) + 560.0*pow(x, 4))/1260.0;
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
