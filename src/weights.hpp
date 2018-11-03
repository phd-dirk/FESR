#ifndef SRC_WEIGHTS_H
#define SRC_WEIGHTS_H

#include "./types.hpp"
#include <cmath>
#include <complex>
#include <functional>
#include <stdexcept>
#include <string>

struct WeightPolynomial {
  int x0;
  int x1;
  int x2;
  int x3;
};

class Weight {
 public:
  Weight(int i) : weightId_(i) {}

  int getId() const { return weightId_; }

  cmplx wD(cmplx x) const {
    switch(weightId_) {
      case 0: return wDUnit(x);
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
      case 13: return wD2(x);
      default:
        throw std::invalid_argument("weight wD(" + std::to_string(weightId_) + ") is not defined");
        return wD00(x);
    }
  }

  cmplx wR(cmplx x) const {
    switch(weightId_) {
      case 0: return wRUnit();
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
      case 13: return wR2(x);
      default:
        throw std::invalid_argument("weight wR(" + std::to_string(weightId_) + ") is not defined");
        wR00(x);
    }
  }

  WeightPolynomial poli() const {
    switch(weightId_) {
      case 0: return wRUnitPoli();
      case 1: return wR00Poli();
      case 13: return wR2Poli();
      default:
        throw std::invalid_argument("weight wR(" + std::to_string(weightId_) + ") is not defined");
        wR00Poli();
    }
  }

  static cmplx wTau(cmplx x) {
    return wD00(x);
  }
  // Polunomial weights
  static WeightPolynomial wRUnitPoli() {
    return { 1, 0, 0, 0 };
  }
  static cmplx wRUnit() {
    return 1.0;
  }
  static cmplx wDUnit(cmplx x) {
    return 1.0-2.0*x;
  }
  static WeightPolynomial wR2Poli() {
    return { 1, 0, -1, 0 };
  }
  static cmplx wR2(cmplx x) {
    return 1.0-pow(x, 2);
  }
  static cmplx wD2(cmplx x) {
    return 1.0 - 2.0*x + 2.0*pow(x, 3)/3.0;
  }
  // weightId: 1
  static WeightPolynomial wR00Poli() {
    return { 1, 0, -3, 2 };
  }
  static cmplx wR00(cmplx x) {
    return pow((1. - x), 2)*(1. + 2.*x);
  }
  static cmplx wD00(cmplx x) {
    return pow((1. - x), 3)*(1. + x);
  }
  // weightId: 2
  static cmplx wR01(cmplx x) {
    return pow(1.0 - x, 2)*(1.0 + 2.0*x)*x;
  }
  static cmplx wD01(cmplx x) {
    return pow(1.0 - x, 3)*(3.0 + 9.0*x + 8.0*pow(x, 2))/10.0;
  }
  // weightId: 3
  static cmplx wR02(cmplx x) {
    return pow(1.0 - x, 2)*(1.0 + 2.0*x)*pow(x, 2);
  }
  static cmplx wD02(cmplx x) {
    return 2.0*pow(1.0 - x, 3)*(1.0 + 3.0*x + 6.0*pow(x, 2) + 5.0*pow(x, 3))/15.0;
  }
  // weightId: 4
  static cmplx wR03(cmplx x) {
    return pow(1.0 - x, 2)*(1.0 + 2.0*x)*pow(x, 3);
  }
  static cmplx wD03(cmplx x) {
    return pow(1.0 - x, 3)*(1.0 + 3.0*x + 6.0*pow(x, 2) + 10.0*pow(x, 3) + 8.0*pow(x, 4))/14.0;
  }
  static cmplx wR10(cmplx x) {
    return pow(1.0 - x, 3)*(1.0 + 2.0*x);
  }
  static cmplx wD10(cmplx x) {
    return pow(1.0 - x, 4)*(7.0 + 8.0*x)/10.0;
  }
  static cmplx wR11(cmplx x) {
    return pow(1.0 - x, 3)*(1.0 + 2.0*x)*x;
  }
  static cmplx wD11(cmplx x) {
    return pow(1.0 - x, 4)*pow(1.0 + 2.0*x, 2)/6.0;
  }
  static cmplx wR12(cmplx x) {
    return pow(1.0 - x, 3)*(1.0 + 2.0*x)*pow(x, 2);
  }
  static cmplx wD12(cmplx x) {
    return pow(1.0 - x, 4)*(13.0 + 52.0*x + 130.0*pow(x, 2) + 120.0*pow(x, 3))/210.0;
  }
  static cmplx wR13(cmplx x) {
    return pow(1.0 - x, 3)*(1.0 + 2.0*x)*pow(x, 3);
  }
  static cmplx wD13(cmplx x) {
    return pow(1.0 - x, 4)*(2.0 + 8.0*x + 20.0*pow(x, 2) + 40.0*pow(x, 3) + 35.0*pow(x, 4))/70.0;
  }
  static cmplx wR20(cmplx x) {
    return pow(1.0 - x, 4)*(1.0 + 2.0*x);
  }
  static cmplx wD20(cmplx x) {
    return 2.0*pow(1.0 - x, 5)*(4.0 + 5.0*x)/15.0;
  }
  static cmplx wR21(cmplx x) {
    return pow(1.0 - x, 4)*(1.0 + 2.0*x)*x;
  }
  static cmplx wD21(cmplx x) {
    return pow(1.0 - x, 5)*(11.0 + 55.0*x + 60.0*pow(x, 2))/105.0;
  }
  static cmplx wR22(cmplx x) {
    return pow(1.0 - x, 4)*(1.0 + 2.0*x)*pow(x, 2);
  }
  static cmplx wD22(cmplx x) {
    return pow(1.0 - x, 5)*(1.0 + 5.0*x + 15.0*pow(x, 2) + 15.0*pow(x, 3))/30.0;
  }
  static cmplx wR23(cmplx x) {
    return pow(1.0 - x, 4)*(1.0 + 2.0*x)*pow(x, 3);
  }
  static cmplx wD23(cmplx x) {
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
