#ifndef SRC_ALPHAS_H
#define SRC_ALPHAS_H

#include "./constants.hpp"
#include "numerics.hpp"
#include <complex>
#include <iostream>
#include <functional>
#include <cmath>

using std::complex;
using std::function;
using std::abs;
using std::cout;
using std::endl;
using std::sqrt;
using std::atan;
using std::log;

class AlphaS: public Numerics {
 public:
  AlphaS(const Constants &constants) : Numerics(constants), const_(constants) {}

  complex<double> operator ()(const complex<double> &q2, const double &p2, const double &ap) {
    return zarg(q2, p2, ap);
  }

  complex<double> operator ()(const complex<double> &q2,const complex<double> &p2, const complex<double> &ap) {
    return zarg(q2, p2, ap);
  }

  complex<double> zarg(const complex<double> &q2, const complex<double> &p2, const complex<double> &ap) {
    // integrate beta function and find root
    // from mathematica coefficients.nb
    auto f = [&ap, &p2, &q2](complex<double> a2) {
      return 0.2222222222222/ap - 0.2222222222222/a2 - 0.3730002895803*atan(0.19576224733469 - 2.7775209170642*ap) + 0.3730002895803*atan(0.19576224733469 - 2.777520917064*a2) +
      0.3950617283951*log(ap) - 0.2727626771781*log(0.353968700519 + 1.*ap) - 0.0611495256085*log(0.13459153249825 - 0.14096185280333*ap + 1.*pow(ap,2)) -
      0.3950617283951*log(a2) + 0.2727626771781*log(0.353968700519 + 1.*a2) + 0.0611495256085*log(0.13459153249825 - 0.14096185280333*a2 + 1.*pow(a2,2)) - log(p2/q2);
    };
    auto df = [](complex<double> a2) {
      return -1.03601610638/(1.0 + pow(0.19576224733469 - 2.777520917064*a2,2)) + 0.2222222222222/pow(a2,2) - 0.3950617283951/a2 + 0.2727626771781/(0.353968700519 + 1.*a2) +
      (0.0611495256085*(-0.14096185280333 + 2.*a2))/(0.13459153249825 - 0.14096185280333*a2 + 1.*pow(a2,2));
    };

    return newtonRaphson(f, df, complex<double>(0.01, 0.01), 1e-15);
  }

 private:
  Constants const_;
};

#endif
