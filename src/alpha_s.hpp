#ifndef SRC_ALPHAS_H
#define SRC_ALPHAS_H

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

// find root of function
inline complex<double> newtonRaphson(
    function<complex<double>(complex<double>)> f,
    function<complex<double>(complex<double>)> df,
    complex<double> x0,
    double acc ) {
  complex<double> xNext = x0;
  while( abs(f(xNext)) > acc ) {
    xNext = xNext - f(xNext)/df(xNext);
  }
  return xNext;
}

inline complex<double> alpha_s(complex<double> s) {
  complex<double> I(0.0, 1.0);
  double a1 = 0.10162679736189885;
  double mu1 = sqrt(3.1570893124000001);
  complex<double> mu2 = s;

  // integrate beta function and find root
  // from mathematica coefficients.nb
  auto f = [&a1, &mu1, &mu2](complex<double> a2) {
    return 0.2222222222222/a1 - 0.2222222222222/a2 - 0.3730002895803*atan(0.19576224733469 - 2.7775209170642*a1) + 0.3730002895803*atan(0.19576224733469 - 2.777520917064*a2) +
    0.3950617283951*log(a1) - 0.2727626771781*log(0.353968700519 + 1.*a1) - 0.0611495256085*log(0.13459153249825 - 0.14096185280333*a1 + 1.*pow(a1,2)) -
    0.3950617283951*log(a2) + 0.2727626771781*log(0.353968700519 + 1.*a2) + 0.0611495256085*log(0.13459153249825 - 0.14096185280333*a2 + 1.*pow(a2,2)) - log(mu1/mu2);
  };
  auto df = [](complex<double> a2) {
      return -1.03601610638/(1.0 + pow(0.19576224733469 - 2.777520917064*a2,2)) + 0.2222222222222/pow(a2,2) - 0.3950617283951/a2 + 0.2727626771781/(0.353968700519 + 1.*a2) + 
      (0.0611495256085*(-0.14096185280333 + 2.*a2))/(0.13459153249825 - 0.14096185280333*a2 + 1.*pow(a2,2));
  };

  return newtonRaphson(f, df, complex<double>(0.01, 0.01), 1e-14);
}

// inline complex<double> runningMass(complex<double> mu) {
//   complex<double> I(0.0, 1.0);
//   double a1 = 0.10162679736189885;
//   double m1 = 
//   auto f = []
// }


#endif
