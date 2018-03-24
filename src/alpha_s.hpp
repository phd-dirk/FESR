#ifndef SRC_ALPHAS_H
#define SRC_ALPHAS_H

#include "./constants.hpp"
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

inline complex<double> alpha_s(const complex<double> &s, const double &astau) {
  complex<double> I(0.0, 1.0);
  // double a1 = 0.10162679736189885;
  Constants constants(3, 3);
  double a1 = astau/constants.kPi;
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

/*
  Calculates the ratio m(q^2)/m(p^2) from integrating the RG-equation in the complex
  q^2 plane from a given a(p^2) at p(^2)
 */
inline complex<double> runMassRatio(const complex<double> &q2,
                                    const complex<double> &p2, const double &astau) {
  complex<double> I(0.0, 1.0);
  complex<double> ap = alpha_s(sqrt(p2), astau);
  complex<double> aq = alpha_s(sqrt(q2), astau);

  auto f = [](complex<double> ap, complex<double> aq) {
    return 0.25000289589113*atan(0.195762247334686 - 2.77752091706421*aq)
    - 0.25000289589113*atan(0.195762247334686 - 2.77752091706421*ap)
    - 0.444444444444444*log(aq)
    - 0.144635029127131*log(0.353968700519028 + 1.*aq)
    - 0.174068624534112*log(0.134591532498253 - 0.140961852803328*aq
                            + 1.*pow(aq, 2))
    + 0.444444444444444*log(ap) + 0.144635029127131*log(0.353968700519028 + 1.*ap)
    + 0.174068624534112*log(0.134591532498253 - 0.140961852803328*ap
                            + 1.*pow(ap, 2));
  };

  return exp(f(ap, aq));
}


#endif
