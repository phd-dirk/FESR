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
  AlphaS(const Constants &constants, const int &order) : Numerics(constants), order_(order), const_(constants) {}

  complex<double> operator ()(const complex<double> &q2, const double &p2, const double &ap) {
    return zarg(q2, p2, ap);
  }

  complex<double> operator ()(const complex<double> &q2,const complex<double> &p2, const complex<double> &ap) {
    return zarg(q2, p2, ap);
  }

  // Calculates a(q^2) from integrating the RG-equation in the complex q^2 plane from a given a(p^2) at p^2
  complex<double> zarg(const complex<double> &q2, const complex<double> &p2, const complex<double> &ap) {
    // integrate beta function and find root
    // from mathematica coefficients.nb
    function<complex<double>(complex<double>)> f, df;

    if (order_ < 5) {
      // BETA 4th order --------------------------------------------------
      f = [&ap, &p2, &q2](complex<double> a2) {
        return 0.2222222222222/ap - 0.2222222222222/a2 - 0.3730002895803*atan(0.19576224733469 - 2.7775209170642*ap) + 0.3730002895803*atan(0.19576224733469 - 2.777520917064*a2) +
        0.3950617283951*log(ap) - 0.2727626771781*log(0.353968700519 + 1.*ap) - 0.0611495256085*log(0.13459153249825 - 0.14096185280333*ap + 1.*pow(ap,2)) -
        0.3950617283951*log(a2) + 0.2727626771781*log(0.353968700519 + 1.*a2) + 0.0611495256085*log(0.13459153249825 - 0.14096185280333*a2 + 1.*pow(a2,2)) - log(sqrt(p2)/sqrt(q2));
      };
      df = [](complex<double> a2) {
        return -1.03601610638/(1.0 + pow(0.19576224733469 - 2.777520917064*a2,2)) + 0.2222222222222/pow(a2,2) - 0.3950617283951/a2 + 0.2727626771781/(0.353968700519 + 1.*a2) +
        (0.0611495256085*(-0.14096185280333 + 2.*a2))/(0.13459153249825 - 0.14096185280333*a2 + 1.*pow(a2,2));
      };
    } else {
      // BETA 5th order --------------------------------------------------
      f = [&ap, &p2, &q2](complex<double> aq) {
        return 0.2222222222222/ap - 0.2222222222222/aq - 0.3245495122567*atan(0.4691036864089 - 3.258646223822*ap)
        + 0.07445667845942*atan(1.5480055209947 + 4.69913810141*ap) + 0.3245495122567*atan(0.4691036864089 - 3.258646223822*aq)
        - 0.07445667845942*atan(1.5480055209947 + 4.69913810141*aq) + 0.3950617283951*log(ap)
        - 0.02467519993787*log(0.11489632695305 - 0.28791323401694*ap + 1.*pow(ap,2))
        - 0.1728556642597*log(0.1538056626836 + 0.6588465746645*ap + 1.*pow(ap,2)) - 0.3950617283951*log(aq)
        + 0.02467519993787*log(0.11489632695305 - 0.2879132340169*aq + 1.*pow(aq,2))
        + 0.1728556642597*log(0.1538056626836 + 0.6588465746645*aq + 1.*pow(aq,2)) - log(p2/q2);
      };
      df = [](complex<double> aq) {
        return -1.057592042559/(1.0 + pow(0.4691036864089 - 3.258646223822*aq,2)) + 0.2222222222222/pow(aq,2) - 0.3950617283951/aq
        + (0.02467519993787*(-0.2879132340169 + 2.*aq))/(0.11489632695305 - 0.2879132340169*aq + 1.*pow(aq,2))
        + (0.1728556642597*(0.6588465746645 + 2.*aq))/(0.1538056626836 + 0.6588465746645*aq + 1.*pow(aq,2))
        - 0.3498822146531/(1.0 + pow(1.5480055209947 + 4.69913810141*aq,2));
      };
    }

    return newtonRaphson(f, df, complex<double>(0.01, 0.01), 1e-15);
  }

 private:
  int order_;
  Constants const_;
};

#endif
