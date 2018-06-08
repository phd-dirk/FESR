#ifndef SRC_MQ_RUN_H
#define SRC_MQ_RUN_H

#include "constants.hpp"
#include "alpha_s.hpp"
#include <complex>
#include <functional>

using std::complex;
using std::function;

class MQRun {
 public:
  MQRun(const Constants &constants) : const_(constants), amu_(constants) {}

  complex<double> operator ()(const complex<double> &q2, const complex<double> &p2, const double &atau) const {
    return runMassRatio(q2, p2, atau);
  }

  /*
    Calculates the ratio m(q^2)/m(p^2) from integrating the RG-equation in the complex
    q^2 plane from a given a(p^2) at p(^2)
  */
  complex<double> runMassRatio(const complex<double> &q2,
                               const complex<double> &p2, const double &atau) const {
    complex<double> I(0.0, 1.0);
    complex<double> ap = amu_(p2, const_.kSTau, atau);
    complex<double> aq = amu_(q2, const_.kSTau, atau);

    auto f = [](complex<double> ap, complex<double> aq) {
      return 0.25000289589113*atan(0.195762247334686 - 2.77752091706421*ap)
      - 0.25000289589113*atan(0.195762247334686 - 2.77752091706421*aq)
      - 0.444444444444444*log(ap)
      - 0.144635029127131*log(0.353968700519028 + 1.*ap)
      - 0.174068624534112*log(0.134591532498253 - 0.140961852803328*ap
                              + 1.*pow(ap, 2))
      + 0.444444444444444*log(aq) + 0.144635029127131*log(0.353968700519028 + 1.*aq)
      + 0.174068624534112*log(0.134591532498253 - 0.140961852803328*aq
                              + 1.*pow(aq, 2));
    };

    return exp(f(ap, aq));
  }
 private:
  Constants const_;
  AlphaS amu_;
};

#endif
