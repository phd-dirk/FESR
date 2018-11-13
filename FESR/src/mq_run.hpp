#ifndef SRC_MQ_RUN_H
#define SRC_MQ_RUN_H

#include "types.hpp"
#include "alpha_s.hpp"

class MQ {
 public:
  /*
    Calculates the ratio m(q^2)/m(p^2) from integrating the RG-equation in the complex
    q^2 plane from a given a(p^2) at p(^2)
  */
  static complex<double> run(
    const complex<double> &q2, const complex<double> &p2,
    const double &sTau, const double &atau
  ) {
    cmplx ap = AlphaS::run(p2, sTau, atau);
    cmplx aq = AlphaS::run(q2, sTau, atau);

    auto f = [](cmplx ap, cmplx aq) {
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
};

#endif
