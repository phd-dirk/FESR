#ifndef SRC_MQ_RUN_H
#define SRC_MQ_RUN_H

#include "types.hpp"
#include "alpha_s.hpp"

class MQRun {
 public:
  MQRun(const double &sTau) : sTau_(sTau) {}

  /*
    Calculates the ratio m(q^2)/m(p^2) from integrating the RG-equation in the complex
    q^2 plane from a given a(p^2) at p(^2)
  */
  cmplx operator ()(const cmplx &q2, const cmplx &p2, const double &atau) const {
    cmplx ap = AlphaS::run(p2, sTau_, atau);
    cmplx aq = AlphaS::run(q2, sTau_, atau);

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

 private:
  const double sTau_;
};

#endif
