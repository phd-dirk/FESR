#ifndef SRC_THEORETICAL_MOMENTS_H
#define SRC_THEORETICAL_MOMENTS_H

#include "./adler_function.hpp"

class TheoreticalMoments: public AdlerFunction {
 public:
  TheoreticalMoments(const int &order, const vector<double> &s0s,
                     function<complex<double>(complex<double>)> weight, Constants constants) :
    AdlerFunction(order, constants), const_(constants), s0s(s0s), weight(weight) {}

  double operator ()(const int &i, const double &astau, const double &aGGinv,
                     const double &rhoVpA, const double &c8VpA, const double &order) {
    // factor 2 for V PLUS A
    double s0 = s0s[i];
    // d0 VpA
    double rTauTh = cIntVpAD0FO(s0, weight, astau, order)
      // d4 VpA
      + cIntVpAD4FO(s0, weight, astau, aGGinv, order)
      // d68 VpA
       + D68CInt(s0, weight, rhoVpA, c8VpA)
       // + D0CInt(s0, weight, astau, 0)*const_.deltaEW
      // deltaP
       + 3.*deltaP(s0, wR00);
    return pow(const_.kVud, 2)*const_.kSEW*rTauTh;
  }

  double cIntVpAD0FO(const double &s0,
                     function<complex<double>(complex<double>)> weight,
                     const double &astau, const int &order) {
    return 2.*D0CInt(s0, weight, astau, order);
  }

  double cIntVpAD4FO(const double &s0,
                     function<complex<double>(complex<double>)> weight,
                     const double &astau, const double &aGGinv, const int &order) {
    return  D4CInt(s0, weight, astau, aGGinv, order, 1) + D4CInt(s0, weight, astau, aGGinv, order, -1);
  }

  double del0(const double &s0,
              function<complex<double>(complex<double>)> weight,
              const double &astau, const int &order) {
    return (cIntVpAD0FO(s0, weight, astau, order)
            - cIntVpAD0FO(s0, weight, astau, 0)
            )/3.0;
  }

  double del2(const double &s0,
              function<complex<double>(complex<double>)> weight,
              const double &astau, const int &order) {
    return (D2CInt(s0, weight, astau, order, 1)
            + D2CInt(s0, weight, astau, order, -1)
            )/3.;
  }

  double del4(const double &s0,
              function<complex<double>(complex<double>)> weight,
              const double &astau, const double &aGGinv, const int &order) {
    return cIntVpAD4FO(s0, weight, astau, aGGinv, order)/3.;
  }

  double del68(const double &s0, function<complex<double>(complex<double>)> weight,
               const double &rhoVpA, const double &c8VpA) {
    return D68CInt(s0, weight, rhoVpA, c8VpA)/3.;
  }


 private:
  Constants const_;
  vector<double> s0s;
  function<complex<double>(complex<double>)> weight;
};

#endif
