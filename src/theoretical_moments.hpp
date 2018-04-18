#ifndef SRC_THEORETICAL_MOMENTS_H
#define SRC_THEORETICAL_MOMENTS_H

#include "./adler_function.hpp"
#include "json.hpp"

using json = nlohmann::json;

class TheoreticalMoments: public AdlerFunction {
 public:
  TheoreticalMoments(const int &order, const vector<double> &s0s, const Weight &weight,
                     const json &configuration, const Constants &constants) :
    AdlerFunction(order, constants), const_(constants), config_(configuration),
    s0s(s0s), weight_(weight) {}

  double operator ()(const int &i, const double &astau, const double &aGGinv,
                     const double &rhoVpA, const double &c8VpA, const double &order) {
    double s0 = s0s[i];

    double rTauTh = 0.;
    // D0
    if ( config_["adler"]["D0"] )
      rTauTh += cIntVpAD0FO(s0, weight_, astau, order);
    // D4
    if ( config_["adler"]["D4"] )
      rTauTh += cIntVpAD4FO(s0, weight_, astau, aGGinv, order);
    // D68
    if ( config_["adler"]["D68"] )
      rTauTh += D68CInt(s0, weight_, rhoVpA, c8VpA);
    // PionPole
    if ( config_["adler"]["PionPole"] )
      rTauTh += 3.*deltaP(s0, weight_);

    return pow(const_.kVud, 2)*const_.kSEW*rTauTh;
  }

  double cIntVpAD0FO(const double &s0, const Weight &weight,
                     const double &astau, const int &order) {
    return 2.*D0CInt(s0, weight, astau, order);
  }

  double cIntVpAD4FO(const double &s0, const Weight &weight,
                     const double &astau, const double &aGGinv, const int &order) {
    return  D4CInt(s0, weight, astau, aGGinv, order, 1) + D4CInt(s0, weight, astau, aGGinv, order, -1);
  }

  double del0(const double &s0, const Weight &weight,
              const double &astau, const int &order) {
    return (cIntVpAD0FO(s0, weight, astau, order)
            - cIntVpAD0FO(s0, weight, astau, 0)
            )/3.0;
  }

  double del2(const double &s0, const Weight &weight,
              const double &astau, const int &order) {
    return (D2CInt(s0, weight, astau, order, 1)
            + D2CInt(s0, weight, astau, order, -1)
            )/3.;
  }

  double del4(const double &s0, const Weight &weight,
              const double &astau, const double &aGGinv, const int &order) {
    return cIntVpAD4FO(s0, weight, astau, aGGinv, order)/3.;
  }

  double del68(const double &s0, const Weight &weight,
               const double &rhoVpA, const double &c8VpA) {
    return D68CInt(s0, weight, rhoVpA, c8VpA)/3.;
  }

  void log(const double &astau, const int &order) {
    cout << "Delta_0: \t" << del0(s0s[0], weight_, astau, order) << endl;
  }


 private:
  Constants const_;
  json config_;
  vector<double> s0s;
  Weight weight_;
};

#endif
