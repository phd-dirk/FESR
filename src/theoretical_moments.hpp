#ifndef SRC_THEORETICAL_MOMENTS_H
#define SRC_THEORETICAL_MOMENTS_H

#include "./adler_function.hpp"
#include "json.hpp"

using json = nlohmann::json;

class TheoreticalMoments: public AdlerFunction {
 public:
  TheoreticalMoments(const json &config) :
    AdlerFunction(config["parameters"]["alphaLoops"], Constants(config)), config_(config), const_(Constants(config)),
    s0s_(config["parameters"]["s0Set"].get<vector<double>>()), weight_(Weight(config["parameters"]["weight"].get<int>())) {
  }

  double operator ()(const int &i, const double &astau, const double &aGGinv,
                     const double &rhoVpA, const double &c8VpA, const double &order) const {
    double s0 = s0s_[i];

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
                     const double &astau, const int &order) const {
    return 2.*D0CInt(s0, weight, astau, order);
  }

  double cIntVpAD4FO(const double &s0, const Weight &weight,
                     const double &astau, const double &aGGinv, const int &order) const {
    return  D4CInt(s0, weight, astau, aGGinv, order, 1) + D4CInt(s0, weight, astau, aGGinv, order, -1);
  }

  double del0(const double &s0, const Weight &weight,
              const double &astau, const int &order) const {
    return (cIntVpAD0FO(s0, weight, astau, order)
            - cIntVpAD0FO(s0, weight, astau, 0)
            )/3.0;
  }

  double del2(const double &s0, const Weight &weight,
              const double &astau, const int &order) const {
    return (D2CInt(s0, weight, astau, order, 1)
            + D2CInt(s0, weight, astau, order, -1)
            )/3.;
  }

  double del4(const double &s0, const Weight &weight,
              const double &astau, const double &aGGinv, const int &order) const {
    return cIntVpAD4FO(s0, weight, astau, aGGinv, order)/3.;
  }

  double del6(const double &s0, const Weight &weight,
               const double &rhoVpA) const {
    return D68CInt(s0, weight, rhoVpA, 0.0)/3.;
  }

  double del8(const double &s0, const Weight &weight,
              const double &c8VpA) const {
    return D68CInt(s0, weight, 0.0, c8VpA)/3.;
  }

  double del68(const double &s0, const Weight &weight,
               const double &rhoVpA, const double &c8VpA) {
    return D68CInt(s0, weight, rhoVpA, c8VpA)/3.;
  }

  void log(const double &astau, const double &aGGinv, const double &rhoVpA, const double &c8VpA, const int &order) const {
    double s0 = config_["constants"]["sTau"];
    cout << "thMom: \t" << operator() (0, astau, aGGinv, rhoVpA, c8VpA, order) << endl;
    cout << "Delta^(0): \t" << del0(s0, weight_, astau, order) << endl;
    cout << "Delta^(4): \t" << del4(s0, weight_, astau, aGGinv, order) << endl;
    cout << "Delta^(6): \t" << del6(s0, weight_, rhoVpA) << endl;
    cout << "Delta^(8): \t" << del8(s0, weight_, c8VpA) << endl;
    cout << "Delta_P+S: \t" << deltaP(s0, weight_) << endl;
  }


 private:
  json config_;
  Constants const_;
  vector<double> s0s_;
  Weight weight_;
};

#endif
