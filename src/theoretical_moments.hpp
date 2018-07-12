#ifndef SRC_THEORETICAL_MOMENTS_H
#define SRC_THEORETICAL_MOMENTS_H

#include "./types.hpp"
#include "./configuration.hpp"
#include "./adler_function.hpp"


class TheoreticalMoments: public AdlerFunction {
 public:
  TheoreticalMoments(const Configuration &config) :
    AdlerFunction(config), config_(config), inputs_(config_.inputs) {}

  vec operator ()(const double &astau,
                     const double &aGGinv, const double &rhoVpA, const double &c8VpA,
                     const double &order) const {
    return thMoms(astau, aGGinv, rhoVpA, c8VpA, order);
  }

  double thMom(const double &s0, const Weight &w, const double &astau,
        const double &aGGinv, const double &rhoVpA, const double &c8VpA,
        const double &order) const {
    double rTauTh = 0.;
    // D0
    if ( config_.OPE.D0 ) {
      // check if FOPT or CIPT
      if ( config_.OPE.scheme == "FO" ) {
        rTauTh += cIntVpAD0FO(s0, w, config_.sTau, astau, order);
      }
      if ( config_.OPE.scheme == "CI") {
        rTauTh += cIntVpAD0CI(s0, w, config_.sTau, astau, order);
      }
    }
    // D4
    if ( config_.OPE.D4 )
      rTauTh += cIntVpAD4FO(s0, w, config_.sTau, astau, aGGinv);
    // D68
    if ( config_.OPE.D68 )
      rTauTh += D68CInt(s0, w, rhoVpA, c8VpA);
    // PionPole
    if ( config_.OPE.PionPole )
      rTauTh += 3.*deltaP(s0, w);

    return pow(config_.kVud, 2)*config_.kSEW*rTauTh;
  }

  vec thMoms(const double &astau, const double &aGGinv, const double &rhoVpA,
              const double &c8VpA, const double &order) const {
    vec moms;
    for(auto const& input: inputs_) {
      vec s0s = input.s0s;
      Weight w = input.weight;
      for(auto const& s0: s0s) {
        moms.push_back(thMom(s0, w, astau, aGGinv, rhoVpA, c8VpA, order));
      }
    }
    return moms;
  }

  double cIntVpAD0FO(const double &s0, const Weight &weight, const double &sTau,
                     const double &astau, const int &order) const {
    return 2.*D0CIntFO(s0, weight, sTau, astau, order);
  }
  double cIntVpAD0CI(const double &s0, const Weight &weight, const double &sTau,
                     const double &astau, const int &order) const {
    return 2.*D0CIntCI(s0, weight, sTau, astau, order);
  }


  double cIntVpAD4FO(const double &s0, const Weight &weight, const double &sTau,
                     const double &astau, const double &aGGinv) const {
    return  D4CInt(s0, weight, sTau, astau, aGGinv, 1) + D4CInt(s0, weight, sTau, astau, aGGinv, -1);
  }
  double del0(const double &s0, const Weight &weight,
              const double &sTau, const double &astau, const int &order) const {
    return (cIntVpAD0FO(s0, weight, sTau, astau, order)
            - cIntVpAD0FO(s0, weight, sTau, astau, 0)
            )/3.0;
  }

  double del2(const double &s0, const Weight &weight,
              const double &astau, const int &order) const {
    return (D2CInt(s0, weight, astau, order, 1)
            + D2CInt(s0, weight, astau, order, -1)
            )/3.;
  }

  double del4(const double &s0, const Weight &weight, const double &sTau,
              const double &astau, const double &aGGinv) const {
    return cIntVpAD4FO(s0, weight, sTau, astau, aGGinv)/3.;
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
    // cout << "thMom: \t" << operator() (config_.sTau, astau, aGGinv, rhoVpA, c8VpA, order) << endl;
    // cout << "Delta^(0): \t" << del0(config_.sTau, config_.weight, config_.sTau, astau, order) << endl;
    // cout << "Delta^(4): \t" << del4(config_.sTau, config_.weight, config_.sTau, astau, aGGinv) << endl;
    // cout << "Delta^(6): \t" << del6(config_.sTau, config_.weight, rhoVpA) << endl;
    // cout << "Delta^(8): \t" << del8(config_.sTau, config_.weight, c8VpA) << endl;
    // cout << "Delta_P+S: \t" << deltaP(config_.sTau, config_.weight) << endl;
  }


 private:
  Configuration config_;
  std::vector<Input> inputs_;
};

#endif
