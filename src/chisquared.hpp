#ifndef SRC_CHISQUARED_H
#define SRC_CHISQUARED_H

#include "./types.hpp"
#include "./configuration.hpp"
#include "./experimentalMoments.hpp"
#include "./theoretical_moments.hpp"
#include "./weights.hpp"
#include "./utils.hpp"
#include "./data.hpp"
#include <nlohmann/json.hpp>
#include <functional>
#include <iostream>
#include <thread>
#include <future>
#include <chrono>
#include <boost/numeric/ublas/matrix.hpp>

using std::function;

class Chi2 {
 public:
  Chi2(Configuration config);
  Chi2(
    const int &nc,
    const int &nf,
    const double &order,
    const double &sTau,
    const double &mPiM,
    const double &fPi,
    const double &dFPi,
    const double &f1P,
    const double &m1P,
    const double &g1P,
    const double &f2P,
    const double &m2P,
    const double &g2P,
    const std::vector<double> &mq,
    const double &be,
    const double &dBe,
    const double &vud,
    const double &dVud,
    const double &SEW,
    const double &dSEW,
    const double &asTauFix,
    const std::vector<double> &qqMTau,
    const std::vector<Input> &inputs,
    const ThMomContribs &thMomContribs,
    const double &RVANormalization
  );

  double operator ()( const double *xx) const;
  double operator ()(
    const double &astau,
    const double &aGGinv,
    const double &rho,
    const double &c8,
    const double &c10,
    const double &c12,
    const double &vDelta,
    const double &vGamma,
    const double &vAlpha,
    const double &vBeta,
    const double &aDelta,
    const double &aGamma,
    const double &aAlpha,
    const double &aBeta,
    const double &order
  ) const;

  static double chi2(
    const double &astau,
    const double &aGGinv,
    const double &rho,
    const double &c8,
    const double &c10,
    const double &c12,
    const double &order,
    const double &deV,
    const double &gaV,
    const double &alV,
    const double &beV,
    const double &deA,
    const double &gaA,
    const double &alA,
    const double &beA,
    const double &sTau,
    const double &mPiM,
    const double &fPi,
    const double &f1P,
    const double &m1P,
    const double &g1P,
    const double &f2P,
    const double &m2P,
    const double &g2P,
    const matrix<double> &c,
    const std::vector<double> &mq,
    const double &vud,
    const double &SEW,
    const Condensates &condensates,
    const std::vector<Input> &inputs,
    const ThMomContribs &thMomContribs,
    const ExpMoms &expMoms
  );

  static std::vector<double> calcThMoms(
    const double &astau,
    const double &aGGinv,
    const double &rhoVpA,
    const double &c8,
    const double &c10,
    const double &c12,
    const double &order,
    const double &deV,
    const double &gaV,
    const double &alV,
    const double &beV,
    const double &deA,
    const double &gaA,
    const double &alA,
    const double &beA,
    const double &sTau,
    const double &mPiM,
    const double &fPi,
    const double &f1P,
    const double &m1P,
    const double &g1P,
    const double &f2P,
    const double &m2P,
    const double &g2P,
    const matrix<double> &c,
    const std::vector<double> &mq,
    const double &vud,
    const double &SEW,
    const Condensates &condensates,
    const std::vector<Input> &inputs,
    const ThMomContribs &thMomContribs
  );

  double chi2_spec_end(const double *xx) const;

  const ExpMoms expMom_;
  matrix<double> c_;
  double order_, sTau_, mPiM_, fPi_, vud_, SEW_,
    f1P_, m1P_, g1P_, f2P_, m2P_, g2P_;
  std::vector<double> mq_;
  Condensates condensates_;
  std::vector<Input> inputs_;
  ThMomContribs thMomContribs_;
};

#endif
