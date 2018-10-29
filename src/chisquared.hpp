#ifndef SRC_CHISQUARED_H
#define SRC_CHISQUARED_H

#include "./types.hpp"
#include "./configuration.hpp"
#include "./experimentalMoments.hpp"
#include "./theoretical_moments.hpp"
#include "./weights.hpp"
#include "./utils.hpp"
#include <nlohmann/json.hpp>
#include <functional>
#include <iostream>
#include <thread>
#include <future>
#include <chrono>

using std::function;

class Chi2 {
 public:
  Chi2(Configuration config);

  double operator ()( const double *xx) const;
  double operator ()(
    const double &astau,
    const double &aGGinv,
    const double &rho,
    const double &c8,
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
    const double &c8VpA,
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

  const Configuration config_;
  const ExpMoms expMom_;
  matrix<double> c_;
};

#endif
