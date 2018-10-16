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

class Chisquared {
 public:
  Chisquared(Configuration config);

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
    const double &aBeta
  ) const;

private:
  double chi2(
    const double &astau,
    const double &aGGinv,
    const double &rho,
    const double &c8,
    const double &deV,
    const double &gaV,
    const double &alV,
    const double &beV,
    const double &deA,
    const double &gaA,
    const double &alA,
    const double &beA
  ) const;

  vec calcThMoms(
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
    const double &beA
  ) const;

  const Configuration config_;
  const ExpMoms expMom_;
};

#endif
