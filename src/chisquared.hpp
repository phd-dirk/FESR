#ifndef SRC_CHISQUARED_H
#define SRC_CHISQUARED_H

#include "./types.hpp"
#include "./configuration.hpp"
#include "./experimentalMoments.hpp"
#include "./theoretical_moments.hpp"
#include "./numerics.hpp"
#include "./weights.hpp"
#include "./utils.hpp"
#include <nlohmann/json.hpp>
#include <functional>
#include <iostream>
#include <thread>
#include <future>
#include <chrono>

using std::function;
using std::cout;
using std::endl;

class Chisquared: Numerics {
 public:
  Chisquared(Configuration config);

  double operator ()( const double *xx) const;
  double operator ()(
    std::vector<Input> inputs,
    const double &astau,
    const double &aGGinv,
    const double &rho,
    const double &c8,
    const double &vKappa,
    const double &vGamma,
    const double &vAlpha,
    const double &vBeta,
    const double &aKappa,
    const double &aGamma,
    const double &aAlpha,
    const double &aBeta
  ) const;

  mat invCovMat_;
private:
  double chi2(
    std::vector<Input> inputs,
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
    std::vector<Input> inputs,
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

  void initInvCovMat();

  const Configuration config_;
  const std::vector<Input> inputs_;
  const int order_;
  int momCount_;
  const ExperimentalMoments expMom_;
};

#endif
