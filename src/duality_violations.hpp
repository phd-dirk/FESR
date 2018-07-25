#ifndef SRC_DUALITY_VIOLATIONS_H
#define SRC_DUALITY_VIOLATIONS_H

#include "./types.hpp"
#include "./numerics.hpp"
#include "./weights.hpp"

class DualityViolations : Numerics
{
 public:
  double DVMomentVpA(cDbl &s0, const Weight &w,
                cDbl &vKappa, cDbl &vGamma, cDbl &vAlpha, cDbl &vBeta,
                cDbl &aKappa, cDbl &aGamma, cDbl &aAlpha, cDbl &aBeta) const
  {
    return moment(s0, w, vKappa, vGamma, vAlpha, vBeta)
        + moment(s0, w, aKappa, aGamma, aAlpha, aBeta);
  }
  double moment(cDbl &s0, const Weight &w, cDbl &kappa,
                cDbl &gamma, const cDbl &alpha, const cDbl &beta) const
  {
    func f = [&](double s) -> double {
      return w.wR(s).real()/s0*model(s, kappa, gamma, alpha, beta);
    };
    return -M_PI*semiInfInt(f, s0);
  }

  double model(double s, double kappa, double gamma, double alpha, double beta)
  const {
    return kappa*exp(-gamma*s)*sin(alpha+beta*s);
  }
};

#endif
