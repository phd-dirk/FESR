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

  double model(double s, double kappa, double gamma, double alpha, double beta) const
  {
    cout << exp(-kappa - gamma*s)*sin(alpha + beta*s) << "\t" << s << "\t" << kappa << "\t" << gamma << "\t" << alpha << "\t" << beta << endl;
    // 3.8787504179628975e-05	3514.4831785351275	3	0.0018963380072094527	-2.2000000000000002	3.8999999999999999
    return exp(-kappa-gamma*s)*sin(alpha+beta*s);
  }
};

#endif
