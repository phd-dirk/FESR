#ifndef SRC_DUALITY_VIOLATIONS_H
#define SRC_DUALITY_VIOLATIONS_H

#include "./types.hpp"
#include "./numerics.hpp"
#include "./weights.hpp"

class DualityViolations : Numerics
{
 public:
  double DVMomentVpA(cDbl &s0, const Weight &w,
                cDbl &deV, cDbl &gaV, cDbl &alV, cDbl &beV,
                cDbl &deA, cDbl &gaA, cDbl &alA, cDbl &beA) const
  {
    return -12.0*pow(M_PI, 2)*(
      cintDVp_VA(s0, w, deV, gaV, alV, beV)
      +cintDVp_VA(s0, w, deA, gaA, alA, beA)
    );
  }

  double cintDVp_VA(cDbl &s0, const Weight &w, cDbl &de, cDbl &ga, cDbl &al, cDbl &be) const
  {
    return w.poli().x0*intDVp0(s0, de, ga, al, be)
        + w.poli().x1*intDVp0(s0, de, ga, al, be)
        + w.poli().x2*intDVp2(s0, de, ga, al, be)
        + w.poli().x3*intDVp3(s0, de, ga, al, be);
  }

  double rhoDV(cDbl &s, cDbl &delta, cDbl &gamma, cDbl &alpha, cDbl &beta) const
  {
    return exp(-delta-gamma*s)*sin(alpha+beta*s);
  }

  double intDVp0(cDbl &s0, cDbl &delta, cDbl &gamma, cDbl &alpha, cDbl &beta) const
  {
    return (exp(-delta - gamma*s0)*(beta*cos(alpha + beta*s0) + gamma*sin(alpha + beta*s0)))/(pow(beta,2) + pow(gamma,2));
  }

  double intDVp1(cDbl &s0, cDbl &delta, cDbl &gamma, cDbl &alpha, cDbl &beta) const
  {
    return (exp(-delta - gamma*s0)*(beta*(2*gamma + (pow(beta,2) + pow(gamma,2))*s0)*cos(alpha + beta*s0) + (pow(gamma,2) + pow(gamma,3)*s0 + pow(beta,2)*(-1 + gamma*s0))*sin(alpha + beta*s0)))/
        pow(pow(beta,2) + pow(gamma,2),2);
  }

  double intDVp2(cDbl &s0, cDbl &delta, cDbl &gamma, cDbl &alpha, cDbl &beta) const
  {
    return (exp(-delta - gamma*s0)*(beta*(-2*pow(beta,2) + 6*pow(gamma,2) + 4*gamma*(pow(beta,2) + pow(gamma,2))*s0 + pow(pow(beta,2) + pow(gamma,2),2)*pow(s0,2))*cos(alpha + beta*s0) +(2*gamma*(-3*pow(beta,2) + pow(gamma,2)) + 2*(-pow(beta,4) + pow(gamma,4))*s0 + gamma*pow(pow(beta,2) + pow(gamma,2),2)*pow(s0,2))*sin(alpha + beta*s0)))/
        pow(pow(beta,2) + pow(gamma,2),3);
  }

  double intDVp3(cDbl &s0, cDbl &delta, cDbl &gamma, cDbl &alpha, cDbl &beta) const
  {
    return (exp(-delta - gamma*s0)*(beta*(pow(beta,6)*pow(s0,3) + 3*pow(beta,4)*s0*(-2 + 2*gamma*s0 + pow(gamma,2)*pow(s0,2)) + 
                                          3*pow(beta,2)*gamma*(-8 + 4*gamma*s0 + 4*pow(gamma,2)*pow(s0,2) + pow(gamma,3)*pow(s0,3)) + 
                                              pow(gamma,3)*(24 + 18*gamma*s0 + 6*pow(gamma,2)*pow(s0,2) + pow(gamma,3)*pow(s0,3)))*cos(alpha + beta*s0) + 
                                        (pow(beta,6)*pow(s0,2)*(-3 + gamma*s0) + 3*pow(beta,4)*(2 - 6*gamma*s0 - pow(gamma,2)*pow(s0,2) + pow(gamma,3)*pow(s0,3)) + 
                                         3*pow(beta,2)*pow(gamma,2)*(-12 - 4*gamma*s0 + pow(gamma,2)*pow(s0,2) + pow(gamma,3)*pow(s0,3)) + 
                                         pow(gamma,4)*(6 + 6*gamma*s0 + 3*pow(gamma,2)*pow(s0,2) + pow(gamma,3)*pow(s0,3)))*sin(alpha + beta*s0)))/pow(pow(beta,2) + pow(gamma,2),4);
  }
};

#endif
