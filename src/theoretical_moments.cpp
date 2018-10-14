#include "./theoretical_moments.hpp"

TheoreticalMoments::TheoreticalMoments(const Configuration &config)
  : AdlerFunction(config), config_(config), inputs_(config_.inputs_) {}

double TheoreticalMoments::thMom(
  cDbl &s0, const Weight &w,
  cDbl &astau, cDbl &aGGinv, cDbl &rhoVpA, cDbl &c8VpA, cDbl &order,
  cDbl &deV, cDbl &gaV, cDbl &alV, cDbl &beV,
  cDbl &deA, cDbl &gaA, cDbl &alA, cDbl &beA) const
{
  double rTauTh = 0.;
  // D0
  if ( config_.thMomContribs_.D0 ) {
    // check if FOPT or CIPT
    if ( config_.thMomContribs_.scheme == "FO" ) {
      rTauTh += cIntVpAD0FO(s0, w, config_.sTau_, astau, order);
    }
    if ( config_.thMomContribs_.scheme == "CI") {
      rTauTh += cIntVpAD0CI(s0, w, config_.sTau_, astau, order);
    }
  }
  // D4
  if ( config_.thMomContribs_.D4 )
    rTauTh += cIntVpAD4FO(s0, w, config_.sTau_, astau, aGGinv);
  // D68
  if ( config_.thMomContribs_.D68 )
    rTauTh += D68CInt(s0, w, rhoVpA, c8VpA);
  // DV
  if ( config_.thMomContribs_.DV )
    rTauTh += DVMomentVpA(s0, w, deV, gaV, alV, beV, deA, gaA, alA, beA);
  // PionPole
  if ( config_.thMomContribs_.PionPole )
    rTauTh += 3.*deltaP(s0, w);

  return pow(config_.vud_, 2)*config_.kSEW*rTauTh;
}

double TheoreticalMoments::cIntVpAD0FO(
  const double &s0, const Weight &weight, const double &sTau,
  const double &astau, const int &order
) const {
  return 2.0*D0CIntFO(s0, weight, sTau, astau, order);
}

double TheoreticalMoments::cIntVpAD0CI(
  const double &s0, const Weight &weight, const double &sTau,
  const double &astau, const int &order
) const {
  return 2.*D0CIntCI(s0, weight, sTau, astau, order);
}

double TheoreticalMoments::cIntVpAD4FO(
  const double &s0, const Weight &weight, const double &sTau,
  const double &astau, const double &aGGinv
) const {
  return  D4CInt(s0, weight, sTau, astau, aGGinv, 1) + D4CInt(s0, weight, sTau, astau, aGGinv, -1);
}

double TheoreticalMoments::del0(
  const double &s0, const Weight &weight,
  const double &sTau, const double &astau, const int &order
) const {
  return (cIntVpAD0FO(s0, weight, sTau, astau, order)
          - cIntVpAD0FO(s0, weight, sTau, astau, 0)
  )/3.0;
}
double TheoreticalMoments::del4(
  const double &s0, const Weight &weight, const double &sTau,
  const double &astau, const double &aGGinv
) const {
  return cIntVpAD4FO(s0, weight, sTau, astau, aGGinv)/3.;
}
double TheoreticalMoments::del6(
  const double &s0, const Weight &weight,
  const double &rhoVpA
) const {
  return D68CInt(s0, weight, rhoVpA, 0.0)/3.;
}
double TheoreticalMoments::del8(
  const double &s0, const Weight &weight,
  const double &c8VpA
) const {
  return D68CInt(s0, weight, 0.0, c8VpA)/3.;
}
double TheoreticalMoments::del68(
  const double &s0, const Weight &weight,
  const double &rhoVpA, const double &c8VpA
) const {
  return D68CInt(s0, weight, rhoVpA, c8VpA)/3.;
}
