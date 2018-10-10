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

  return pow(config_.kVud, 2)*config_.kSEW*rTauTh;
}
