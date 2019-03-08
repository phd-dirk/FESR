#include "./theoretical_moments.hpp"

ThMoms::ThMoms(const Configuration &config)
{
  thMomContribs_ = config.thMomContribs_;
  vud_ = config.vud_;
  SEW_ = config.SEW_;
  mq_ = config.mq_;
  condensates_ = config.condensates_;
}
ThMoms::ThMoms() {}

double ThMoms::calc(
  const double &s0, const Weight &w,
  const double &astau, const double &aGGinv, const double &c6,
  const double &c8, const double &c10, const double &c12,
  const double &order, const double &sTau,
  const double &deV, const double &gaV, const double &alV, const double &beV,
  const double &deA, const double &gaA, const double &alA, const double &beA,
  const double &mPiM, const double &fPi,
  const double &f1P, const double &m1P, const double &g1P,
  const double &f2P, const double &m2P, const double &g2P,
  const matrix<double> &c, const std::vector<double> &mq,
  const Condensates &condensates, const double &vud, const double &SEW,
  const ThMomContribs &thMomContribs
) {
  double rTauTh = 0.;

  // D0
  if ( thMomContribs.D0 ) {
    // check if FOPT or CIPT
    if ( thMomContribs.scheme == "FO" ) {
      rTauTh += cIntVpAD0FO(s0, w, sTau, astau, c, order);
    }
    if ( thMomContribs.scheme == "CI") {
      // rTauTh += cIntVpAD0CI(s0, w, sTau, astau, c, order);
    }
  }

  // D4
  if ( thMomContribs.D4 )
    rTauTh += cIntVpAD4(s0, w, sTau, astau, aGGinv, mq, condensates);

  // D6
  if ( thMomContribs.D6 )
    rTauTh += OPE::D6CInt(s0, w, c6);

  // D8
  if ( thMomContribs.D8 )
    rTauTh += OPE::D8CInt(s0, w, c8);

  // D10
  if ( thMomContribs.D10 )
    rTauTh += OPE::D10CInt(s0, w, c10);

  // D10
  if ( thMomContribs.D12 )
    rTauTh += OPE::D12CInt(s0, w, c12);

  // DV
  if ( thMomContribs.DV )
    rTauTh += DV::cIntVpA(s0, w, deV, gaV, alV, beV, deA, gaA, alA, beA);

  // PionPole
  if ( thMomContribs.PionPole )
    rTauTh += 3.*PSPheno::deltaP(
      s0, w, sTau, mPiM,
      fPi, f1P, m1P, g1P, f2P, m2P, g2P
    );

  return pow(vud, 2)*SEW*rTauTh;
}

double ThMoms::cIntVpAD0FO(
  const double &s0, const Weight &weight, const double &sTau,
  const double &astau, const matrix<double> &c, const int &order
) {
  return 2.0*OPE::D0CIntFO(s0, weight, sTau, astau, c, order);
}

double ThMoms::cIntVpAD0CI(
  const double &s0, const Weight &weight, const double &sTau,
  const double &astau, const matrix<double> &c, const int &order
) const {
  return 2.*OPE::D0CIntCI(s0, weight, sTau, astau, c, order);
}

double ThMoms::cIntVpAD4(
  const double &s0, const Weight &weight, const double &sTau,
  const double &astau, const double &aGGinv,
  const std::vector<double> mq, Condensates condensates
) {
  return  OPE::D4CInt(s0, weight, sTau, astau, aGGinv, 1, mq, condensates)
    + OPE::D4CInt(s0, weight, sTau, astau, aGGinv, -1, mq, condensates);
}

double ThMoms::del0(
  const double &s0, const Weight &weight, const double &sTau,
  const double &astau, const matrix<double> &c, const int &order
) {
  return (
    cIntVpAD0FO(s0, weight, sTau, astau, c, order)
    - cIntVpAD0FO(s0, weight, sTau, astau, c, 0)
  )/3.0;
}
double ThMoms::del4(
  const double &s0, const Weight &weight, const double &sTau,
  const double &astau, const double &aGGInv,
  const std::vector<double> mq, Condensates condensates
) {
  return cIntVpAD4(s0, weight, sTau, astau, aGGInv, mq, condensates)/3.0;
}
double ThMoms::del6(
  const double &s0, const Weight &weight, const double &c6
) {
  return OPE::D6CInt(s0, weight, c6)/3.0;
}
double ThMoms::del8(
  const double &s0, const Weight &weight, const double &c8
) {
  return OPE::D8CInt(s0, weight, c8)/3.0;
}
double ThMoms::del10(
  const double &s0, const Weight &weight, const double &c10
) {
  return OPE::D10CInt(s0, weight, c10)/3.0;
}
double ThMoms::del12(
  const double &s0, const Weight &weight, const double &c12
) {
  return OPE::D10CInt(s0, weight, c12)/3.0;
}

void ThMoms::logDeltas(
  const double &s0, const Weight &weight, const double &sTau,
  const double &astau, const double &aGGInv, const matrix<double> &c,
  const std::vector<double> &mq, const Condensates &condensates, const int &order,
  const double &c6, const double &c8, const double &c10, const double &c12
) {
  std::cout << "Delta contributions: " << std::endl;
  std::cout << "del^(0) \t" << ThMoms::del0(s0, weight, sTau, astau, c, order) << std::endl;
  std::cout << "del^(4) \t" << ThMoms::del4(s0, weight, sTau, astau, aGGInv, mq, condensates) << std::endl;
  std::cout << "del^(6) \t" << ThMoms::del6(s0, weight, c6) << std::endl;
  std::cout << "del^(8) \t" << ThMoms::del8(s0, weight, c8) << std::endl;
  std::cout << "del^(10) \t" << ThMoms::del10(s0, weight, c10) << std::endl;
  std::cout << "del^(12) \t" << ThMoms::del12(s0, weight, c12) << std::endl;
  std::cout << std::endl;
}
