#include "./theoretical_moments.hpp"

ThMoms::ThMoms(const Configuration &config)
  : OPE(config), inputs_(config.inputs_)
{
  thMomContribs_ = config.thMomContribs_;
  vud_ = config.vud_;
  SEW_ = config.SEW_;
  mq_ = config.mq_;
  condensates_ = config.condensates_;
}
ThMoms::ThMoms(
  const int &nc,
  const int &nf,
  const std::vector<double> &mq,
  const Condensates &condensates,
  const double &sTau,
  const double &pionMinusMass,
  const double &fPi,
  const double &f1P,
  const double &m1P,
  const double &g1P,
  const double &f2P,
  const double &m2P,
  const double &g2P,
  const std::vector<Input> &inputs,
  const ThMomContribs &thMomContribs,
  const double &vud,
  const double &SEW
) : OPE
    (
      nc,
      nf,
      mq,
      condensates,
      sTau,
      pionMinusMass,
      fPi,
      f1P,
      m1P,
      g1P,
      f2P,
      m2P,
      g2P
    ), inputs_(inputs)
{
  thMomContribs_ = thMomContribs;
  vud_ = vud;
  SEW_ = SEW;
  nc_ = nc;
  nf_ = nf;
  mq_ = mq;
  condensates_ = condensates;
}

double ThMoms::operator() (
  cDbl &s0, const Weight &w,
  cDbl &astau, cDbl &aGGinv, cDbl &rhoVpA, cDbl &c8VpA, cDbl &order,
  cDbl &deV, cDbl &gaV, cDbl &alV, cDbl &beV,
  cDbl &deA, cDbl &gaA, cDbl &alA, cDbl &beA
) const
{
  double rTauTh = 0.;
  matrix<double> c = Configuration::adlerCoefficients(
    nf_, Configuration::betaCoefficients(nc_, nf_)
  );

  // D0
  if ( thMomContribs_.D0 ) {
    // check if FOPT or CIPT
    if ( thMomContribs_.scheme == "FO" ) {
      rTauTh += cIntVpAD0FO(s0, w, sTau_, astau, c, order);
    }
    // if ( thMomContribs_.scheme == "CI") {
    //   rTauTh += cIntVpAD0CI(s0, w, sTau_, astau, order);
    // }
  }
  // // D4
  // if ( thMomContribs_.D4 )
  //   rTauTh += cIntVpAD4FO(s0, w, sTau_, astau, aGGinv);
  // // D68
  // if ( thMomContribs_.D68 )
  //   rTauTh += D68CInt(s0, w, rhoVpA, c8VpA);
  // // DV
  // if ( thMomContribs_.DV )
  //   rTauTh += DVMomentVpA(s0, w, deV, gaV, alV, beV, deA, gaA, alA, beA);
  // // PionPole
  // if ( thMomContribs_.PionPole )
  //   rTauTh += 3.*deltaP(s0, w);

  return pow(vud_, 2)*SEW_*rTauTh;
}

double ThMoms::cIntVpAD0FO(
  const double &s0, const Weight &weight, const double &sTau,
  const double &astau, const matrix<double> &c, const int &order
) {
  return 2.0*D0CIntFO(s0, weight, sTau, astau, c, order);
}

double ThMoms::cIntVpAD0CI(
  const double &s0, const Weight &weight, const double &sTau,
  const double &astau, const int &order
) const {
  return 2.*D0CIntCI(s0, weight, sTau, astau, order);
}

double ThMoms::cIntVpAD4FO(
  const double &s0, const Weight &weight, const double &sTau,
  const double &astau, const double &aGGinv
) const {
  return  D4CInt(s0, weight, sTau, astau, aGGinv, 1, mq_, condensates_)
    + D4CInt(s0, weight, sTau, astau, aGGinv, -1, mq_, condensates_);
}

double ThMoms::del0(
  const double &s0, const Weight &weight, const double &sTau,
  const double &astau, const matrix<double> &c, const int &order
) const {
  return (
    cIntVpAD0FO(s0, weight, sTau, astau, c, order)
    - cIntVpAD0FO(s0, weight, sTau, astau, c, 0)
  )/3.0;
}
double ThMoms::del4(
  const double &s0, const Weight &weight, const double &sTau,
  const double &astau, const double &aGGinv
) const {
  return cIntVpAD4FO(s0, weight, sTau, astau, aGGinv)/3.;
}
double ThMoms::del6(
  const double &s0, const Weight &weight,
  const double &rhoVpA
) const {
  return D68CInt(s0, weight, rhoVpA, 0.0)/3.;
}
double ThMoms::del8(
  const double &s0, const Weight &weight,
  const double &c8VpA
) const {
  return D68CInt(s0, weight, 0.0, c8VpA)/3.;
}
double ThMoms::del68(
  const double &s0, const Weight &weight,
  const double &rhoVpA, const double &c8VpA
) const {
  return D68CInt(s0, weight, rhoVpA, c8VpA)/3.;
}
