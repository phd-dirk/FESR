#ifndef SRC_THEORETICAL_MOMENTS_H
#define SRC_THEORETICAL_MOMENTS_H

#include "./configuration.hpp"
#include "./ope.hpp"
#include "./duality_violations.hpp"
#include "./condensates.hpp"

class ThMoms:
    public OPE, DualityViolations
{
 public:
  ThMoms(const Configuration &config);
  ThMoms(
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
  );

  double operator() (
    cDbl &s0, const Weight &w, cDbl &astau,
    cDbl &aGGinv, cDbl &rhoVpA, cDbl &c8VpA, cDbl &order,
    cDbl &deV, cDbl &gaV, cDbl &alV, cDbl &beV,
    cDbl &deA, cDbl &gaA, cDbl &alA, cDbl &beA
  ) const;

  static double cIntVpAD0FO(
    const double &s0, const Weight &weight, const double &sTau,
    const double &astau, const matrix<double> &c, const int &order
  );
  double cIntVpAD0CI(
    const double &s0, const Weight &weight, const double &sTau,
    const double &astau, const int &order
  ) const;
  double cIntVpAD4FO(
    const double &s0, const Weight &weight, const double &sTau,
    const double &astau, const double &aGGinv
  ) const;

  double del0(
    const double &s0, const Weight &weight, const double &sTau,
    const double &astau, const matrix<double> &c, const int &order
  ) const;
  double del4(
    const double &s0, const Weight &weight, const double &sTau,
    const double &astau, const double &aGGinv
  ) const;
  double del6(
    const double &s0, const Weight &weight,
    const double &rhoVpA
  ) const;
  double del8(
    const double &s0, const Weight &weight,
    const double &c8VpA) const;
  double del68(
    const double &s0, const Weight &weight,
    const double &rhoVpA, const double &c8VpA
  ) const;

 private:
  std::vector<Input> inputs_;
  ThMomContribs thMomContribs_;
  double nc_, nf_, vud_, SEW_;
  std::vector<double> mq_;
  Condensates condensates_;
};

#endif
