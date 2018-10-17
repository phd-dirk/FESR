#ifndef SRC_THEORETICAL_MOMENTS_H
#define SRC_THEORETICAL_MOMENTS_H

#include "./configuration.hpp"
#include "./adler_function.hpp"
#include "./duality_violations.hpp"

class TheoreticalMoments:
    public AdlerFunction, DualityViolations
{
 public:
  TheoreticalMoments(const Configuration &config);
  TheoreticalMoments();

  double thMom(
    cDbl &s0, const Weight &w, cDbl &astau,
    cDbl &aGGinv, cDbl &rhoVpA, cDbl &c8VpA, cDbl &order,
    cDbl &deV, cDbl &gaV, cDbl &alV, cDbl &beV,
    cDbl &deA, cDbl &gaA, cDbl &alA, cDbl &beA
  ) const;

  double cIntVpAD0FO(
    const double &s0, const Weight &weight, const double &sTau,
    const double &astau, const int &order) const;
  double cIntVpAD0CI(
    const double &s0, const Weight &weight, const double &sTau,
    const double &astau, const int &order) const;
  double cIntVpAD4FO(
    const double &s0, const Weight &weight, const double &sTau,
    const double &astau, const double &aGGinv
  ) const;

  double del0(
    const double &s0, const Weight &weight,
    const double &sTau, const double &astau, const int &order
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
  double vud_, SEW_;

};

#endif
