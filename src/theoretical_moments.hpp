#ifndef SRC_THEORETICAL_MOMENTS_H
#define SRC_THEORETICAL_MOMENTS_H

#include "./configuration.hpp"
#include "./ope.hpp"
#include "./duality_violations.hpp"
#include "./condensates.hpp"
#include "./pseudoscalar_pheno.hpp"

class ThMoms
{
 public:
  ThMoms(const Configuration &config);
  ThMoms();

  static double calc(
    const double &s0, const Weight &w,
    const double &astau, const double &aGGinv, const double &rhoVpA,
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
  );

  static double cIntVpAD0FO(
    const double &s0, const Weight &weight, const double &sTau,
    const double &astau, const matrix<double> &c, const int &order
  );
  double cIntVpAD0CI(
    const double &s0, const Weight &weight, const double &sTau,
    const double &astau, const matrix<double> &c, const int &order
  ) const;
  static double cIntVpAD4(
    const double &s0, const Weight &weight, const double &sTau,
    const double &astau, const double &aGGinv,
    const std::vector<double> mq, Condensates condensates
  );

  static double del0(
    const double &s0, const Weight &weight, const double &sTau,
    const double &astau, const matrix<double> &c, const int &order
  );
  static double del4(
    const double &s0, const Weight &weight, const double &sTau,
    const double &astau, const double &aGGinv,
    const std::vector<double> mq, Condensates condensates
  );
  static double del6(
    const double &s0, const Weight &weight, const double &c6
  );
  static double del8(
    const double &s0, const Weight &weight, const double &c8
  );
  static double del10(
    const double &s0, const Weight &weight, const double &c10
  );
  static double del12(
    const double &s0, const Weight &weight, const double &c12
  );

  static void logDeltas(
    const double &s0, const Weight &weight, const double &sTau,
    const double &astau, const double &aGGInv, const matrix<double> &c,
    const std::vector<double> &mq, const Condensates &condensates, const int &order,
    const double &c6, const double &c8, const double &c10, const double &c12
  );

 private:
  ThMomContribs thMomContribs_;
  double vud_, SEW_;
  std::vector<double> mq_;
  Condensates condensates_;
};

#endif
