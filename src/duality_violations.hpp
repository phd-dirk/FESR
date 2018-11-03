#ifndef SRC_DUALITY_VIOLATIONS_H
#define SRC_DUALITY_VIOLATIONS_H

#include "./numerics.hpp"
#include "./weights.hpp"

class DV
{
 public:
  static double cIntVpA(
    const double &s0, const Weight &w,
    const double &deV, const double &gaV, const double &alV, const double &beV,
    const double &deA, const double &gaA, const double &alA, const double &beA
  );

  static double cIntVA(
    const double &s0, const Weight &w, const double &de, const double &ga,
    const double &al, const double &be
  );

  static double rhoDV(
    const double &s, const double &delta, const double &gamma,
    const double &alpha, const double &beta
  );

  static double intP0(
    const double &s0, const double &delta, const double &gamma,
    const double &alpha, const double &beta
  );
  static double intP1(
    const double &s0, const double &delta, const double &gamma,
    const double &alpha, const double &beta
  );
  static double intP2(
    const double &s0, const double &delta, const double &gamma,
    const double &alpha, const double &beta
  );
  static double intP3(
    const double &s0, const double &delta, const double &gamma,
    const double &alpha, const double &beta
  );
};

#endif
