#ifndef SRC_PSEUDOSCALAR_PHENO_H
#define SRC_PSEUDOSCALAR_PHENO_H

#include "./weights.hpp"
#include "./numerics.hpp"

class PSPheno {
 public:
  // Pseudoscalar contribution from pion pion pole and excited resonances
  static double deltaP(
    const double &s0, const Weight &weight, const double &sTau, const double &mPiM,
    const double &fPi, const double &f1P, const double &m1P, const double &g1P,
    const double &f2P, const double &m2P, const double &g2P
  );
  static double breitwigner(const double &s, const double &mbw, const double &gbw);
};

#endif
