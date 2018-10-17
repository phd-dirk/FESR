#ifndef SRC_ADLER_FUNCTION_H
#define SRC_ADLER_FUNCTION_H

#include "./configuration.hpp"
#include "./numerics.hpp"
#include "./alpha_s.hpp"
#include "./mq_run.hpp"
#include "./alpha_s.hpp"
#include "./weights.hpp"
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <exception>
#include <boost/numeric/ublas/matrix.hpp>

using boost::numeric::ublas::matrix;

class AdlerFunction : Numerics {
public:
  AdlerFunction(const Configuration &config);

  complex<double> D0(
    const complex<double> &s, const cmplx &mu2, const double &sTau,
    const double &astau, const double &order
  ) const;

  // Contour integral for D_V+A(D=0) in FOPT
  double D0CIntFO(
    const double &s0, const Weight weight, const double &sTau,
    const double &astau, const double &order
  ) const;
  double D0CIntCI(
    const double &s0, const Weight weight, const double &sTau,
    const double &astau, const double &order
  ) const;

  complex<double> D4(
    const cmplx &s, const cmplx &mu2, const double &sTau,
    const double &astau, const double &aGGinv,
    const int &order, const int &r
  ) const;
  double D4CInt(
    const double &s0, const Weight &weight, const double &sTau,
    const double &astau, const double &aGGinv, const int &r
  ) const;
  complex<double> D68(
    const complex<double> &s, const double &rho, const double &c8
  ) const;
  double D68CInt(
    const double &s0, const Weight &weight, const double &rho, const double &c8
  ) const;

  // Pseudoscalar contribution from pion pion pole and excited resonances
  double deltaP(const double &s0, const Weight &weight) const;
  double breitwigner(const double &s, const double &mbw, const double &gbw) const;

 private:
  AlphaS amuRun_;
  MQRun mqRun_;
  matrix<double> c_;
  std::vector<double> mq_;
  std::vector<double> qqInv_;
  double sTau_;
  double pionMinusMass_;
  double fPi_, f1P_, m1P_, g1P_, f2P_, m2P_, g2P_;
}; // end AdlerFunction

#endif
