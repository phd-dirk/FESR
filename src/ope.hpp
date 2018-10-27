#ifndef SRC_OPE_H
#define SRC_OPE_H

#include "./configuration.hpp"
#include "./numerics.hpp"
#include "./alpha_s.hpp"
#include "./mq_run.hpp"
#include "./alpha_s.hpp"
#include "./weights.hpp"
#include "./condensates.hpp"
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <exception>
#include <boost/numeric/ublas/matrix.hpp>

using boost::numeric::ublas::matrix;

class OPE : Numerics {
public:
  static complex<double> D0(
    const complex<double> &s, const cmplx &mu2, const double &sTau,
    const double &astau, const matrix<double> &c, const double &order
  );

  // Contour integral for D_V+A(D=0) in FOPT
  static double D0CIntFO(
    const double &s0, const Weight weight, const double &sTau,
    const double &astau, const matrix<double> &c, const double &order
  );
  static double D0CIntCI(
    const double &s0, const Weight weight, const double &sTau, const double &astau,
    const matrix<double> &c, const double &order
  );

  static complex<double> D4(
    const complex<double> &s, const complex<double> &mu2, const double &sTau,
    const double &astau, const double &aGGinv, const int &order, const int &r,
    const std::vector<double> &mq, const Condensates &condensates
  );
  static double D4CInt(
    const double &s0, const Weight &weight, const double &sTau,
    const double &astau, const double &aGGinv, const int &r,
    const std::vector<double> &mq, const Condensates &condensates
  );
  static complex<double> D68(
    const complex<double> &s, const double &rho, const double &c8
  );
  static double D68CInt(
    const double &s0, const Weight &weight, const double &rho, const double &c8
  );
};

#endif
