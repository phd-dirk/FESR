#ifndef SRC_ADLER_FUNCTION_H
#define SRC_ADLER_FUNCTION_H

#include "./types.hpp"
#include "./configuration.hpp"
#include "./numerics.hpp"
#include "./alpha_s.hpp"
#include "./mq_run.hpp"
#include "./alpha_s.hpp"
#include "./weights.hpp"
#include <stdexcept>
#include <cmath>
#include <iostream>

class AdlerFunction : Numerics {
public:
  AdlerFunction(const Configuration &config) :
    Numerics(), config_(config), amu_(), mq_(config.sTau) {}

  cmplx D0(const cmplx &s, const cmplx &mu2, const double &sTau, const double &astau,
           const double &order) const;
  cmplx D0CI(const cmplx &s, const cmplx &mu2, const double &sTau, const double &astau,
            const double &order) const;
  // Contour integral for D_V+A(D=0) in FOPT
  double D0CIntFO(const double &s0, const Weight weight, const double &sTau,
                  const double &astau, const double &order) const;
  double D0CIntCI(const double &s0, const Weight weight, const double &sTau,
                  const double &astau, const double &order) const;

  cmplx D2(const cmplx &s, const cmplx mu2, const double &astau, const int &order,
           const int &r) const;
  double D2CInt(const double &s0, const Weight &weight, const double &astau, const int &order, const int &r) const {
    cmplxFunc f =
      [&](cmplx s) -> cmplx {
      cmplx mu2 = s0;
      return weight.wD(s)*D2(s0*s, mu2, astau, order, r);
    };

    return (3*M_PI*complexContourIntegral(f)).real();
  };

  cmplx D4(const cmplx &s, const cmplx &mu2, const double &sTau,
           const double &astau, const double &aGGinv,
           const int &order, const int &r) const;
  double D4CInt(const double &s0, const Weight &weight, const double &sTau,
                const double &astau, const double &aGGinv, const int &r) const;

  cmplx D68(const cmplx &s, const double &rho, const double &c8) const {
    return 0.03*rho/pow(s, 3) + 0.04*c8/pow(s, 4);
  }
  double D68CInt(const double &s0, const Weight &weight, const double &rhoVpA,
                 const double &c8VpA) const {
    cmplxFunc f =
      [&](cmplx s) -> cmplx {
      return weight.wD(s)*D68(s0*s, rhoVpA, c8VpA);
    };

    return (3*M_PI*complexContourIntegral(f)).real();
  };

  // Pseudoscalar contribution from pion pion pole and excited resonances
  double deltaP(const double &s0, const Weight &weight) const {
    double spi = pow(config_.kPionMinusMass, 2);
    double pionPole = -4.*pow(config_.kFPi, 2)/s0*spi/(config_.sTau + 2.*spi)
      *weight.wR(spi/s0).real();
    double xth = 9.*spi/s0;

    function<double(double)> f = [&](const double &s) -> double {
      double x = s0*s;
      double rhores =  2./pow(x, 2)
      *(pow(config_.kF1P, 2)*pow(config_.kM1P, 4)*breitwigner(x, config_.kM1P, config_.kG1P)
        + pow(config_.kF2P, 2)*pow(config_.kM2P, 4)*breitwigner(x, config_.kM2P, config_.kG2P) );

      return weight.wR(s).real()*2.*x/(config_.sTau + 2.*x)*rhores;
    };

    return 4.*pow(M_PI, 2)*( pionPole - adaptiveIntegrate(f, xth, 1.));
  }

  double breitwigner(const double &s, const double &mbw, const double &gbw) const {
    return mbw*gbw/M_PI/(pow(s - pow(mbw, 2), 2)+ pow(mbw, 2)*pow(gbw, 2));
  }

 private:
  Configuration config_;
  AlphaS amu_;
  MQRun mq_;
}; // end AdlerFunction

#endif
