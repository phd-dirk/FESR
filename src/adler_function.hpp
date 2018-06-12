#ifndef SRC_ADLER_FUNCTION_H
#define SRC_ADLER_FUNCTION_H

#include "./types.hpp"
#include "./numerics.hpp"
#include "./alpha_s.hpp"
#include "./mq_run.hpp"
#include "./constants.hpp"
#include "./alpha_s.hpp"
#include "./weights.hpp"
#include <complex>
#include <stdexcept>
#include <functional>
#include <cmath>
#include <iostream>
#include <vector>

using std::invalid_argument;
using std::vector;
using std::pow;
using std::log;
using std::exp;
using std::function;
using std::cout;
using std::endl;
using namespace std::complex_literals;

class AdlerFunction : Numerics {
public:
  AdlerFunction(const Constants &constants) :
    Numerics(constants), const_(constants), amu_(constants), mq_(constants) {
  }

  complex<double> D0(const complex<double> &s, const complex<double> &mu2,
                     const double &astau, const double &order) const;
  complex<double> D0CI(const complex<double> &s, const complex<double> &mu2,
                     const double &astau, const double &order) const;
  // Contour integral for D_V+A(D=0) in FOPT
  double D0CIntFO(const double &s0, const Weight weight, const double &astau, const double &order) const;
  double D0CIntCI(const double &s0, const Weight weight, const double &astau, const double &order) const;

  complex<double> D2(const complex<double> &s,
                     const complex<double> mu2, const double &astau,
                     const int &order, const int &r) const;
  double D2CInt(const double &s0, const Weight &weight, const double &astau, const int &order, const int &r) const {
    function<complex<double>(complex<double>)> f =
      [&](complex<double> s) -> complex<double> {
      complex<double> mu2 = s0;
      return weight.wD(s)*D2(s0*s, mu2, astau, order, r);
    };

    return (3*const_.kPi*complexContourIntegral(f)).real();
  };

  complex<double> D4(const complex<double> &s, const complex<double> &mu2,
                     const double &astau, const double &aGGinv,
                     const int &order, const int &r) const;
  double D4CInt(const double &s0, const Weight &weight, const double &astau,
                const double &aGGinv, const int &r) const;

  complex<double> D68(const complex<double> &s, const double &rhoVpA, const double &c8VpA) const {
    return 3.e-2*rhoVpA/pow(s, 3) + 4.e-2*c8VpA/pow(s, 4);
  }
  double D68CInt(const double &s0, const Weight &weight, const double &rhoVpA, const double &c8VpA) const {
    function<complex<double>(complex<double>)> f =
      [&](complex<double> s) -> complex<double> {
      return weight.wD(s)*D68(s0*s, rhoVpA, c8VpA);
    };

    return (3*const_.kPi*complexContourIntegral(f)).real();
  };

  // Pseudoscalar contribution from pion pion pole and excited resonances
  double deltaP(const double &s0, const Weight &weight) const {
    double sTau = pow(const_.kMTau, 2);
    double spi = pow(const_.kPionMinusMass, 2);
    double pionPole = -4.*pow(const_.kFPi, 2)/s0*spi/(sTau + 2.*spi)
      *weight.wR(spi/s0).real();
    double xth = 9.*spi/s0;

    function<double(double)> f = [&](const double &s) -> double {
      double x = s0*s;
      double rhores =  2./pow(x, 2)
      *(pow(const_.kF1P, 2)*pow(const_.kM1P, 4)*breitwigner(x, const_.kM1P, const_.kG1P)
        + pow(const_.kF2P, 2)*pow(const_.kM2P, 4)*breitwigner(x, const_.kM2P, const_.kG2P) );

      return weight.wR(s).real()*2.*x/(sTau + 2.*x)*rhores;
    };

    return 4.*pow(const_.kPi, 2)*( pionPole - adaptiveIntegrate(f, xth, 1.));
  }

  double breitwigner(const double &s, const double &mbw, const double &gbw) const {
    return mbw*gbw/const_.kPi/(pow(s - pow(mbw, 2), 2)+ pow(mbw, 2)*pow(gbw, 2));
  }

 private:
  Constants const_;
  AlphaS amu_;
  MQRun mq_;
}; // end AdlerFunction

#endif
