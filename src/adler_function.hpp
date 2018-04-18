#ifndef SRC_ADLER_FUNCTION_H
#define SRC_ADLER_FUNCTION_H

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
using std::complex;
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
  AdlerFunction(const int &order, const Constants &constants) :
    Numerics(constants), const_(constants), order_(order), amu_(constants, order), mq_(constants, order) {
    if (order > 5) { throw invalid_argument("order cannot be higher than 5");
    };
  }

  complex<double> D0(const complex<double> &s, const complex<double> &mu2,
                     const double &astau, const double &order);
  double D0CInt(const double &s0, const Weight weight, const double &astau, const double &order);

  complex<double> D2(const complex<double> &s,
                     const complex<double> mu2, const double &astau,
                     const int &order, const int &r);
  double D2CInt(const double &s0, const Weight &weight, const double &astau, const int &order, const int &r) {
    function<complex<double>(complex<double>)> f =
      [&](complex<double> s) -> complex<double> {
      complex<double> mu2 = s0;
      return weight.wD(s)*D2(s0*s, mu2, astau, order, r);
    };

    return (3*const_.kPi*complexContourIntegral(f)).real();
  };

  complex<double> D4(const complex<double> &s, const complex<double> &mu2,
                     const double &astau, const double &aGGinv,
                     const int &order, const int &r);
  double D4CInt(const double &s0, const Weight &weight, const double &astau,
                const double &aGGinv, const int &order, const int &r);

  complex<double> D68(const complex<double> &s, const double &rhoVpA, const double &c8VpA) {
    return 3.e-2*rhoVpA/pow(s, 3) + 4.e-2*c8VpA/pow(s, 4);
  }
  double D68CInt(const double &s0, const Weight &weight, const double &rhoVpA, const double &c8VpA) {
    function<complex<double>(complex<double>)> f =
      [&](complex<double> s) -> complex<double> {
      return weight.wD(s)*D68(s0*s, rhoVpA, c8VpA);
    };

    return (3*const_.kPi*complexContourIntegral(f)).real();
  };

  // Pseudoscalar contribution from pion pion pole and excited resonances
  double deltaP(const double &s0, const Weight &weight) {
    double spi = pow(const_.kPionMinusMass, 2);
    double pionPole = -4.*pow(const_.kFPi, 2)/s0*spi/(const_.kSTau + 2.*spi)
      *weight.wR(spi/s0).real();
    double xth = 9.*spi/s0;

    function<double(double)> f = [&](const double &s) -> double {
      double x = s0*s;
      double rhores =  2./pow(x, 2)
      *(pow(const_.kF1P, 2)*pow(const_.kM1P, 4)*breitwigner(x, const_.kM1P, const_.kG1P)
        + pow(const_.kF2P, 2)*pow(const_.kM2P, 4)*breitwigner(x, const_.kM2P, const_.kG2P) );

      return weight.wR(s).real()*2.*x/(const_.kSTau + 2.*x)*rhores;
    };

    return 4.*pow(const_.kPi, 2)*( pionPole - adaptiveIntegrate(f, xth, 1.));
  }

  double breitwigner(const double &s, const double &mbw, const double &gbw) {
    return mbw*gbw/const_.kPi/(pow(s - pow(mbw, 2), 2)+ pow(mbw, 2)*pow(gbw, 2));
  }

 private:
  Constants const_;
  int order_;
  AlphaS amu_;
  MQRun mq_;
}; // end AdlerFunction

#endif
