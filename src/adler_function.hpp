#ifndef SRC_ADLER_FUNCTION_H
#define SRC_ADLER_FUNCTION_H

#include "./numerics.hpp"
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

class AdlerFunction : public Numerics {
public:
  AdlerFunction(const int &order, Constants constants) :
    Numerics(1e-13, 0, constants), const_(constants), order_(order) {
    if (order > 5) { throw invalid_argument("order cannot be higher than 5");
    };
  }

  complex<double> D0(const complex<double> &s, const complex<double> &mu2,
                     const double &astau, const double &order) {
    // ATTENTION: alphaMu(mu)  is only equal to Matthias zarg() within a certain range around mu^2 ~ 3.
    complex<double> L = log(-s/mu2);
    complex<double> amu = alpha_s(sqrt(mu2), astau);
    complex<double> sum(0., 0.);
    for (int n = 1; n <= order; n++) {
      for (int k = 1; k <= n; k++) {
        sum += pow(amu, n)*(double)k*const_.c_[n][k]*pow(L,k-1);
      }
    }

    return 1/4./pow(const_.kPi, 2)*(const_.c_[0][1] + sum);
  }
  double D0CInt(const double &s0, function<complex<double>(complex<double>)> weight,
                const double &astau, const double &order) {
    function<complex<double>(complex<double>)> f =
      [&](complex<double> s) -> complex<double> {
      complex<double> mu2 = s0;
      return weight(s)*D0(s0*s, mu2, astau, order);
    };

    return (3*const_.kPi*complexContourIntegral(s0, f)).real();
  };

  complex<double> D2(const complex<double> &s,
                     const complex<double> mu2, const double &astau,
                     const int &order, const int &r);
  double D2CInt(double s0,
                function<complex<double>(complex<double>)> weight,
                const double &astau, const int &order, const int &r) {
    function<complex<double>(complex<double>)> f =
      [&](complex<double> s) -> complex<double> {
      complex<double> mu2 = s0;
      return weight(s)*D2(s0*s, mu2, astau, order, r);
    };

    return (3*const_.kPi*complexContourIntegral(s0, f)).real();
  };

  complex<double> D4(const complex<double> &s, const complex<double> &mu2,
                     const double &astau, const double &aGGinv,
                     const int &order, const int &r);
  double D4CInt(double s0,
                function<complex<double>(complex<double>)> weight, const double &astau,
                const double &aGGinv, const int &order, const int &r);

  complex<double> D68(const complex<double> &s, const double &rhoVpA, const double &c8VpA) {
    return 3.e-2*rhoVpA/pow(s, 3) + 4.e-2*c8VpA/pow(s, 4);
  }
  double D68CInt(const double &s0, function<complex<double>(complex<double>)> weight,
                 const double &rhoVpA, const double &c8VpA) {
    function<complex<double>(complex<double>)> f =
      [&](complex<double> s) -> complex<double> {
      return weight(s)*D68(s0*s, rhoVpA, c8VpA);
    };

    return (3*const_.kPi*complexContourIntegral(s0, f)).real();
  };

  // Pseudoscalar contribution from pion pion pole and excited resonances
  double deltaP(const double &s0, function<complex<double>(complex<double>)> weightRho) {
    double spi = pow(const_.kPionMinusMass, 2);
    double pionPole = -4.*pow(const_.kFPi, 2)/s0*spi/(const_.kSTau + 2.*spi)
      *weightRho(spi/s0).real();
    double xth = 9.*spi/s0;

    function<double(double)> f = [&](const double &s) -> double {
      double x = s0*s;
      double rhores =  2./pow(x, 2)
      *(pow(const_.kF1P, 2)*pow(const_.kM1P, 4)*breitwigner(x, const_.kM1P, const_.kG1P)
        + pow(const_.kF2P, 2)*pow(const_.kM2P, 4)*breitwigner(x, const_.kM2P, const_.kG2P) );

      return weightRho(s).real()*2.*x/(const_.kSTau + 2.*x)*rhores;
    };

    return 4.*pow(const_.kPi, 2)*( pionPole - integrate(f, xth, 1.));
  }

  double breitwigner(const double &s, const double &mbw, const double &gbw) {
    return mbw*gbw/const_.kPi/(pow(s - pow(mbw, 2), 2)+ pow(mbw, 2)*pow(gbw, 2));
  }

 private:
  Constants const_;
  int order_;
}; // end AdlerFunction

#endif