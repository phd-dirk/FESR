#ifndef SRC_THEORETICAL_MOMENTS_H
#define SRC_THEORETICAL_MOMENTS_H

#include "./constants.hpp"
#include "./numerics.hpp"
#include "./alpha_s.hpp"
#include "./weights.hpp"
#include <vector>
#include <functional>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <iostream>


using std::vector;
using std::pow;
using std::log;
using std::exp;
using std::invalid_argument;
using std::complex;
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

  complex<double> D2(const complex<double> &s, const complex<double> mu2, const double &astau,
                     const int &i, const int &j, const int &r) {

    // double ePLT3 = 1.e2; // guesstimate

    complex<double> L = log(-s/mu2);
    complex<double> amu = alpha_s(sqrt(mu2), astau);
    complex<double> rmq = runMassRatio(mu2, const_.kSTau, astau);

    complex<double> m2a = (pow(const_.mq[i], 2) + pow(const_.mq[j], 2))*pow(rmq, 2);
    complex<double> m2b = r*const_.mq[i]*const_.mq[j]*pow(rmq, 2);
    complex<double> m2c = (pow(const_.mq[1], 2) + pow(const_.mq[2], 2) + pow(const_.mq[3], 2))*pow(rmq, 2);

    complex<double> sum = 0.;

    if ( order_ > -1 )
      sum += m2a;
    if ( order_ > 0 )
      sum += ((13./3.-2.*L)*m2a + 2./3.*m2b)*amu;
    if ( order_ > 1 )
      sum += ((4.25*pow(L, 2) - 26.*L + 23077./432. + 179./54.*const_.zeta_[3] - 520./27.*const_.zeta_[5])*m2a
              + (-17./6.*L + 769./54. - 55./27.*const_.zeta_[3] - 5./27.*const_.zeta_[5])*m2b
              + (-32./9. + 8./3.*const_.zeta_[3])*m2c
              )*pow(amu, 2);
    // if ( order_ > 2 )
    //   sum += ((-221./24.*pow(L, 3) + 1153./12.*pow(L, 2) + (-46253/108. - 1787./108.*zeta_[3] + 3380./27.*zeta_[5])*L
    //            + 3909929./5184. - pow(kPi, 4)/36. - 1541./648.*zeta_[3] + 26.5*pow(zeta_[3], 2) - 54265./108.*zeta_[5]
    //            + 79835./648.*zeta_[7])*m2a
    //           + (221./24.*pow(L, 2) + (-10831./108. + 715./54.*zeta_[3] + 65./54.*zeta_[5])*L + 4421./54. + ePLT3
    //              - 715./54.*zeta_[3] - 65./54.*zeta_[5])*m2b
    //           + ((208./9. - 52./3.*zeta_[3])*L - 2222./27. + 1592./27.*zeta_[3] + 4.*pow(zeta_[3], 2) - 80./27.*zeta_[5])*m2c
    //           )*pow(amu, 3);

    return 3.*sum/(4*pow(const_.kPi, 2)*s);
  }
  double D2CInt(double s0, function<complex<double>(complex<double>)> weight, const double &astau,
                const int &i, const int &j, const int &r) {
    function<complex<double>(complex<double>)> f =
      [&](complex<double> s) -> complex<double> {
      complex<double> mu2 = s0;
      return weight(s)*D2(s0*s, mu2, astau, i, j, r);
    };

    return (3*const_.kPi*complexContourIntegral(s0, f)).real();
  };

  complex<double> D4(const complex<double> &s, const complex<double> &mu2,
                     const double &astau, const double &aGGinv, const int i,
                     const int &j, const int &r) {

    complex<double> L = log(-s/mu2);
    complex<double> amu = alpha_s(sqrt(mu2), astau);

    auto gluonCondensate = [&]() {
      complex<double> sum(0.0, 0.0);
      double pLT3 = 0; // the coefficient is not yet known

      if (order_ > 0)
        sum += 1./6.;
      if (order_ > 1)
        sum += -11./108.*amu;
      if (order_ > 2)
        sum += (11./48. * log(-s/mu2) + pLT3/6. - 8773./15552.)*pow(amu, 2);

      // cout << "amu : \t" << amu << endl;
      // cout << "order : \t" << order_ << endl;
      // cout << "s : \t" << s << endl;
      // cout << "mu2 : \t" << mu2 << endl;
      // cout << "aGGinv : \t" << aGGinv << endl;

      return sum*aGGinv/pow(s, 2);
    };

    auto quarkCondensate = [&]() {
      double pLT3 = 0., qLT3 = 0., rLT3 = 0., tLT3 = 0.;

      double mqqa = const_.mq[i]*const_.qqinv[i] + const_.mq[j]*const_.qqinv[j];
      double mqqb = r*(const_.mq[i]*const_.qqinv[j] + const_.mq[j]*const_.qqinv[i]);
      double mqqs = const_.mq[0]*const_.qqinv[0] + const_.mq[1]*const_.qqinv[1] + const_.mq[2]*const_.qqinv[2];

      complex<double> sum(0., 0.);
      if (order_ > -1)
        sum += 2.*mqqa;
      if (order_ > 0)
        sum += (-2.*mqqa + 8./3.*mqqb + 8./27.*mqqs)*amu;
      if (order_ > 1)
        sum += ((4.5*L - 131./12.)*mqqa
                + (-6.*L + 68./3.)*mqqb
                + (-2./3.*L - 176./243. + 8./3.*const_.zeta_[3])*mqqs
                )*pow(amu, 2);
      if (order_ > 2)
        sum += ((-81./8.*pow(L, 2) + 457./8.*L + 2.*qLT3)*mqqa
                + (27./2.*pow(L, 2) - 338/3.*L + 8./3.*tLT3)*mqqb
                + (1.5*pow(L, 2) + (56./27. -12.*const_.zeta_[3])*L
                   + 50407./17496. + 8./27.*pLT3 + rLT3 -20./27.*const_.zeta_[3])*mqqs
                )*pow(amu, 3);
      return sum/pow(s, 2);
    };

    auto m4 = [&]() {
      complex<double> rmq = runMassRatio(mu2, const_.kSTau, astau);

      complex<double> m4a = (pow(const_.mq[i], 4) + pow(const_.mq[j], 4))*pow(rmq, 4);
      complex<double> m4b = r*(const_.mq[i]*pow(const_.mq[j], 3) + const_.mq[j]*pow(const_.mq[i], 3))*pow(rmq, 4);
      complex<double> m4c = (pow(const_.mq[i], 2)*pow(const_.mq[j], 2))*pow(rmq, 4);
      complex<double> m4d = (pow(const_.mq[1], 4) + pow(const_.mq[2], 4) + pow(const_.mq[3], 4))*pow(rmq, 4);

      complex<double> sum(0., 0.);

      if ( order_ > 1)
        sum += (-6./7./amu + 1.5*L - 0.25)*m4a - 8./7.*m4b + 3.*m4c - m4d/14.;
      if ( order_ > 2)
        sum += ((-3.*pow(L, 2) + 74./7.*L)*m4a + 32./7.*L*m4b - 12.*L*m4c + 2./7.*L*m4d)*amu;

      return sum/pow(const_.kPi*s, 2);
    };

    return gluonCondensate() + quarkCondensate() + m4();
  }
  double D4CInt(double s0, function<complex<double>(complex<double>)> weight, const double &astau,
                const double &aGGinv, const int &i, const int &j, const int &r) {
    function<complex<double>(complex<double>)> f =
      [&](complex<double> s) -> complex<double> {
      complex<double> mu2 = s0;
      return weight(s)*D4(s0*s, mu2, astau, aGGinv, i, j, r);
    };

    return (3*const_.kPi*complexContourIntegral(s0, f)).real();
  };

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

class TheoreticalMoments: public AdlerFunction {
 public:
  TheoreticalMoments(const int &order, const vector<double> &s0s,
                     function<complex<double>(complex<double>)> weight, Constants constants) :
    AdlerFunction(order, constants), const_(constants), s0s(s0s), weight(weight) {}

  double operator ()(const int &i, const double &astau, const double &aGGinv,
                     const double &rhoVpA, const double &c8VpA, const double &order) {
    // factor 2 for V PLUS A
    double s0 = s0s[i];
    // d0 VpA
    double rTauTh = cIntVpAD0FO(s0, weight, astau, order)
      // d4 VpA
      + D4CInt(s0, weight, astau, aGGinv, 0, 1, 1 )
      + D4CInt(s0, weight, astau, aGGinv, 0, 1, -1 )
      // d68 VpA
      + D68CInt(s0, weight, rhoVpA, c8VpA)
      // deltaP
      + 3.*deltaP(s0, wR00);
    return pow(const_.kVud, 2)*const_.kSEW*rTauTh;
  }

  double cIntVpAD0FO(const double &s0,
                     function<complex<double>(complex<double>)> weight,
                     const double &astau, const double &order) {
    return 2.*D0CInt(s0, weight, astau, order);
  }

  double del0(const double &s0,
                   function<complex<double>(complex<double>)> weight,
                   const double &astau, const int &order) {
    return ( cIntVpAD0FO(s0, weight, astau, order) - cIntVpAD0FO(s0, weight, astau, order) ) / 3.0;
  }

  double del68(const double &s0, function<complex<double>(complex<double>)> weight,
               const double &rhoVpA, const double &c8VpA) {
    return D68CInt(s0, weight, rhoVpA, c8VpA)/3.;
  }

  double del4(const double &s0,
              function<complex<double>(complex<double>)> weight,
              const double &astau, const double &aGGinv) {
    return (D4CInt(s0, weight, astau, aGGinv, 0, 1, 1)
            + D4CInt(s0, weight, astau, aGGinv, 0, 1, -1)
            )/3.;
  }

 private:
  Constants const_;
  vector<double> s0s;
  function<complex<double>(complex<double>)> weight;
};

#endif
