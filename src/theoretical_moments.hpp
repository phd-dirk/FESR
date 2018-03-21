#ifndef SRC_THEORETICAL_MOMENTS_H
#define SRC_THEORETICAL_MOMENTS_H

#include "./constants.hpp"
#include "./numerics.hpp"
#include "./alpha_s.hpp"
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
  AdlerFunction(const int &nc, const int &nf, const int &order) : Numerics(1e-10, 0), nc_(nc), nf_(nf), order_(order),
                                             beta_(5), c_(6, vector<double>(6)) {
    if (order > 5) { throw invalid_argument("order cannot be higher than 5"); };

    beta_[1] = 1./6.*(11.*nc_ - 2.*nf_);
    beta_[2] = 51./4. - 19./12.*nf_;  // rgm06
    beta_[3] = 2857./64. - 5033./576.*nf_ + 325./1728.*pow(nf_, 2);  // rgm06
    beta_[4] = 149753./768. + 891./32.*zeta_[3]  // rgm06
        -(1078361./20736. + 1627./864.*zeta_[3])*nf_
        + (50065./20736. + 809./1296.*zeta_[3])*pow(nf_, 2)
        + 1093./93312.*pow(nf_, 3);

    c_[0][0] = -5./3.; c_[0][1] = 1;  // rgm06
    c_[1][1] = 1.; c_[1][2] = 0.;  //  rgm06
    // c_[2][1] = 1.63982120489698476474;  // Matthias cVA0_const.f90
    c_[2][1] = 365./24. - 11.*zeta_[3] - (11./12. - 2./3.*zeta_[3])*nf_;
    c_[2][2] = -beta_[1]*c_[1][1]/4.; c_[2][3] = 0.;  // rgm06
    // c_[3][1] = 6.37101448310094071138;  // Matthias cVA0_const.f90
    c_[3][1] = 87029./288. - 1103./4.*zeta_[3] + 275./6.*zeta_[5]
        +(-7847./216. + 262./9.*zeta_[3] - 25./9.*zeta_[5])*nf_
        + (151./162.-19./27.*zeta_[3])*pow(nf_, 2);  // rgm06
    c_[3][2] = -1./4.*(beta_[2]*c_[1][1]+2*beta_[1]*c_[2][1]);
    c_[3][3] = pow(beta_[1], 2)/12.*c_[1][1]; c_[3][4] = 0.;  // rgm06
    c_[4][1] = 78631453./20736. - 1704247./432.*zeta_[3]
        + 4185./8.*pow(zeta_[3], 2) + 34165./96.*zeta_[5]
        - 1995./16.*zeta_[7];  // Diogo PHD
    c_[4][2] = -1./4.*(beta_[3]*c_[1][1]+2*beta_[2]*c_[2][1]+3*beta_[1]*c_[3][1]);
    c_[4][3] = beta_[1]/24.*(5.*beta_[2]*c_[1][1]+6*beta_[1]*c_[2][1]);  // rgm-6
    c_[4][4] = -pow(beta_[1], 3)/32.*c_[1][1]; c_[4][5] = 0.;  // rgm06
    c_[5][1] = 283.;
    c_[5][2] = 1./4.*(-beta_[4]*c_[1][1] - 2.*beta_[3]*c_[2][1]-3.
                      *beta_[2]*c_[3][1]-4.*beta_[1]*c_[4][1]);
    c_[5][3] = 1./24.*(12.*c_[3][1]*pow(beta_[1], 2)+6.*beta_[1]*beta_[3]*c_[1][1]
                       +14.*beta_[2]*beta_[1]*c_[2][1]+3.*pow(beta_[2], 2)*c_[1][1]);
    c_[5][4] = 1./96.*(-12*pow(beta_[1], 3)*c_[2][1]
                       -13.*beta_[2]*pow(beta_[1], 2)*c_[1][1]);
    c_[5][5] = 1./80.*pow(beta_[1], 4)*c_[1][1];beta_[1] = 11./2. - 1./3.*nf_;  // rgm06
  }

  complex<double> D0(const complex<double> &s, const complex<double> &mu2, double astau) {
    // ATTENTION: alphaMu(mu)  is only equal to Matthias zarg() within a certain range around mu^2 ~ 3.
    complex<double> L = log(-s/mu2);
    complex<double> amu = alpha_s(sqrt(mu2), astau);
    complex<double> sum(0., 0.);
    for (int n = 1; n <= order_; n++) {
      for (int k = 1; k <= n; k++) {
        sum += pow(amu, n)*(double)k*c_[n][k]*pow(L,k-1);
      }
    }

    return 1/4./pow(kPi, 2)*(c_[0][1] + sum);
  }
  double D0CInt(const double &s0, function<complex<double>(complex<double>)> weight, const double &astau) {
    function<complex<double>(complex<double>)> f =
      [&](complex<double> s) -> complex<double> {
      complex<double> mu2 = s0;
      return weight(s)*D0(s0*s, mu2, astau);
    };

    return (3*kPi*complexContourIntegral(s0, f)).real();
  };

  complex<double> D2(const complex<double> &s, const complex<double> mu2, const double &astau,
                     const int &i, const int &j, const int &r) {

    double ePLT3 = 1.e2; // guesstimate

    complex<double> L = log(-s/mu2);
    complex<double> amu = alpha_s(s, astau);
    complex<double> rmq = runMassRatio(mu2, kSTau, astau);

    complex<double> m2a = (pow(mq[i], 2) + pow(mq[j], 2))*pow(rmq, 2);
    complex<double> m2b = r*mq[i]*mq[j]*pow(rmq, 2);
    complex<double> m2c = (pow(mq[1], 2) + pow(mq[2], 2) + pow(mq[3], 2))*pow(rmq, 2);

    complex<double> sum = 0.;

    if ( order_ > -1 )
      sum += m2a;
    if ( order_ > 0 )
      sum += ((13./3.-2.*L)*m2a+2./3.*m2b)*amu;
    if ( order_ > 1 )
      sum += ((4.25*pow(L, 2) - 26.*L + 23077./432. + 179./54.*zeta_[3] - 520./27.*zeta_[5])*m2a
              + (-17./6.*L + 769./54. - 55./27.*zeta_[3] - 5./27.*zeta_[5])*m2b
              + (-32./9. + 8./3.*zeta_[3])*m2c
              )*pow(amu, 2);
    if ( order_ > 2 )
      sum += ((-221./24.*pow(L, 3) + 1153./12.*pow(L, 2) + (-46253/108. - 1787./108.*zeta_[3] + 3380./27.*zeta_[5])*L
               + 3909929./5184. - pow(kPi, 4)/36. - 1541./648.*zeta_[3] + 26.5*pow(zeta_[3], 2) - 54265./108.*zeta_[5]
               + 79835./648.*zeta_[7])*m2a
              + (221./24.*pow(L, 2) + (-10831./108. + 715./54.*zeta_[3] + 65./53.*zeta_[5])*L + 4421./54. + ePLT3
                 - 715./54.*zeta_[3] - 65./54.*zeta_[5])*m2b
              + ((208./9. - 52./3.*zeta_[3])*L - 2222./27. + 1592./27.*zeta_[3] + 4.*pow(zeta_[3], 2) = 80./27.*zeta_[5])*m2c
              )*pow(amu, 4);

    return 3.*sum/(4*pow(kPi, 2)*s);
  }

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

      double mqqa = mq[i]*qqinv[i] + mq[j]*qqinv[j];
      double mqqb = r*(mq[i]*qqinv[j] + mq[j]*qqinv[i]);
      double mqqs = mq[0]*qqinv[0] + mq[1]*qqinv[1] + mq[2]*qqinv[2];

      complex<double> sum(0., 0.);
      if (order_ > -1)
        sum += 2.*mqqa;
      if (order_ > 0)
        sum += (-2.*mqqa + 8./3.*mqqb + 8./27.*mqqs)*amu;
      if (order_ > 1)
        sum += ((4.5*L - 131./12.)*mqqa
                + (-6.*L + 68./3.)*mqqb
                + (-2./3.*L - 176./243. + 8./3.*zeta_[3])*mqqs
                )*pow(amu, 2);
      if (order_ > 2)
        sum += ((-81./8.*pow(L, 2) + 457./8.*L + 2.*qLT3)*mqqa
                + (27./2.*pow(L, 2) - 338/3.*L + 8./3.*tLT3)*mqqb
                + (1.5*pow(L, 2) + (56./27. -12.*zeta_[3])*L
                   + 50407./17496. + 8./27.*pLT3 + rLT3 -20./27.*zeta_[3])*mqqs
                )*pow(amu, 3);
      return sum/pow(s, 2);
    };

    auto m4 = [&]() {
      complex<double> rmq = runMassRatio(mu2, kSTau, astau);

      complex<double> m4a = (pow(mq[i], 4) + pow(mq[j], 4))*pow(rmq, 4);
      complex<double> m4b = r*(mq[i]*pow(mq[j], 3) + mq[j]*pow(mq[i], 3))*pow(rmq, 4);
      complex<double> m4c = (pow(mq[i], 2)*pow(mq[j], 2))*pow(rmq, 4);
      complex<double> m4d = (pow(mq[1], 4) + pow(mq[2], 4) + pow(mq[3], 4))*pow(rmq, 4);

      complex<double> sum(0., 0.);

      if ( order_ > 1)
        sum += (-6./7./amu + 1.5*L - 0.25)*m4a - 8./7.*m4b + 3.*m4c - m4d/14.;
      if ( order_ > 2)
        sum += ((-3.*pow(L, 2) + 74./7.*L)*m4a + 32./7.*L*m4b - 12.*L*m4c + 2./7.*L*m4d)*amu;

      return sum/pow(kPi*s, 2);
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

    return (3*kPi*complexContourIntegral(s0, f)).real();
  };

  complex<double> D68(const complex<double> &s, const double &rhoVpA, const double &c8VpA) {
    return 3.e-2*rhoVpA/pow(s, 3) + 4.e-2*c8VpA/pow(s, 4);
  }
  double D68CInt(double s0, function<complex<double>(complex<double>)> weight,
                 const double &rhoVpA, const double &c8VpA) {
    function<complex<double>(complex<double>)> f =
      [&](complex<double> s) -> complex<double> {
      return weight(s)*D68(s0*s, rhoVpA, c8VpA);
    };

    return (3*kPi*complexContourIntegral(s0, f)).real();
  };

  double getC(int i, int j) {
    return c_[i][j];
  }

  vector<vector<double>> c_;
  vector<double> beta_;
private:
  int nc_;
  int nf_;
  int order_;
}; // end AdlerFunction

class TheoreticalMoments: public AdlerFunction {
 public:
  TheoreticalMoments(const int &nc, const int &nf, const int &order,
                     const vector<double> &s0s, function<complex<double>(complex<double>)> weight) :
      AdlerFunction(nc, nf, order), s0s(s0s), weight(weight) {}

  double operator ()(const int &i, const double &astau, const double &aGGinv) {
    // factor 2 for V PLUS A
    double s0 = s0s[i];
    // d0 VpA
    double rTauTh = 2.*D0CInt(s0, weight, astau)
      // d4 VpA
      + D4CInt(s0, weight, astau, aGGinv, 0, 1, 1 )
      + D4CInt(s0, weight, astau, aGGinv, 0, 1, -1 );
    return pow(kVud, 2)*kSEW*rTauTh;
  }

 private:
  vector<double> s0s;
  function<complex<double>(complex<double>)> weight;
};

#endif
