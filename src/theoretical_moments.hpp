#ifndef SRC_THEORETICAL_MOMENTS_H
#define SRC_THEORETICAL_MOMENTS_H

#include "../lib/cRunDec/CRunDec.h"
#include "./constants.hpp"
#include "./numerics.hpp"
#include <vector>
#include <functional>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <iostream>


using std::vector;
using std::pow;
using std::log;
using std::invalid_argument;
using std::complex;
using std::cout;
using std::endl;
using std::bind;
using std::placeholders::_1;


class AdlerFunction : public Numerics {
public:
  AdlerFunction(int nc, int nf, int order) : Numerics(1e-7, 1e-7), nc_(nc), nf_(nf), order_(order),
                                             beta_(5), c_(6, vector<double>(6)),
                                             cRunDec_(CRunDec(nf)) {
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

  complex<double> D0(complex<double> s, double mu) {
    // ATTENTION: alphaMu(mu)  is only equal to Matthias zarg() within a certain range around mu^2 ~ 3.
    complex<double> L = log(-s/pow(mu, 2));
    complex<double> amu(alphaMu(pow(mu, 2))/kPi, 0.);
    complex<double> sum(0., 0.);
    for (int n = 1; n <= order_; n++) {
      for (int k = 1; k <= n; k++) {
        sum += pow(amu, n)*(double)k*c_[n][k]*pow(L,k-1);
      }
    }
    return 1/4./pow(kPi, 2)*(c_[0][1] + sum);
  }

  complex<double> contourIntegral(double mu) {
    auto D0Real = [mu, this](double t) {
      return D0(complex<double>(t, t), mu).real();
    };
    auto D0Imag = [mu, this](double t) {
      return D0(complex<double>(t, t), mu).imag();
    };

    double cintReal = integrate(D0Real);
    double cintImag = integrate(D0Imag);

    return complex<double>(cintReal, cintImag);
  };

  double alphaMu(double mu) {

    return cRunDec_.AlphasExact(kAlphaTau, kSTau, mu, 4);
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
  CRunDec cRunDec_;
};

CRunDec cRunDec(5);
double jo = cRunDec.AlphasExact(0.108, 1. ,2., 4);

#endif
