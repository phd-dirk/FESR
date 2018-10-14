#include "./adler_function.hpp"

typedef std::function<complex<double>(complex<double>)> cmplxFunc;

complex<double> AdlerFunction::D0(
  const complex<double> &s, const complex<double> &mu2, const double &sTau,
  const double &astau, const double &order
) const {
  complex<double> L = log(-s/mu2);
  complex<double> amu = amu_(mu2, sTau, astau/M_PI);

  complex<double> sum(0., 0.);
  for (int n = 1; n <= order; n++) {
    for (int k = 1; k <= n; k++) {
      // cout << "c(" << n << ", " << k << ") \t" << config_.c[n][k] << endl;
      complex<double> powL = ( k-1 == 0 ) ? 1.0 : pow(L, k-1);
      sum += pow(amu, n)*(double)k*config_.c[n][k]*powL;
    }
  }

  return 1/4./pow(M_PI, 2)*(config_.c[0][1] + sum);
}

double AdlerFunction::D0CIntFO(
  const double &s0, const Weight weight,
  const double &sTau, const double &astau,
  const double &order
) const {
   cmplxFunc f =
     [&](cmplx s) -> cmplx {
       cmplx mu2(s0, 0.);
       return weight.wD(s)*D0(s0*s, mu2, sTau, astau, order);
     };
  return (3*M_PI*complexContourIntegral(f)).real();
};

double AdlerFunction::D0CIntCI(
  const double &s0, const Weight weight, const double &sTau, const double &astau,
  const double &order
) const {
  cmplxFunc f =
    [&](complex<double> x) -> complex<double> {
    return weight.wD(x)*D0(s0*x, -x*s0, sTau, astau, order);
  };

  return (3*M_PI*gaussIntegration(f)).real();
};

cmplx AdlerFunction::D4(
  const cmplx &s, const cmplx &mu2, const double &sTau, const double &astau,
  const double &aGGinv, const int &order, const int &r
) const {
  const int i = 0, j = 1;
  cmplx L = log(-s/mu2);
  cmplx amu = amu_(mu2, sTau, astau/M_PI);

  cmplx gluonCondensate(0.0, 0.0);
  const double pLT3 = 0; // the coefficient is not yet known

  if (order > 0)
    gluonCondensate += 1./6.;
  if (order > 1)
    gluonCondensate += -11./108.*amu;
  if (order > 2)
    gluonCondensate += (11./48. * log(-s/mu2) + pLT3/6. - 8773./15552.)*pow(amu, 2);

  gluonCondensate *= aGGinv/pow(s, 2);

  complex<double> quarkCondensate(0.0, 0.0);
  const double qLT3 = 0., rLT3 = 0., tLT3 = 0.;

  const double mqqa = config_.mq[i]*config_.qqinv[i] + config_.mq[j]*config_.qqinv[j];
  const double mqqb = r*(config_.mq[i]*config_.qqinv[j] + config_.mq[j]*config_.qqinv[i]);
  const double mqqs = config_.mq[0]*config_.qqinv[0] + config_.mq[1]*config_.qqinv[1] + config_.mq[2]*config_.qqinv[2];

  if (order > -1)
    quarkCondensate += 2.*mqqa;
  if (order > 0)
    quarkCondensate += (-2.*mqqa + 8./3.*mqqb + 8./27.*mqqs)*amu;
  if (order > 1)
    quarkCondensate += (
      (4.5*L - 131./12.)*mqqa
      + (-6.*L + 68./3.)*mqqb
      + (-2./3.*L - 176./243. + 8./3.*config_.zeta[3])*mqqs
    )*pow(amu, 2);
  if (order > 2)
    quarkCondensate += (
      (-81./8.*pow(L, 2) + 457./8.*L + 2.*qLT3)*mqqa
      + (27./2.*pow(L, 2) - 338/3.*L + 8./3.*tLT3)*mqqb
      + (
        1.5*pow(L, 2) + (56./27. -12.*config_.zeta[3])*L
        + 50407./17496. + 8./27.*pLT3 + rLT3 -20./27.*config_.zeta[3]
      )*mqqs
    )*pow(amu, 3);
  quarkCondensate /= pow(s, 2);

  // m4
  cmplx m4(0.0, 0.0);
  cmplx rmq = mq_(mu2, config_.sTau_, astau/M_PI);

  cmplx m4a = (pow(config_.mq[i], 4) + pow(config_.mq[j], 4))*pow(rmq, 4);
  cmplx m4b = r*(config_.mq[i]*pow(config_.mq[j], 3) + config_.mq[j]*pow(config_.mq[i], 3))*pow(rmq, 4);
  cmplx m4c = (pow(config_.mq[i], 2)*pow(config_.mq[j], 2))*pow(rmq, 4);
  cmplx m4d = (pow(config_.mq[0], 4) + pow(config_.mq[1], 4) + pow(config_.mq[2], 4))*pow(rmq, 4);

  if ( order > 1)
    m4 += (-6./7./amu + 1.5*L - 0.25)*m4a - 8./7.*m4b + 3.*m4c - m4d/14.;
  if ( order > 2)
    m4 += ((-3.*pow(L, 2) + 74./7.*L)*m4a + 32./7.*L*m4b - 12.*L*m4c + 2./7.*L*m4d)*amu;

  m4 /= pow(M_PI*s, 2);

  return gluonCondensate + quarkCondensate + m4;
}

double AdlerFunction::D4CInt(
  const double &s0, const Weight &weight, const double &sTau,
  const double &astau, const double &aGGinv, const int &r
) const {
  cmplxFunc fTest =
    [&](cmplx s) -> cmplx const {
    return weight.wD(s)/pow(s, 2);
  };

  double norder = 2;
  if ( abs(complexContourIntegral(fTest).real()) < 2.e-14) {
    norder = 3;
  }

  cmplxFunc f =
    [&](cmplx s) -> cmplx {
    cmplx mu2(s0, 0.);
    return weight.wD(s)*D4(s0*s, mu2, sTau, astau, aGGinv, norder, r);
  };

  return (3*M_PI*complexContourIntegral(f)).real();
};

complex<double> AdlerFunction::D68(
  const complex<double> &s, const double &rho, const double &c8
) const {
  return 0.03*rho/pow(s, 3) + 0.04*c8/pow(s, 4);
}

double AdlerFunction::D68CInt(
  const double &s0, const Weight &weight, const double &rho, const double &c8
) const {
  cmplxFunc f =
    [&](cmplx s) -> cmplx {
      return weight.wD(s)*D68(s0*s, rho, c8);
    };

  return (3*M_PI*complexContourIntegral(f)).real();
};

double AdlerFunction::deltaP(const double &s0, const Weight &weight) const {
  double spi = pow(config_.kPionMinusMass, 2);
  double pionPole = -4.*pow(config_.fPi_, 2)/s0*spi/(config_.sTau_ + 2.*spi)
    *weight.wR(spi/s0).real();
  double xth = 9.*spi/s0;

  function<double(double)> f =
    [&](const double &s) -> double {
      double x = s0*s;
      double rhores =  2./pow(x, 2)
        *(
          pow(config_.kF1P, 2)*pow(config_.kM1P, 4)*breitwigner(x, config_.kM1P, config_.kG1P)
          + pow(config_.kF2P, 2)*pow(config_.kM2P, 4)*breitwigner(x, config_.kM2P, config_.kG2P)
        );
      return weight.wR(s).real()*2.*x/(config_.sTau_ + 2.*x)*rhores;
    };

  return 4.*pow(M_PI, 2)*( pionPole - adaptiveIntegrate(f, xth, 1.));
}

double AdlerFunction::breitwigner(
  const double &s, const double &mbw, const double &gbw
) const {
  return mbw*gbw/M_PI/(pow(s - pow(mbw, 2), 2)+ pow(mbw, 2)*pow(gbw, 2));
}
