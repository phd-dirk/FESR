#include "./adler_function.hpp"

cmplx AdlerFunction::D0(const cmplx &s, const cmplx &mu2, const double &sTau,
                        const double &astau, const double &order) const {
  cmplx L = log(-s/mu2);

  cmplx amu = amu_(mu2, sTau, astau/M_PI);
  // cout << "mu: " << mu << "\t" << "s: " << s << endl;
  // cout << "astau: " << astau/M_PI << endl;
  // cout << "L: " << L << "\t" << "amu: " << amu << endl;

  cmplx sum(0., 0.);
  for (int n = 1; n <= order; n++) {
    for (int k = 1; k <= n; k++) {
      // cout << "c(" << n << ", " << k << ") \t" << const_.c_[n][k] << endl;
      sum += pow(amu, n)*(double)k*config_.c[n][k]*pow(L,k-1);
      // cout << sum  << "\t" << pow(L,k-1) << endl;
    }
  }

  return 1/4./pow(M_PI, 2)*(config_.c[0][1] + sum);
}

cmplx AdlerFunction::D0CI(const cmplx &s, const cmplx &mu, const double &sTau, const double &astau, const double &order) const {
  cmplx amu = amu_(mu, sTau, astau/M_PI);

  cmplx sum(0., 0.);
  for (int n = 1; n <= order; n++) {
    int k = 1;
    sum += pow(amu, n)*(double)k*config_.c[n][k];
  }

  return 1/4./pow(M_PI, 2)*(config_.c[0][1] + sum);
}

double AdlerFunction::D0CIntFO(const double &s0, const Weight weight,
                               const double &sTau, const double &astau,
                               const double &order) const {
   cmplxFunc f =
    [&](cmplx s) -> cmplx {
    cmplx mu2(s0, 0.);

    return weight.wD(s)*D0(s0*s, mu2, sTau, astau, order);
  };

  return (3*M_PI*complexContourIntegral(f)).real();
};

double AdlerFunction::D0CIntCI(const double &s0, const Weight weight,
                               const double &sTau, const double &astau,
                               const double &order) const {
  function<complex<double>(complex<double>)> f =
    [&](complex<double> x) -> complex<double> {
    complex<double> xmu2 = -x*s0;
    return weight.wD(x)*D0CI(s0*x, xmu2, astau, sTau, order);
  };

  return (3*M_PI*complexContourIntegral(f)).real();
};

cmplx AdlerFunction::D2(const cmplx &s, const cmplx mu, const double &astau,
                        const int &order, const int &r) const {
  // double ePLT3 = 1.e2; // guesstimate
  const int i = 0, j = 1;

  cmplx mu2 = pow(mu, 2);
  cmplx L = log(-s/mu2);
  cmplx amu = amu_(mu, config_.sTau, astau/M_PI);
  cmplx rmq = mq_(mu, config_.sTau, astau/M_PI);

  cmplx m2a = (pow(config_.mq[i], 2) + pow(config_.mq[j], 2))*pow(rmq, 2);
  cmplx m2b = r*config_.mq[i]*config_.mq[j]*pow(rmq, 2);
  cmplx m2c = (pow(config_.mq[0], 2) + pow(config_.mq[1], 2) + pow(config_.mq[2], 2))*pow(rmq, 2);

  cmplx sum = 0.;

  if ( order > -1 )
    sum += m2a;
  if ( order > 0 )
    sum += ((13./3.-2.*L)*m2a + 2./3.*m2b)*amu;
  if ( order > 1 )
    sum += ((4.25*pow(L, 2) - 26.*L + 23077./432. + 179./54.*config_.zeta[3] - 520./27.*config_.zeta[5])*m2a
            + (-17./6.*L + 769./54. - 55./27.*config_.zeta[3] - 5./27.*config_.zeta[5])*m2b
            + (-32./9. + 8./3.*config_.zeta[3])*m2c
            )*pow(amu, 2);
  // if ( order > 2 )
  //   sum += ((-221./24.*pow(L, 3) + 1153./12.*pow(L, 2) + (-46253/108. - 1787./108.*zeta[3] + 3380./27.*zeta[5])*L
  //            + 3909929./5184. - pow(kPi, 4)/36. - 1541./648.*zeta[3] + 26.5*pow(zeta[3], 2) - 54265./108.*zeta[5]
  //            + 79835./648.*zeta[7])*m2a
  //           + (221./24.*pow(L, 2) + (-10831./108. + 715./54.*zeta[3] + 65./54.*zeta[5])*L + 4421./54. + ePLT3
  //              - 715./54.*zeta[3] - 65./54.*zeta[5])*m2b
  //           + ((208./9. - 52./3.*zeta[3])*L - 2222./27. + 1592./27.*zeta[3] + 4.*pow(zeta[3], 2) - 80./27.*zeta[5])*m2c
  //           )*pow(amu, 3);

  return 3.*sum/(4*pow(M_PI, 2)*s);
}

cmplx AdlerFunction::D4(const cmplx &s, const cmplx &mu,
                        const double &astau, const double &aGGinv,
                        const int &order, const int &r) const {

  const int i = 0, j = 1;
  cmplx mu2 = pow(mu, 2);
  cmplx L = log(-s/mu2);
  cmplx amu = amu_(mu, config_.sTau, astau/M_PI);


  cmplx gluonCondensate(0.0, 0.0);
  const double pLT3 = 0; // the coefficient is not yet known

  if (order > 0)
    gluonCondensate += 1./6.;
  if (order > 1)
    gluonCondensate += -11./108.*amu;
  if (order > 2)
    gluonCondensate += (11./48. * log(-s/mu2) + pLT3/6. - 8773./15552.)*pow(amu, 2);

  gluonCondensate *= aGGinv/pow(s, 2);


  cmplx quarkCondensate(0.0, 0.0);
  const double qLT3 = 0., rLT3 = 0., tLT3 = 0.;

  const double mqqa = config_.mq[i]*config_.qqinv[i] + config_.mq[j]*config_.qqinv[j];
  const double mqqb = r*(config_.mq[i]*config_.qqinv[j] + config_.mq[j]*config_.qqinv[i]);
  const double mqqs = config_.mq[0]*config_.qqinv[0] + config_.mq[1]*config_.qqinv[1] + config_.mq[2]*config_.qqinv[2];

  if (order > -1)
    quarkCondensate += 2.*mqqa;
  if (order > 0)
    quarkCondensate += (-2.*mqqa + 8./3.*mqqb + 8./27.*mqqs)*amu;
  if (order > 1)
    quarkCondensate += ((4.5*L - 131./12.)*mqqa
                        + (-6.*L + 68./3.)*mqqb
                        + (-2./3.*L - 176./243. + 8./3.*config_.zeta[3])*mqqs
                        )*pow(amu, 2);
  if (order > 2)
    quarkCondensate += ((-81./8.*pow(L, 2) + 457./8.*L + 2.*qLT3)*mqqa
                        + (27./2.*pow(L, 2) - 338/3.*L + 8./3.*tLT3)*mqqb
                        + (1.5*pow(L, 2) + (56./27. -12.*config_.zeta[3])*L
                           + 50407./17496. + 8./27.*pLT3 + rLT3 -20./27.*config_.zeta[3])*mqqs
                        )*pow(amu, 3);
  quarkCondensate /= pow(s, 2);

  cmplx m4(0.0, 0.0);
  cmplx rmq = mq_(mu, config_.sTau, astau/M_PI);

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
double AdlerFunction::D4CInt(const double &s0, const Weight &weight, const double &sTau,
                             const double &astau, const double &aGGinv, const int &r) const {
  cmplxFunc fTest =
    [&](cmplx s) -> cmplx const {
    return weight.wD(s)/pow(s, 2);
  };

  double norder = 2;
  if ( abs(gaussIntegration(fTest).real()) < 2.e-14) {
    norder = 3;
  }

  cmplxFunc f =
    [&](cmplx s) -> cmplx {
    cmplx mu2(s0, 0.);
    cmplx mu = sqrt(mu2);
    return weight.wD(s)*D4(s0*s, mu, astau, aGGinv, norder, r);
  };

  return (3*M_PI*gaussIntegration(f)).real();
};
