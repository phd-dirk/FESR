#include "./adler_function.hpp"

complex<double> AdlerFunction::D2(const complex<double> &s,
                                  const complex<double> mu2, const double &astau,
                                  const int &order, const int &r) {
  // double ePLT3 = 1.e2; // guesstimate
  const int i = 0, j = 1;

  complex<double> L = log(-s/mu2);
  complex<double> amu = alpha_s(sqrt(mu2), astau);
  complex<double> rmq = runMassRatio(mu2, const_.kSTau, astau);

  complex<double> m2a = (pow(const_.mq[i], 2) + pow(const_.mq[j], 2))*pow(rmq, 2);
  complex<double> m2b = r*const_.mq[i]*const_.mq[j]*pow(rmq, 2);
  complex<double> m2c = (pow(const_.mq[0], 2) + pow(const_.mq[1], 2) + pow(const_.mq[2], 2))*pow(rmq, 2);

  complex<double> sum = 0.;

  if ( order > -1 )
    sum += m2a;
  if ( order > 0 )
    sum += ((13./3.-2.*L)*m2a + 2./3.*m2b)*amu;
  if ( order > 1 )
    sum += ((4.25*pow(L, 2) - 26.*L + 23077./432. + 179./54.*const_.zeta_[3] - 520./27.*const_.zeta_[5])*m2a
            + (-17./6.*L + 769./54. - 55./27.*const_.zeta_[3] - 5./27.*const_.zeta_[5])*m2b
            + (-32./9. + 8./3.*const_.zeta_[3])*m2c
            )*pow(amu, 2);
  // if ( order > 2 )
  //   sum += ((-221./24.*pow(L, 3) + 1153./12.*pow(L, 2) + (-46253/108. - 1787./108.*zeta_[3] + 3380./27.*zeta_[5])*L
  //            + 3909929./5184. - pow(kPi, 4)/36. - 1541./648.*zeta_[3] + 26.5*pow(zeta_[3], 2) - 54265./108.*zeta_[5]
  //            + 79835./648.*zeta_[7])*m2a
  //           + (221./24.*pow(L, 2) + (-10831./108. + 715./54.*zeta_[3] + 65./54.*zeta_[5])*L + 4421./54. + ePLT3
  //              - 715./54.*zeta_[3] - 65./54.*zeta_[5])*m2b
  //           + ((208./9. - 52./3.*zeta_[3])*L - 2222./27. + 1592./27.*zeta_[3] + 4.*pow(zeta_[3], 2) - 80./27.*zeta_[5])*m2c
  //           )*pow(amu, 3);

  return 3.*sum/(4*pow(const_.kPi, 2)*s);
}

complex<double> AdlerFunction::D4(const complex<double> &s, const complex<double> &mu2,
                                  const double &astau, const double &aGGinv,
                                  const int &order, const int &r) {

  const int i = 0, j = 1;
  complex<double> L = log(-s/mu2);
  complex<double> amu = alpha_s(sqrt(mu2), astau);


  complex<double> gluonCondensate(0.0, 0.0);
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

  const double mqqa = const_.mq[i]*const_.qqinv[i] + const_.mq[j]*const_.qqinv[j];
  const double mqqb = r*(const_.mq[i]*const_.qqinv[j] + const_.mq[j]*const_.qqinv[i]);
  const double mqqs = const_.mq[0]*const_.qqinv[0] + const_.mq[1]*const_.qqinv[1] + const_.mq[2]*const_.qqinv[2];

  if (order > -1)
    quarkCondensate += 2.*mqqa;
  if (order > 0)
    quarkCondensate += (-2.*mqqa + 8./3.*mqqb + 8./27.*mqqs)*amu;
  if (order > 1)
    quarkCondensate += ((4.5*L - 131./12.)*mqqa
                        + (-6.*L + 68./3.)*mqqb
                        + (-2./3.*L - 176./243. + 8./3.*const_.zeta_[3])*mqqs
                        )*pow(amu, 2);
  if (order > 2)
    quarkCondensate += ((-81./8.*pow(L, 2) + 457./8.*L + 2.*qLT3)*mqqa
                        + (27./2.*pow(L, 2) - 338/3.*L + 8./3.*tLT3)*mqqb
                        + (1.5*pow(L, 2) + (56./27. -12.*const_.zeta_[3])*L
                           + 50407./17496. + 8./27.*pLT3 + rLT3 -20./27.*const_.zeta_[3])*mqqs
                        )*pow(amu, 3);
  quarkCondensate /= pow(s, 2);

  complex<double> m4(0.0, 0.0);
  complex<double> rmq = runMassRatio(mu2, const_.kSTau, astau);

  complex<double> m4a = (pow(const_.mq[i], 4) + pow(const_.mq[j], 4))*pow(rmq, 4);
  complex<double> m4b = r*(const_.mq[i]*pow(const_.mq[j], 3) + const_.mq[j]*pow(const_.mq[i], 3))*pow(rmq, 4);
  complex<double> m4c = (pow(const_.mq[i], 2)*pow(const_.mq[j], 2))*pow(rmq, 4);
  complex<double> m4d = (pow(const_.mq[0], 4) + pow(const_.mq[1], 4) + pow(const_.mq[2], 4))*pow(rmq, 4);


  if ( order > 1)
    m4 += (-6./7./amu + 1.5*L - 0.25)*m4a - 8./7.*m4b + 3.*m4c - m4d/14.;
  if ( order > 2)
    m4 += ((-3.*pow(L, 2) + 74./7.*L)*m4a + 32./7.*L*m4b - 12.*L*m4c + 2./7.*L*m4d)*amu;

  m4 /= pow(const_.kPi*s, 2);

  return gluonCondensate + quarkCondensate + m4;
}
double AdlerFunction::D4CInt(double s0,
                             function<complex<double>(complex<double>)> weight, const double &astau,
                             const double &aGGinv, const int &order, const int &r) {
  function<complex<double>(complex<double>)> f =
    [&](complex<double> s) -> complex<double> {
    complex<double> mu2 = s0;
    return weight(s)*D4(s0*s, mu2, astau, aGGinv, order, r);
  };

  return (3*const_.kPi*complexContourIntegral(s0, f)).real();
};
