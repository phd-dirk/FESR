#include "./adler_function.hpp"

complex<double> AdlerFunction::D4(const complex<double> &s, const complex<double> &mu2,
                                  const double &astau, const double &aGGinv, const int &r) {

  const int i = 0, j = 1;
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
