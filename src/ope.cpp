#include "./ope.hpp"

typedef std::function<complex<double>(complex<double>)> cmplxFunc;

complex<double> OPE::D0(
  const complex<double> &s, const complex<double> &mu2, const double &sTau,
  const double &astau, const matrix<double> &c, const double &order
) {
  complex<double> L = log(-s/mu2);
  complex<double> amu = AlphaS::run(mu2, sTau, astau/M_PI);

  complex<double> sum(0., 0.);
  for (int n = 1; n <= order; n++) {
    for (int k = 1; k <= n; k++) {
      complex<double> powL = ( k-1 == 0 ) ? 1.0 : pow(L, k-1);
      sum += pow(amu, n)*(double)k*c(n, k)*powL;
    }
  }

  return 1/4./pow(M_PI, 2)*(c(0, 1) + sum);
}

double OPE::D0CIntFO(
  const double &s0, const Weight weight, const double &sTau,
  const double &astau, const matrix<double> &c, const double &order
) {
   cmplxFunc f =
     [&](cmplx s) -> cmplx {
       cmplx mu2(s0, 0.);
       return weight.wD(s)*D0(s0*s, mu2, sTau, astau, c, order);
     };
   return (3*M_PI*Numerics::gaussInt(f)).real();
};

double OPE::D0CIntCI(
  const double &s0, const Weight weight, const double &sTau, const double &astau,
  const matrix<double> &c, const double &order
) {
  cmplxFunc f =
    [&](complex<double> x) -> complex<double> {
    return weight.wD(x)*D0(s0*x, -x*s0, sTau, astau, c, order);
  };

  return (3*M_PI*Numerics::gaussInt(f)).real();
};

complex<double> OPE::D4(
  const complex<double> &s, const complex<double> &mu2, const double &sTau,
  const double &astau, const double &aGGinv, const int &order, const int &r,
  const std::vector<double> &mq, const Condensates &condensates
) {
  const int i = 0, j = 1;
  cmplx L = log(-s/mu2);
  cmplx amu = AlphaS::run(mu2, sTau, astau/M_PI);

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

  const double mqqa = mq[i]*condensates.qqInv_[i] + mq[j]*condensates.qqInv_[j];
  const double mqqb = r*(
    mq[i]*condensates.qqInv_[j] + mq[j]*condensates.qqInv_[i]
  );
  const double mqqs = mq[0]*condensates.qqInv_[0]
    + mq[1]*condensates.qqInv_[1]
    + mq[2]*condensates.qqInv_[2];

  if (order > -1)
    quarkCondensate += 2.*mqqa;
  if (order > 0)
    quarkCondensate += (-2.*mqqa + 8./3.*mqqb + 8./27.*mqqs)*amu;
  if (order > 1)
    quarkCondensate += (
      (4.5*L - 131./12.)*mqqa
      + (-6.*L + 68./3.)*mqqb
      + (-2./3.*L - 176./243. + 8./3.*Numerics::zeta_[3])*mqqs
    )*pow(amu, 2);
  if (order > 2)
    quarkCondensate += (
      (-81./8.*pow(L, 2) + 457./8.*L + 2.*qLT3)*mqqa
      + (27./2.*pow(L, 2) - 338/3.*L + 8./3.*tLT3)*mqqb
      + (
        1.5*pow(L, 2) + (56./27. -12.*Numerics::zeta_[3])*L
        + 50407./17496. + 8./27.*pLT3 + rLT3 -20./27.*Numerics::zeta_[3]
      )*mqqs
    )*pow(amu, 3);
  quarkCondensate /= pow(s, 2);

  // m4
  cmplx m4(0.0, 0.0);
  cmplx rmq = MQ::run(mu2, sTau, sTau, astau/M_PI);

  cmplx m4a = (pow(mq[i], 4) + pow(mq[j], 4))*pow(rmq, 4);
  cmplx m4b = r*(mq[i]*pow(mq[j], 3) + mq[j]*pow(mq[i], 3))*pow(rmq, 4);
  cmplx m4c = (pow(mq[i], 2)*pow(mq[j], 2))*pow(rmq, 4);
  cmplx m4d = (pow(mq[0], 4) + pow(mq[1], 4) + pow(mq[2], 4))*pow(rmq, 4);

  if ( order > 1)
    m4 += (-6./7./amu + 1.5*L - 0.25)*m4a - 8./7.*m4b + 3.*m4c - m4d/14.;
  if ( order > 2)
    m4 += ((-3.*pow(L, 2) + 74./7.*L)*m4a + 32./7.*L*m4b - 12.*L*m4c + 2./7.*L*m4d)*amu;

  m4 /= pow(M_PI*s, 2);

  return gluonCondensate + quarkCondensate + m4;
}

double OPE::D4CInt(
  const double &s0, const Weight &weight, const double &sTau,
  const double &astau, const double &aGGinv, const int &r,
  const std::vector<double> &mq, const Condensates &condensates
) {
  cmplxFunc fTest =
    [&](cmplx s) -> cmplx const {
    return weight.wD(s)/pow(s, 2);
  };

  double norder = 2;
  if ( abs(Numerics::gaussInt(fTest).real()) < 2.e-14) {
    norder = 3;
  }

  cmplxFunc f =
    [&](cmplx s) -> cmplx {
    cmplx mu2(s0, 0.);
    return weight.wD(s)*D4(
      s0*s, mu2, sTau, astau, aGGinv, norder, r,
      mq, condensates
    );
  };

  return (3*M_PI*Numerics::gaussInt(f)).real();
};

// rhoVA = - 10^2 C_6V/A;
complex<double> OPE::D6(
  const complex<double> &s, const double &c6
) {
  return 0.03*c6/pow(s, 3);
}
double OPE::D6CInt(
  const double &s0, const Weight &weight, const double &c6
) {
  cmplxFunc f =
    [&](cmplx s) -> cmplx {
      return weight.wD(s)*D6(s0*s, c6);
    };

  return (3*M_PI*Numerics::gaussInt(f)).real();
};

// c8VA = 10^2 C_8V/A
complex<double> OPE::D8(
  const complex<double> &s, const double &c8
) {
  return 0.04*c8/pow(s, 4);
}
double OPE::D8CInt(
  const double &s0, const Weight &weight, const double &c8
) {
  cmplxFunc f =
    [&](cmplx s) -> cmplx {
      return weight.wD(s)*D8(s0*s, c8);
    };

  return (3*M_PI*Numerics::gaussInt(f)).real();
};

// c10VA = -10^2 C_10VA
complex<double> OPE::D10(
  const complex<double> &s, const double &c10
) {
  return 0.05*c10/pow(s, 5);
}
double OPE::D10CInt(
  const double &s0, const Weight &weight, const double &c10
) {
  cmplxFunc f =
    [&](cmplx s) -> cmplx {
      return weight.wD(s)*D10(s0*s, c10);
    };

  return (3*M_PI*Numerics::gaussInt(f)).real();
};

// c12VA = 10^2 C_12VA
complex<double> OPE::D12(
  const complex<double> &s, const double &c12
) {
  return 0.06*c12/pow(s, 5);
}
double OPE::D12CInt(
  const double &s0, const Weight &weight, const double &c12
) {
  cmplxFunc f =
    [&](cmplx s) -> cmplx {
      return weight.wD(s)*D12(s0*s, c12);
    };

  return (3*M_PI*Numerics::gaussInt(f)).real();
};
