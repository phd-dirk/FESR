#include "./ope.hpp"

typedef std::function<complex<double>(complex<double>)> cmplxFunc;

OPE::OPE(const Configuration &config)
  :Numerics(), condensates_(config.condensates_)
{
  beta_ = Configuration::betaCoefficients(config.nc_, config.nf_);
  c_ = Configuration::adlerCoefficients(config.nf_, beta_);
  mq_ = config.mq_;
  sTau_ = config.sTau_;
  pionMinusMass_ = config.pionMinusMass_;
  fPi_ = config.fPi_;
  f1P_ = config.f1P_; m1P_ = config.m1P_; g1P_=config.g1P_;
  f2P_ = config.f2P_; m2P_ = config.m2P_; g2P_=config.g2P_;
}
OPE::OPE(
  const int &nc,
  const int &nf,
  const std::vector<double> &mq,
  const Condensates &condensates,
  const double &sTau,
  const double &pionMinusMass,
  const double &fPi,
  const double &f1P,
  const double &m1P,
  const double &g1P,
  const double &f2P,
  const double &m2P,
  const double &g2P
) : condensates_(condensates) {
  beta_ = Configuration::betaCoefficients(nc, nf);
  c_ = Configuration::adlerCoefficients(nf, beta_);
  sTau_ = sTau;
  pionMinusMass_ = pionMinusMass;
  fPi_ = fPi;
  f1P_ = f1P; m1P_ = m1P; g1P_ = g1P;
  f2P_ = f2P; m2P_ = m2P; g2P_ = g2P;
}

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
  return (3*M_PI*complexContourIntegral(f)).real();
};

double OPE::D0CIntCI(
  const double &s0, const Weight weight, const double &sTau, const double &astau,
  const double &order
) const {
  cmplxFunc f =
    [&](complex<double> x) -> complex<double> {
    return weight.wD(x)*D0(s0*x, -x*s0, sTau, astau, c_, order);
  };

  return (3*M_PI*gaussIntegration(f)).real();
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

  std::cout << std::endl;
  std::cout << mqqa << "\t" << mqqb << "\t" << mqqs << std::endl;

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
  std::cout << std::endl << "quarkCondensate \t" << quarkCondensate << std::endl;

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
    return weight.wD(s)*D4(
      s0*s, mu2, sTau, astau, aGGinv, norder, r,
      mq_, condensates_
    );
  };

  return (3*M_PI*complexContourIntegral(f)).real();
};

complex<double> OPE::D68(
  const complex<double> &s, const double &rho, const double &c8
) const {
  return 0.03*rho/pow(s, 3) + 0.04*c8/pow(s, 4);
}

double OPE::D68CInt(
  const double &s0, const Weight &weight, const double &rho, const double &c8
) const {
  cmplxFunc f =
    [&](cmplx s) -> cmplx {
      return weight.wD(s)*D68(s0*s, rho, c8);
    };

  return (3*M_PI*complexContourIntegral(f)).real();
};

double OPE::deltaP(const double &s0, const Weight &weight) const {
  double spi = pow(pionMinusMass_, 2);
  double pionPole = -4.*pow(fPi_, 2)/s0*spi/(sTau_ + 2.*spi)
    *weight.wR(spi/s0).real();
  double xth = 9.*spi/s0;

  function<double(double)> f =
    [&](const double &s) -> double {
      double x = s0*s;
      double rhores =  2./pow(x, 2)
        *(
          pow(f1P_, 2)*pow(m1P_, 4)*breitwigner(x, m1P_, g1P_)
          + pow(f2P_, 2)*pow(m2P_, 4)*breitwigner(x, m2P_, g2P_)
        );
      return weight.wR(s).real()*2.*x/(sTau_ + 2.*x)*rhores;
    };

  return 4.*pow(M_PI, 2)*( pionPole - adaptiveIntegrate(f, xth, 1.));
}

double OPE::breitwigner(
  const double &s, const double &mbw, const double &gbw
) const {
  return mbw*gbw/M_PI/(pow(s - pow(mbw, 2), 2)+ pow(mbw, 2)*pow(gbw, 2));
}
