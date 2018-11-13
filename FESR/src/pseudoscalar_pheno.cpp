#include "./pseudoscalar_pheno.hpp"

double PSPheno::deltaP(
  const double &s0, const Weight &weight, const double &sTau, const double &mPiM,
  const double &fPi, const double &f1P, const double &m1P, const double &g1P,
  const double &f2P, const double &m2P, const double &g2P
) {
  double spi = pow(mPiM, 2);
  double pionPole = -4.*pow(fPi, 2)/s0*spi/(sTau + 2.*spi)
    *weight.wR(spi/s0).real();
  double xth = 9.*spi/s0;

  std::function<double(double)> f =
    [&](const double &s) -> double {
      double x = s0*s;
      double rhores =  2./pow(x, 2)
        *(
          pow(f1P, 2)*pow(m1P, 4)*breitwigner(x, m1P, g1P)
          + pow(f2P, 2)*pow(m2P, 4)*breitwigner(x, m2P, g2P)
        );
      return weight.wR(s).real()*2.*x/(sTau + 2.*x)*rhores;
    };

  return 4.*pow(M_PI, 2)*( pionPole - Numerics::adaptiveIntegrate(f, xth, 1.));
}

double PSPheno::breitwigner(
  const double &s, const double &mbw, const double &gbw
) {
  return mbw*gbw/M_PI/(pow(s - pow(mbw, 2), 2)+ pow(mbw, 2)*pow(gbw, 2));
}
