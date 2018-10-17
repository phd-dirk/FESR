#ifndef SRC_CONFIGURATION_HPP
#define SRC_CONFIGURATION_HPP

#include "./types.hpp"
#include "./weights.hpp"
#include "./numerics.hpp"
#include <string>
#include <boost/numeric/ublas/matrix.hpp>

struct Variable {
  bool isFixed;
  double value;
  double stepSize;
};

struct ThMomContribs {
  string scheme;
  bool D0;
  bool D4;
  bool D68;
  bool DV;
  bool PionPole;
};

struct Input {
  Weight weight;
  std::vector<double> s0s;
};

using boost::numeric::ublas::matrix;

class Configuration {
 public:
  Configuration(string configFilePath);
  Configuration(
    const double &be,
    const double &dBe,
    const double &vud,
    const double &dVud,
    const std::vector<Input> &inputs,
    const double &mTau,
    const int &nc,
    const int &nf,
    const int &order,
    const double &RVANormalization,
    const ThMomContribs &thMomContribs
  );

  int dof() const;

  int order_;
  double RVANormalization_;
  std::vector<Input> inputs_;
  uint momCount_ = 0;

  // QCD
  int nf_, nc_;

 // OPE
  ThMomContribs thMomContribs_;
  Variable astau, aGGInv, rhoVpA, c8VpA, deltaV, gammaV, alphaV, betaV, deltaA, gammaA, alphaA, betaA;

  // masses
  double mTau_, sTau_;
  const double mumtau = 2.8e-3, mdmtau = 5.e-3, msmtau = 97e-3;
  const double pionMinusMass_ = 0.13957018; // M_pi^-
  const std::vector<double> mq_ = {mumtau, mdmtau, msmtau};
  const double kTauMass = 1.77682; // PDF 2012

  // RGE
  std::vector<double> beta_;

  // alpha_s
  const double kAsTauBJ = 0.3156;
  const double kATauBJ = kAsTauBJ/M_PI;

  // condensates
  const double uumtau = -pow(0.272, 3);
  const double ddmtau = -pow(0.272, 3);
  const double kappa = 0.8;
  const double ssmtau = kappa*uumtau;

  const std::vector<double> qqInv_ = {
    uumtau + 3.*pow(mumtau, 3)/(7.*pow(M_PI, 2))*(1./kATauBJ - 53./24.),
    ddmtau + 3.*pow(mdmtau, 3)/(7.*pow(M_PI, 2))*(1./kATauBJ - 53./24.),
    ssmtau + 3.*pow(msmtau, 3)/(7.*pow(M_PI, 2))*(1./kATauBJ - 53./24.)
  };

  // Excited resonance parameters
  const double f1P_ = 2.2e-3, m1P_ = 1.3, g1P_ = 0.4;
  const double f2P_ = 0.19e-3, m2P_ = 1.8, g2P_ = 0.21;

  // Various
  double vud_, dVud_;
  const double SEW_ = 1.0198, dSEW_ = 0.0006; // EW radiative corr.
  double be_, dBe_; // HFAG
  const double fPi_ = 92.21e-3; // PDG 2010
  const double dFPi_ = 0.14e-3;
  const double kDRTauVex = 0.0;
  const double deltaEW = 0.001;

  // minuit
  double tolerance;

  static std::vector<double> betaCoefficients(const int &nc, const int &nf);
  static matrix<double> adlerCoefficients(const int &nf, const std::vector<double> &beta);
};

#endif
