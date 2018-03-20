#ifndef SRC_CONSTANTS_H
#define SRC_CONSTANTS_H

#include <vector>
#include <cmath>

using std::pow;
using std::vector;

class Constants {
 public:
  Constants() {}

  // math
  const double kPi = 3.141592653589793;
  const vector<double> zeta_ = {
    0,
    0,
    0,
    1.2020569031595942,
    0,
    1.036927755143369926,
    0,
    1.008349277381922827
  };

  // masses
  const double mumtau = 2.8e-3;
  const double mdmtau = 5.e-3;
  const double msmtau = 97e-3;
  const vector<double> mq = {mumtau, mdmtau, msmtau};

  // alpha_s
  const double kAsTauBJ = 0.3156;
  const double kATauBJ = kAsTauBJ/kPi;

  // condensates
  const double uumtau = -pow(0.272, 3);
  const double ddmtau = -pow(0.272, 3);
  const double kappa = 0.8;
  const double ssmtau = kappa*uumtau;

  const vector<double> qqinv = {
    uumtau + 3.*pow(mumtau, 3)/(7.*pow(kPi, 2))*(1./kATauBJ - 53./24.),
    ddmtau + 3.*pow(mdmtau, 3)/(7.*pow(kPi, 2))*(1./kATauBJ - 53./24.),
    ssmtau + 3.*pow(msmtau, 3)/(7.*pow(kPi, 2))*(1./kATauBJ - 53./24.)
  };

  // alpha
  const double kAlphaTau = 0.31927; // kTauMass^2

  // Various
  const double kVud = 0.97425; // Towner, Hardy 2009
  const double kDVud = 0.00022;
  const double kSEW = 1.0198; // EW radiative corr.
  const double kDSEW = 0.0006;

  constexpr static double kBe = 17.827; // HFAG 2011
  constexpr static double kDBe = 0.04; // HFAG 2011
  constexpr static double kFPi = 92.21e-3; // PDG 2010
  constexpr static double kDFPi = 0.14e-3;
  constexpr static double kTauMass = 1.77682; // PDF 2012
  constexpr static double kPionMinusMass = 0.13957018; // M_pi^-
  constexpr static double kSTau = 3.1570893124; // kTauMass^2
  constexpr static double kDRTauVex = 0.0;
};

#endif
