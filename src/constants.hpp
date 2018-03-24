#ifndef SRC_CONSTANTS_H
#define SRC_CONSTANTS_H

#include <vector>
#include <cmath>

using std::pow;
using std::vector;

class Constants {
 public:
  Constants(const int &nc, const int &nf) : nc_(nc), nf_(nf) {
    beta_[1] = 1./6.*(11.*nc_ - 2.*nf_);
    beta_[2] = 51./4. - 19./12.*nf_;  // rgm06
    beta_[3] = 2857./64. - 5033./576.*nf_ + 325./1728.*pow(nf_, 2);  // rgm06
    beta_[4] = 149753./768. + 891./32.*zeta_[3]  // rgm06
        -(1078361./20736. + 1627./864.*zeta_[3])*nf_
        + (50065./20736. + 809./1296.*zeta_[3])*pow(nf_, 2)
        + 1093./93312.*pow(nf_, 3);

    c_[0][0] = -5./3.; c_[0][1] = 1;  // rgm06
    c_[1][1] = 1.; c_[1][2] = 0.;  //  rgm06
    c_[2][1] = 365./24. - 11.*zeta_[3] - (11./12. - 2./3.*zeta_[3])*nf_;
    c_[2][2] = -beta_[1]*c_[1][1]/4.; c_[2][3] = 0.;  // rgm06
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

  double nc_;
  double nf_;

  // RGE
  vector<double> beta_;
  vector<vector<double>> c_;

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

  // s
  const double kSTau = 3.1570893124; // kTauMass^2

  // masses
  const double mumtau = 2.8e-3;
  const double mdmtau = 5.e-3;
  const double msmtau = 97e-3;
  const double kPionMinusMass = 0.13957018; // M_pi^-
  const double kFPi = 92.21e-3; // PDG 2010
  const vector<double> mq = {mumtau, mdmtau, msmtau};
  const double kTauMass = 1.77682; // PDF 2012

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

  // Excited resonance parameters
  const double kF1P = 2.2e-3, kM1P = 1.3, kG1P = 0.4;
  const double kF2P = 0.19e-3, kM2P = 1.8, kG2P = 0.21;

  // Various
  const double kVud = 0.97425; // Towner, Hardy 2009
  const double kDVud = 0.00022;
  const double kSEW = 1.0198; // EW radiative corr.
  const double kDSEW = 0.0006;
  const double kBe = 17.827; // HFAG 2011
  const double kDBe = 0.04; // HFAG 2011
  const double kDFPi = 0.14e-3;
  const double kDRTauVex = 0.0;
};

#endif
