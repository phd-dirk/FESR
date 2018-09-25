#ifndef SRC_CONFIGURATION_HPP
#define SRC_CONFIGURATION_HPP

#include "./types.hpp"
#include "./weights.hpp"
#include <string>

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


class Configuration {
 public:
  Configuration(string configFilePath) {
    ifstream configFile(configFilePath);
    json jsonConfig;
    configFile >> jsonConfig;

    order = jsonConfig["parameters"]["order"];
    RVANormalization = jsonConfig["parameters"]["RVANormalization"];

    nc = jsonConfig["parameters"]["nc"];
    nf = jsonConfig["parameters"]["nf"];

    thMomContribs = {
      jsonConfig["scheme"], jsonConfig["thMomContribs"]["D0"], jsonConfig["thMomContribs"]["D4"],
      jsonConfig["thMomContribs"]["D68"], jsonConfig["thMomContribs"]["DV"], jsonConfig["thMomContribs"]["PionPole"]
    };

    astau = {
      jsonConfig["variables"]["astau"]["fixed"],
      jsonConfig["variables"]["astau"]["value"],
      jsonConfig["variables"]["astau"]["stepSize"]
    };
    aGGInv = {
      jsonConfig["variables"]["aGGInv"]["fixed"],
      jsonConfig["variables"]["aGGInv"]["value"],
      jsonConfig["variables"]["aGGInv"]["stepSize"]
    };
    rhoVpA = {
      jsonConfig["variables"]["rhoVpA"]["fixed"],
      jsonConfig["variables"]["rhoVpA"]["value"],
      jsonConfig["variables"]["rhoVpA"]["stepSize"]
    };
    c8VpA = {
      jsonConfig["variables"]["c8VpA"]["fixed"],
      jsonConfig["variables"]["c8VpA"]["value"],
      jsonConfig["variables"]["c8VpA"]["stepSize"]
    };
    deltaV = {
      jsonConfig["variables"]["deltaV"]["fixed"],
      jsonConfig["variables"]["deltaV"]["value"],
      jsonConfig["variables"]["deltaV"]["stepSize"]
    };
    gammaV = {
      jsonConfig["variables"]["gammaV"]["fixed"],
      jsonConfig["variables"]["gammaV"]["value"],
      jsonConfig["variables"]["gammaV"]["stepSize"]
    };
    alphaV = {
      jsonConfig["variables"]["alphaV"]["fixed"],
      jsonConfig["variables"]["alphaV"]["value"],
      jsonConfig["variables"]["alphaV"]["stepSize"]
    };
    betaV = {
      jsonConfig["variables"]["betaV"]["fixed"],
      jsonConfig["variables"]["betaV"]["value"],
      jsonConfig["variables"]["betaV"]["stepSize"]
    };
    deltaA = {
      jsonConfig["variables"]["deltaA"]["fixed"],
      jsonConfig["variables"]["deltaA"]["value"],
      jsonConfig["variables"]["deltaA"]["stepSize"]
    };
    gammaA = {
      jsonConfig["variables"]["gammaA"]["fixed"],
      jsonConfig["variables"]["gammaA"]["value"],
      jsonConfig["variables"]["gammaA"]["stepSize"]
    };
    alphaA = {
      jsonConfig["variables"]["alphaA"]["fixed"],
      jsonConfig["variables"]["alphaA"]["value"],
      jsonConfig["variables"]["alphaA"]["stepSize"]
    };
    betaA = {
      jsonConfig["variables"]["betaA"]["fixed"],
      jsonConfig["variables"]["betaA"]["value"],
      jsonConfig["variables"]["betaA"]["stepSize"]
    };

    mTau = jsonConfig["parameters"]["mTau"];
    sTau = pow(mTau, 2);

    be = jsonConfig["parameters"]["be"];
    dBe = jsonConfig["parameters"]["dBe"];

    // add weights & s0s
    for(auto const& input : jsonConfig["parameters"]["input"]) {
      Weight w(input["weight"].get<int>());
      vec s0s = input["s0s"].get<vec>();
      inputs.push_back({ w, s0s });
      momCount += s0s.size();
    }
    initializeBetaCoefficients();
    initializeAdlerCoefficients();
  }

  int dof() const {
    // model: [{w_1, [s1, s2, s3, ...]}, {w_2, [s1, s_2, ...]}, ...]
    // sum_i sum_j s_{ij}, where i is index of weight and j is the index of the corresponding s0
    int dof = 0;
    for(auto const &input: inputs) {
      vec s0s = input.s0s;
      dof += s0s.size();
    }

    if(!astau.isFixed)
      dof--;
    if(!aGGInv.isFixed)
      dof--;
    if(!rhoVpA.isFixed)
      dof--;
    if(!c8VpA.isFixed)
      dof--;

    return dof;
  }

  int order;
  double RVANormalization;
  std::vector<Input> inputs;
  uint momCount = 0;

  // QCD
  double nc;
  double nf;

 // OPE
  ThMomContribs thMomContribs;
  Variable astau, aGGInv, rhoVpA, c8VpA, deltaV, gammaV, alphaV, betaV, deltaA, gammaA, alphaA, betaA;

  // masses
  double mTau;
  double sTau;
  const double mumtau = 2.8e-3;
  const double mdmtau = 5.e-3;
  const double msmtau = 97e-3;
  const double kPionMinusMass = 0.13957018; // M_pi^-
  const double kFPi = 92.21e-3; // PDG 2010
  const vec mq = {mumtau, mdmtau, msmtau};
  const double kTauMass = 1.77682; // PDF 2012

  // RGE
  double beta[5];
  double c[6][6];

  // math
  const vec zeta = {
    0,
    0,
    0,
    1.2020569031595942,
    0,
    1.036927755143369926,
    0,
    1.008349277381922827
  };

  // alpha_s
  const double kAsTauBJ = 0.3156;
  const double kATauBJ = kAsTauBJ/M_PI;

  // condensates
  const double uumtau = -pow(0.272, 3);
  const double ddmtau = -pow(0.272, 3);
  const double kappa = 0.8;
  const double ssmtau = kappa*uumtau;

  const vec qqinv = {
    uumtau + 3.*pow(mumtau, 3)/(7.*pow(M_PI, 2))*(1./kATauBJ - 53./24.),
    ddmtau + 3.*pow(mdmtau, 3)/(7.*pow(M_PI, 2))*(1./kATauBJ - 53./24.),
    ssmtau + 3.*pow(msmtau, 3)/(7.*pow(M_PI, 2))*(1./kATauBJ - 53./24.)
  };

  // Excited resonance parameters
  const double kF1P = 2.2e-3, kM1P = 1.3, kG1P = 0.4;
  const double kF2P = 0.19e-3, kM2P = 1.8, kG2P = 0.21;

  // Various
  const double kVud = 0.97425; // Towner, Hardy 2009
  const double kDVud = 0.00022;
  const double kSEW = 1.0198; // EW radiative corr.
  const double kDSEW = 0.0006;
  double be; // HFAG
  double dBe; // HFAG
  const double kDFPi = 0.14e-3;
  const double kDRTauVex = 0.0;
  const double deltaEW = 0.001;


 private:
  void initializeBetaCoefficients() {
    beta[1] = 1./6.*(11.*nc - 2.*nf);
    beta[2] = 51./4. - 19./12.*nf;  // rgm06
    beta[3] = 2857./64. - 5033./576.*nf + 325./1728.*pow(nf, 2);  // rgm06
    beta[4] = 149753./768. + 891./32.*zeta[3]  // rgm06
      -(1078361./20736. + 1627./864.*zeta[3])*nf
      + (50065./20736. + 809./1296.*zeta[3])*pow(nf, 2)
      + 1093./93312.*pow(nf, 3);
  }

  void initializeAdlerCoefficients() {
    c[0][0] = -5./3.; c[0][1] = 1;  // rgm06
    c[1][1] = 1.; c[1][2] = 0.;  //  rgm06
    c[2][1] = 365./24. - 11.*zeta[3] - (11./12. - 2./3.*zeta[3])*nf;
    c[2][2] = -beta[1]*c[1][1]/4.; c[2][3] = 0.;  // rgm06
    c[3][1] = 87029./288. - 1103./4.*zeta[3] + 275./6.*zeta[5]
        +(-7847./216. + 262./9.*zeta[3] - 25./9.*zeta[5])*nf
        + (151./162.-19./27.*zeta[3])*pow(nf, 2);  // rgm06
    c[3][2] = -1./4.*(beta[2]*c[1][1]+2*beta[1]*c[2][1]);
    c[3][3] = pow(beta[1], 2)/12.*c[1][1]; c[3][4] = 0.;  // rgm06
    c[4][1] = 78631453./20736. - 1704247./432.*zeta[3]
        + 4185./8.*pow(zeta[3], 2) + 34165./96.*zeta[5]
        - 1995./16.*zeta[7];  // Diogo PHD
    c[4][2] = -1./4.*(beta[3]*c[1][1]+2*beta[2]*c[2][1]+3*beta[1]*c[3][1]);
    c[4][3] = beta[1]/24.*(5.*beta[2]*c[1][1]+6*beta[1]*c[2][1]);  // rgm-6
    c[4][4] = -pow(beta[1], 3)/32.*c[1][1]; c[4][5] = 0.;  // rgm06
    c[5][1] = 283.;
    c[5][2] = 1./4.*(-beta[4]*c[1][1] - 2.*beta[3]*c[2][1]-3.
                      *beta[2]*c[3][1]-4.*beta[1]*c[4][1]);
    c[5][3] = 1./24.*(12.*c[3][1]*pow(beta[1], 2)+6.*beta[1]*beta[3]*c[1][1]
                       +14.*beta[2]*beta[1]*c[2][1]+3.*pow(beta[2], 2)*c[1][1]);
    c[5][4] = 1./96.*(-12*pow(beta[1], 3)*c[2][1]
                       -13.*beta[2]*pow(beta[1], 2)*c[1][1]);
    c[5][5] = 1./80.*pow(beta[1], 4)*c[1][1];beta[1] = 11./2. - 1./3.*nf;  // rgm06
  }
};

#endif
