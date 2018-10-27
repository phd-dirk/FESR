#include "./configuration.hpp"

Configuration::Configuration(string configFilePath)
{
  ifstream configFile(configFilePath);
  json jsonConfig;
  configFile >> jsonConfig;

  nf_ = jsonConfig["parameters"]["nf"];
  nc_ = jsonConfig["parameters"]["nc"];

  order_ = jsonConfig["parameters"]["order"];
  RVANormalization_ = jsonConfig["parameters"]["RVANormalization"];

  thMomContribs_ = {
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


  mTau_ = jsonConfig["parameters"]["mTau"];
  sTau_ = pow(mTau_, 2);
  mq_ = {
    jsonConfig["parameters"]["muMTau"],
    jsonConfig["parameters"]["mdMTau"],
    jsonConfig["parameters"]["msMTau"],
  };

  asTau_ = jsonConfig["parameters"]["asTau"];

  be_ = jsonConfig["parameters"]["be"];
  dBe_ = jsonConfig["parameters"]["dBe"];
  vud_ = jsonConfig["parameters"]["vud"];
  dVud_ = jsonConfig["parameters"]["dVud"];
  SEW_ = jsonConfig["parameters"]["SEW"];
  mPiM_ = jsonConfig["parameters"]["mPiM"];
  fPi_ = jsonConfig["parameters"]["fPi"];
  f1P_ = jsonConfig["parameters"]["f1P"];
  m1P_ = jsonConfig["parameters"]["m1P"];
  g1P_ = jsonConfig["parameters"]["g1P"];
  f2P_ = jsonConfig["parameters"]["f2P"];
  m2P_ = jsonConfig["parameters"]["m2P"];
  g2P_ = jsonConfig["parameters"]["g2P"];

  // add weights & s0s
  for(auto const& input : jsonConfig["parameters"]["input"]) {
    Weight w(input["weight"].get<int>());
    vec s0s = input["s0s"].get<vec>();
    inputs_.push_back({ w, s0s });
    momCount_ += s0s.size();
  }

  beta_ = betaCoefficients(jsonConfig["parameters"]["nc"], jsonConfig["parameters"]["nf"]);

  std::vector<double> qqMTau = {
    jsonConfig["condensates"]["uuMTau"],
    jsonConfig["condensates"]["ddMTau"],
    jsonConfig["condensates"]["ssMTau"]
  };
  condensates_ = Condensates(asTau_, qqMTau);

  // minuit
  tolerance = jsonConfig["tolerance"];
}

Configuration::Configuration (
  const double &asTau,
  const double &be,
  const double &dBe,
  const double &vud,
  const double &dVud,
  const double &SEW,
  const double &mPiM,
  const double &fPi,
  const double &f1P,
  const double &m1P,
  const double &g1P,
  const double &f2P,
  const double &m2P,
  const double &g2P,
  const std::vector<double> &mq,
  const std::vector<double> &qqMTau,
  const std::vector<Input> &inputs,
  const double &mTau,
  const int &nc,
  const int &nf,
  const int &order,
  const double &RVANormalization,
  const ThMomContribs &thMomContribs
) {
  inputs_ = inputs;

  nc_ = nc;
  nf_ = nf;

  order_ = order;
  RVANormalization_ = RVANormalization;

  thMomContribs_ = thMomContribs;

  mTau_ = mTau;
  sTau_ = pow(mTau, 2);
  mq_ = mq;

  asTau_ = asTau;

  be_ = be;
  dBe_ = dBe;
  vud_ = vud;
  dVud_ = dVud;
  SEW_ = SEW;
  mPiM_ = mPiM;
  fPi_ = fPi;
  f2P_ = f2P, m2P_ = m2P, g2P_ = g2P;

  beta_ = betaCoefficients(nc, nf);

  condensates_ = Condensates(asTau_, qqMTau);

  // moment count
  for(auto const &input: inputs) {
    momCount_ += input.s0s.size();
  }
}

int Configuration::dof() const {
  // model: [{w_1, [s1, s2, s3, ...]}, {w_2, [s1, s_2, ...]}, ...]
  // sum_i sum_j s_{ij}, where i is index of weight and j is the index of the corresponding s0
  int dof = 0;
  for(auto const &input: inputs_) {
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

// private
std::vector<double> Configuration::betaCoefficients(const int &nc, const int &nf) {
  std::vector<double> beta(5);
  beta[1] = 1./6.*(11.*nc - 2.*nf);
  beta[2] = 51./4. - 19./12.*nf;  // rgm06
  beta[3] = 2857./64. - 5033./576.*nf + 325./1728.*pow(nf, 2);  // rgm06
  beta[4] = 149753./768. + 891./32.*Numerics::zeta_[3]  // rgm06
    -(1078361./20736. + 1627./864.*Numerics::zeta_[3])*nf
    + (50065./20736. + 809./1296.*Numerics::zeta_[3])*pow(nf, 2)
    + 1093./93312.*pow(nf, 3);

  return beta;
}

matrix<double> Configuration::adlerCoefficients(const int &nf, const std::vector<double> &beta) {
  matrix<double> c(6, 6);

  c(0, 0) = -5./3.; c(0, 1) = 1;  // rgm06
  c(1, 1) = 1.; c(1, 2) = 0.;  //  rgm06
  c(2, 1) = 365./24. - 11.*Numerics::zeta_[3] - (11./12. - 2./3.*Numerics::zeta_[3])*nf;
  c(2, 2) = -beta[1]*c(1, 1)/4.; c(2, 3) = 0.;  // rgm06
  c(3, 1) = 87029./288. - 1103./4.*Numerics::zeta_[3] + 275./6.*Numerics::zeta_[5]
    +(-7847./216. + 262./9.*Numerics::zeta_[3] - 25./9.*Numerics::zeta_[5])*nf
    + (151./162.-19./27.*Numerics::zeta_[3])*pow(nf, 2);  // rgm06
  c(3, 2) = -1./4.*(beta[2]*c(1, 1)+2*beta[1]*c(2, 1));
  c(3, 3) = pow(beta[1], 2)/12.*c(1, 1); c(3, 4) = 0.;  // rgm06
  c(4, 1) = 78631453./20736. - 1704247./432.*Numerics::zeta_[3]
    + 4185./8.*pow(Numerics::zeta_[3], 2) + 34165./96.*Numerics::zeta_[5]
    - 1995./16.*Numerics::zeta_[7];  // Diogo PHD
  c(4, 2) = -1./4.*(beta[3]*c(1, 1)+2*beta[2]*c(2, 1)+3*beta[1]*c(3, 1));
  c(4, 3) = beta[1]/24.*(5.*beta[2]*c(1, 1)+6*beta[1]*c(2, 1));  // rgm-6
  c(4, 4) = -pow(beta[1], 3)/32.*c(1, 1); c(4, 5) = 0.;  // rgm06
  c(5, 1) = 283.;
  c(5, 2) = 1./4.*(-beta[4]*c(1, 1) - 2.*beta[3]*c(2, 1)-3.
                   *beta[2]*c(3, 1)-4.*beta[1]*c(4, 1));
  c(5, 3) = 1./24.*(12.*c(3, 1)*pow(beta[1], 2)+6.*beta[1]*beta[3]*c(1, 1)
                    +14.*beta[2]*beta[1]*c(2, 1)+3.*pow(beta[2], 2)*c(1, 1));
  c(5, 4) = 1./96.*(-12*pow(beta[1], 3)*c(2, 1)
                    -13.*beta[2]*pow(beta[1], 2)*c(1, 1));
  c(5, 5) = 1./80.*pow(beta[1], 4)*c(1, 1);  // rgm06
  return c;
}

