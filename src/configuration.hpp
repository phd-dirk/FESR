#ifndef SRC_CONFIGURATION_HPP
#define SRC_CONFIGURATION_HPP

#include "./types.hpp"
#include "./weights.hpp"
#include "./numerics.hpp"
#include "./condensates.hpp"
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
    const double &asTau,
    const double &be,
    const double &dBe,
    const double &vud,
    const double &dVud,
    const double &SEW,
    const double &pionMinusMass_,
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
  double pionMinusMass_; // M_pi^-
  std::vector<double> mq_;
  const double kTauMass = 1.77682; // PDF 2012

  // alpha_s(mTau)
  double asTau_;

  // condensates
  Condensates condensates_;

  // RGE
  std::vector<double> beta_;

  // Excited resonance parameters
  double f1P_, m1P_, g1P_;
  double f2P_, m2P_, g2P_;

  // Various
  double vud_, dVud_;
  double SEW_, dSEW_ = 0.0006; // EW radiative corr.
  double be_, dBe_; // HFAG
  double fPi_;
  const double dFPi_ = 0.14e-3;
  const double kDRTauVex = 0.0;
  const double deltaEW = 0.001;

  // minuit
  double tolerance;

  static std::vector<double> betaCoefficients(const int &nc, const int &nf);
  static matrix<double> adlerCoefficients(const int &nf, const std::vector<double> &beta);

  void log();
};

#endif
