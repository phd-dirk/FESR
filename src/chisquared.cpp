#include "./chisquared.hpp"

// Public
Chi2::Chi2(Configuration config)
  :config_(config),
   expMom_(
     ExpMoms("/Users/knowledge/Developer/PhD/FESR/aleph.json", config)
   ) {}

double Chi2::operator() ( const double *xx) const {
  // init fit parameters
  double astau = xx[0];
  double aGGinv = xx[1];
  double rhoD6VpA = xx[2];
  double c8D8VpA = xx[3];
  double deV = xx[4];
  double gaV = xx[5];
  double alV = xx[6];
  double beV = xx[7];
  double deA = xx[8];
  double gaA = xx[9];
  double alA = xx[10];
  double beA = xx[11];

  return chi2(
    astau, aGGinv, rhoD6VpA, c8D8VpA,
    deV, gaV, alV, beV,
    deA, gaA, alA, beA,
    config_.order_, config_.sTau_, config_.mPiM_, config_.fPi_,
    config_.f1P_, config_.m1P_, config_.g1P_,
    config_.f2P_, config_.m2P_, config_.g2P_,
    config_.nc_, config_.nf_, config_.mq_, config_.vud_, config_.SEW_, config_.condensates_,
    config_.inputs_, config_.thMomContribs_, expMom_
  );
}

double Chi2::operator ()(
  const double &astau,
  const double &aGGinv,
  const double &rho,
  const double &c8,
  const double &deV,
  const double &gaV,
  const double &alV,
  const double &beV,
  const double &deA,
  const double &gaA,
  const double &alA,
  const double &beA,
  const double &order
) const {
  return chi2(
    astau, aGGinv, rho, c8,
    deV, gaV, alV, beV,
    deA, gaA, alA, alA,
    order, config_.sTau_, config_.mPiM_, config_.fPi_,
    config_.f1P_, config_.m1P_, config_.g1P_,
    config_.f2P_, config_.m2P_, config_.g2P_,
    config_.nc_, config_.nf_, config_.mq_, config_.vud_, config_.SEW_, config_.condensates_,
    config_.inputs_, config_.thMomContribs_, expMom_
  );
}

double Chi2::chi2(
  const double &astau,
  const double &aGGinv,
  const double &rho,
  const double &c8,
  const double &deV,
  const double &gaV,
  const double &alV,
  const double &beV,
  const double &deA,
  const double &gaA,
  const double &alA,
  const double &beA,
  const double &order,
  const double &sTau,
  const double &mPiM,
  const double &fPi,
  const double &f1P,
  const double &m1P,
  const double &g1P,
  const double &f2P,
  const double &m2P,
  const double &g2P,
  const int &nc,
  const int &nf,
  const std::vector<double> &mq,
  const double &vud,
  const double &SEW,
  const Condensates &condensates,
  const std::vector<Input> &inputs,
  const ThMomContribs &thMomContribs,
  const ExpMoms &expMoms
) {
  double chi = 0;

  vec thMoms = Chi2::calcThMoms(
    astau, aGGinv, rho, c8, order,
    deV, gaV, alV, beV, deA, gaA, alA, beA,
    sTau, mPiM, fPi,f1P, m1P, g1P, f2P, m2P, g2P,
    nc, nf, mq, vud, SEW, condensates, inputs, thMomContribs
  );

  double momCount = 0;
  for(auto const &input: inputs) {
    momCount += input.s0s.size();
  }

  vec momDiff(momCount);
  for(uint i = 0; i < momCount; i++) {
    momDiff[i] = expMoms()[i] - thMoms[i];
  }


  // std::vector<double> chi2n(momCount);
  double itr = 1;
  for(uint k = 0; k < momCount; k++) {
    // chi = 0.0;
    for(uint l = 0; l < momCount; l++) {
      chi += momDiff[k] * expMoms.invCovMat_(k, l) * momDiff[l];
    }
    // chi2n[k] = chi;
  }


  // chi = 0.0;
  // for(uint k = 0; k < momCount; k++) {
  //   chi += chi2n[k];
  // }

  return chi;
}

std::vector<double> Chi2::calcThMoms(
  const double &astau,
  const double &aGGinv,
  const double &rhoVpA,
  const double &c8VpA,
  const double &order,
  const double &deV,
  const double &gaV,
  const double &alV,
  const double &beV,
  const double &deA,
  const double &gaA,
  const double &alA,
  const double &beA,
  const double &sTau,
  const double &mPiM,
  const double &fPi,
  const double &f1P,
  const double &m1P,
  const double &g1P,
  const double &f2P,
  const double &m2P,
  const double &g2P,
  const int &nc,
  const int &nf,
  const std::vector<double> &mq,
  const double &vud,
  const double &SEW,
  const Condensates &condensates,
  const std::vector<Input> &inputs,
  const ThMomContribs &thMomContribs
) {
  double momCount = 0;
  for(auto const &input: inputs) {
    momCount += input.s0s.size();
  }

  std::vector<double> thMoms(momCount);
  std::vector<std::future<double>> ftrs(momCount);
  ThMoms th(nc, nf, mq, condensates, inputs, thMomContribs, vud, SEW);

  int i = 0;
  for(auto const& input: inputs) {
    vec s0s = input.s0s;
    Weight w = input.weight;
    for(auto const& s0: s0s) {
      ftrs[i] = std::async(
        &ThMoms::operator(),
        &th, s0, w, astau, aGGinv, rhoVpA, c8VpA, order, sTau,
        deV, gaV, alV, beV, deA, gaA, alA, beA,
        mPiM, fPi,f1P, m1P, g1P, f2P, m2P, g2P
      );
      i++;
    }
  }

  for(uint i=0; i<momCount; i++) {
    thMoms[i] = ftrs[i].get();
  }

  return thMoms;
}
