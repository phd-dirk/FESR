#include "./chisquared.hpp"

// Public
Chi2::Chi2(Configuration config)
  : expMom_(
     ExpMoms("/Users/knowledge/Developer/PhD/FESR/aleph.json", config)
   ) {
  order_ = config.order_;
  sTau_ = config.sTau_;
  mPiM_ = config.mPiM_;
  fPi_ = config.fPi_;
  f1P_ = config.f1P_;
  m1P_ = config.m1P_;
  g1P_ = config.g1P_;
  f2P_ = config.f2P_;
  m2P_ = config.m2P_;
  g2P_ = config.g2P_;
  c_ = Configuration::adlerCoefficients(
    config.nf_, Configuration::betaCoefficients( config.nc_, config.nf_ )
  );
  mq_ = config.mq_;
  vud_ = config.vud_;
  SEW_ = config.SEW_;
  condensates_ = config.condensates_;
  inputs_ = config.inputs_;
  thMomContribs_ = config.thMomContribs_;
}

Chi2::Chi2(
  const int &nc,
  const int &nf,
  const double &order,
  const double &sTau,
  const double &mPiM,
  const double &fPi,
  const double &dFPi,
  const double &f1P,
  const double &m1P,
  const double &g1P,
  const double &f2P,
  const double &m2P,
  const double &g2P,
  const std::vector<double> &mq,
  const double &be,
  const double &dBe,
  const double &vud,
  const double &dVud,
  const double &SEW,
  const double &dSEW,
  const double &asTauFix,
  const std::vector<double> &qqMTau,
  const std::vector<Input> &inputs,
  const ThMomContribs &thMomContribs,
  const double &RVANormalization
): expMom_(
  "/Users/knowledge/Developer/PhD/FESR/aleph.json",
  inputs,
  sTau,
  be,
  dBe,
  vud,
  dVud,
  SEW,
  dSEW,
  fPi,
  dFPi,
  mPiM,
  RVANormalization
) {
  order_ = order;
  sTau_ = sTau;
  mPiM_ = mPiM;
  fPi_ = fPi;
  f1P_ = f1P;
  m1P_ = m1P;
  g1P_ = g1P;
  f2P_ = f2P;
  m2P_ = m2P;
  g2P_ = g2P;
  c_ = Configuration::adlerCoefficients(
    nf, Configuration::betaCoefficients( nc, nf )
  );
  mq_ = mq;
  vud_ = vud;
  SEW_ = SEW;
  condensates_ = Condensates( asTauFix, qqMTau, mq );
  inputs_ = inputs;
  thMomContribs_ = thMomContribs;
}

double Chi2::operator() ( const double *xx) const {
  // init fit parameters
  double astau = xx[0];
  double aGGinv = xx[1];
  double c6 = xx[2];
  double c8 = xx[3];
  double c10 = xx[4];
  double c12 = xx[5];
  double deV = xx[6];
  double gaV = xx[7];
  double alV = xx[8];
  double beV = xx[9];
  double deA = xx[10];
  double gaA = xx[11];
  double alA = xx[12];
  double beA = xx[13];

  return chi2(
    astau, aGGinv, c6, c8, c10, c12,
    deV, gaV, alV, beV,
    deA, gaA, alA, beA,
    order_, sTau_, mPiM_, fPi_,
    f1P_, m1P_, g1P_,
    f2P_, m2P_, g2P_,
    c_, mq_, vud_, SEW_, condensates_,
    inputs_, thMomContribs_, expMom_
  );
}

double Chi2::operator ()(
  const double &astau,
  const double &aGGinv,
  const double &rho,
  const double &c8,
  const double &c10,
  const double &c12,
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
    astau, aGGinv, rho, c8, c10, c12,
    deV, gaV, alV, beV,
    deA, gaA, alA, alA,
    order, sTau_, mPiM_, fPi_,
    f1P_, m1P_, g1P_,
    f2P_, m2P_, g2P_,
    c_, mq_, vud_, SEW_, condensates_,
    inputs_, thMomContribs_, expMom_
  );
}

double Chi2::chi2(
  const double &astau,
  const double &aGGinv,
  const double &rho,
  const double &c8,
  const double &c10,
  const double &c12,
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
  const matrix<double> &c,
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
    astau, aGGinv, rho, c8, c10, c12, order,
    deV, gaV, alV, beV, deA, gaA, alA, beA,
    sTau, mPiM, fPi,f1P, m1P, g1P, f2P, m2P, g2P,
    c, mq, vud, SEW, condensates, inputs, thMomContribs
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
  const double &c8,
  const double &c10,
  const double &c12,
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
  const matrix<double> &c,
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
  ThMoms th;

  int i = 0;
  for(auto const& input: inputs) {
    vec s0s = input.s0s;
    Weight w = input.weight;
    for(auto const& s0: s0s) {
      ftrs[i] = std::async(
        &ThMoms::calc,
        s0, w, astau, aGGinv, rhoVpA, c8, c10, c12, order, sTau,
        deV, gaV, alV, beV, deA, gaA, alA, beA,
        mPiM, fPi,f1P, m1P, g1P, f2P, m2P, g2P, c, mq, condensates,
        vud, SEW, thMomContribs
      );
      i++;
    }
  }

  for(uint i=0; i<momCount; i++) {
    thMoms[i] = ftrs[i].get();
  }

  return thMoms;
}

double Chi2::chi2_spec_end(const double *xx, const int num_bins) const {
  double rho = xx[0];
  double rho_first = xx[1];
  double rho_second = xx[2];

  double stau = pow(1.77686, 2);
  double sew = 1.0198;
  double vud = 0.97420;
  double be = 17.815;

  const Data data = Data("/Users/knowledge/Developer/PhD/FESR/aleph.json", 0.99743669);
  const int dataSize = data.sbins.size();

  boost::numeric::ublas::matrix<double> corMat(num_bins, num_bins);
  for(int i=0; i<num_bins; i++) {
    for(int j=0; j<num_bins; j++) {
      corMat(i, j) = data.corerrs(i+dataSize-num_bins-1, j+dataSize-num_bins-1);
    }
  }
  boost::numeric::ublas::matrix<double> covMat(num_bins, num_bins);
  for(int i=0; i<num_bins; i++) {
    for(int j=0; j<num_bins; j++) {
      covMat(i, j) = corMat(i, j)
        *data.derrs[i+dataSize-num_bins-1]/data.dsbins[i+dataSize-num_bins-1]
        *data.derrs[j+dataSize-num_bins-1]/data.dsbins[j+dataSize-num_bins-1]/100.0;
    }
  }
  // std::cout << covMat << std::endl;

  auto wTau = [&](auto s)
  {
    return pow(1.0 - s/stau, 2)*(1 + 2*s/stau);
  };

  auto func = [&](double s, double rho, double rhoFirst, double rhoSecond)
  {
    return 2*wTau(s)/stau*3*pow(vud, 2)*sew*be*(
      rho
      + rhoFirst*(s - stau)
      + rhoSecond*1.0/2.0*pow(s - stau, 2)
    );
  };

  boost::numeric::ublas::matrix<double> invCovMat(num_bins, num_bins);
  Numerics::invertMatrix(covMat, invCovMat);

  double chi = 0.0;
  for(int i=0; i<num_bins; i++) {
    for(int j=0; j<num_bins; j++) {
      double xi = data.sbins[i+dataSize-num_bins-1];
      double yi = data.sfm2s[i+dataSize-num_bins-1]/data.dsbins[i+dataSize-num_bins-1];
      double xj = data.sbins[j+dataSize-num_bins-1];
      double yj = data.sfm2s[j+dataSize-num_bins-1]/data.dsbins[j+dataSize-num_bins-1];
      chi += (yi - func(xi, rho, rho_first, rho_second))
        *invCovMat(i, j)*(yj - func(xj, rho, rho_first, rho_second));
    }
  }

  return chi;
}
