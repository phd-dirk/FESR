#include "./experimentalMoments.hpp"
#include "./numerics.hpp"


ExpMoms::ExpMoms(const string &filename, const Configuration &config)
  :data_(Data(filename, config.RVANormalization_))
{
  inputs_ = config.inputs_;
  sTau_ = config.sTau_;
  be_ = config.be_;
  dBe_ = config.dBe_;
  vud_ = config.vud_;
  dVud_ = config.dVud_;
  SEW_ = config.SEW_;
  dSEW_ = config.dSEW_;
  fPi_ = config.fPi_;
  dFPi_ = config.dFPi_;
  mPiM_ = config.mPiM_;
  momCount_ = config.momCount_;

  // cache experimental moments
  initExpMoms();

  // cache covariance matrix
  covMat_ = mat(config.momCount_, config.momCount_);
  invCovMat_ = mat(momCount_, momCount_);
  initCovMat();
  initInvCovMat(covMat_);
};

ExpMoms::ExpMoms(
  const string &filename,
  const std::vector<Input> &inputs,
  const double &sTau,
  const double &be,
  const double &dBe,
  const double &vud,
  const double &dVud,
  const double &SEW,
  const double &dSEW,
  const double &fPi,
  const double &dFPi,
  const double &mPiM,
  const double &RVANormalization
) : data_(Data(filename, RVANormalization)) {
  inputs_ = inputs;
  sTau_ = sTau;
  be_ = be;
  dBe_ = dBe;
  vud_ = vud;
  dVud_ = dVud;
  SEW_ = SEW;
  dSEW_ = dSEW;
  fPi_ = fPi;
  dFPi_ = dFPi;
  mPiM_ = mPiM;

  momCount_ = 0;
  for(const auto &input : inputs) {
    momCount_ += input.s0s.size();
  }

  // cache experimental moments
  initExpMoms();

  // cache covariance matrix
  covMat_ = mat(momCount_, momCount_);
  invCovMat_ = mat(momCount_, momCount_);
  initCovMat();
  initInvCovMat(covMat_);
}

vec ExpMoms::operator () () const {
  return expMoms;
}

void ExpMoms::initExpMoms() {
  for(auto const &input: inputs_) {
    vec s0s = input.s0s;
    Weight w = input.weight;
    for(auto const &s0: s0s) {
      expMoms.push_back(expPlusPionMom(s0, w));
    }
  }
}

double ExpMoms::expMom(const double &s0, const Weight &w) const {
  double mom = 0;
  for(int j = 0; j <= closestBinToS0(s0); j++) {
    mom += sTau_/s0/be_*data_.sfm2s[j]
      *wRatio(s0, w, j);
  }
  return mom;
}

int ExpMoms::closestBinToS0(const double &s0) const {
  int i = data_.dsbins.size()-1;
  if (s0 > data_.sbins.back())
    return i;

  // -1.e-6: if exactly between two bins choose smaller one as closest
  // e.g. s0 = 2.1; sbins[71] = 2.05; sbins[72] = 2.15 => closest bin: 71
  while(abs(s0-data_.sbins[i]-1.e-6) > data_.dsbins[i]/2.)
    i--;

  return i;
}

double ExpMoms::wRatio(
  const double &s0, const Weight &w, const int &bin
) const {
  double binRight = data_.sbins[bin]+data_.dsbins[bin]/2.;
  double binLeft = data_.sbins[bin]-data_.dsbins[bin]/2.;
  return (
    s0/sTau_*(
      (w.wD(binLeft/s0) - w.wD(binRight/s0))/
      (w.wTau(binLeft/sTau_) - w.wTau(binRight/sTau_))
    )
  ).real();
}


double ExpMoms::pionPoleMoment(
  const double &s0, const Weight &w
) const
{
  double axialMoment = 0;
  double pseudoMoment = 0;
  axialMoment += piFac()/s0*w.wR(pow(mPiM_, 2)/s0).real();
  pseudoMoment += axialMoment*(-2.*pow(mPiM_, 2)/(sTau_ + 2.*pow(mPiM_, 2)));
  return axialMoment + pseudoMoment;
}

double ExpMoms::expPlusPionMom(const double &s0, const Weight &w) const {
  return expMom(s0, w) + pionPoleMoment(s0, w);
}


mat ExpMoms::errMat() const {
  mat errMat(data_.binCount+2, data_.binCount+2);
  for (int i = 0; i < data_.binCount+2; i++) {
    for (int j = 0; j < data_.binCount+2; j++) {
      if (i < data_.binCount && j < data_.binCount) {
        errMat(i, j) = data_.corerrs(i, j)*data_.derrs[i]*data_.derrs[j]/1.e2;
      } else {
        errMat(i, j) = 0.;
      }
    }
  }
  errMat(data_.binCount, data_.binCount) = pow(dBe_, 2);
  errMat(data_.binCount+1, data_.binCount+1) = pow(dPiFac(), 2);
  return errMat;
}

mat ExpMoms::jacMat() const {
  mat jac(data_.binCount+2, momCount_);
  int xMom = 0;
  for(auto const &input: inputs_) {
    vec s0s = input.s0s;
    Weight w = input.weight;

    for(auto const &s0: s0s) {
      for (int j = 0; j < data_.binCount; j++) {
        if (j <= closestBinToS0(s0)) {
          jac(j, xMom) = sTau_/s0/be_*wRatio(s0, w, j);
        } else {
          jac(j, xMom) = 0.;
        }
      }
      jac(data_.binCount, xMom) = (
        pionPoleMoment(s0, w) - expPlusPionMom(s0, w)
      )/be_;
      jac(data_.binCount+1, xMom) = pionPoleMoment(s0, w)/piFac();
      xMom++;
    }

  }
  return jac;
}

void ExpMoms::initCovMat() {
  mat jac = jacMat();
  mat err = errMat();
  for(int i = 0; i < momCount_; i++) {
    for(int j = 0; j < momCount_; j++) {
      covMat_(i, j) = 0.0;
      for (int k = 0; k < data_.binCount+2; k++) {
        for (int l = 0; l < data_.binCount+2; l++) {
          covMat_(i, j) += jac(k, i)*err(k, l)*jac(l, j);
        }
      }
    }
  }
}

void ExpMoms::initInvCovMat(mat covMat) {
  // Remove correlations with R_tau,V+A in Aleph fit
  for (int i = 1; i < momCount_; i++) {
    covMat(0, i) = 0.;
    covMat(i, 0) = 0.;
  }

  // employ uncertainity of R_VA = 3.4718(72) (HFLAV 2017)
  covMat(0, 0) = pow(0.0072, 2);

  Numerics::invertMatrix(covMat, invCovMat_);
  // Numerics::invMat(covMat, invCovMat_);

  // invCovMat_(0, 0) = 19290.123456790123;
  // invCovMat_(0, 1) = 0.0;
  // invCovMat_(0, 2) = 0.0;
  // invCovMat_(0, 3) = 0.0;
  // invCovMat_(0, 4) = 0.0;
  // invCovMat_(0, 5) = 0.0;
  // invCovMat_(0, 6) = 0.0;
  // invCovMat_(0, 7) = 0.0;
  // invCovMat_(0, 8) = 0.0;

  // invCovMat_(1, 0) = 0.0;
  // invCovMat_(1, 1) = 4380565.4450900145;
  // invCovMat_(1, 2) = -19527009.072698355;
  // invCovMat_(1, 3) = 40984663.028208137;
  // invCovMat_(1, 4) = -84638945.086208180;
  // invCovMat_(1, 5) = 97316900.741565526;
  // invCovMat_(1, 6) = -49343953.720028751;
  // invCovMat_(1, 7) = 11674972.251288334;
  // invCovMat_(1, 8) = -847642.72515251662;

  // invCovMat_(2, 0) = 0.0;
  // invCovMat_(2, 1) = -19527009.072729781;
  // invCovMat_(2, 2) = 98150487.413654625;
  // invCovMat_(2, 3) = -231227676.50983095;
  // invCovMat_(2, 4) = 523320149.90033841;
  // invCovMat_(2, 5) = -619639670.88996792;
  // invCovMat_(2, 6) = 321400104.48728836;
  // invCovMat_(2, 7) = -79056190.238436818;
  // invCovMat_(2, 8) = 6598930.0692821946;

  // invCovMat_(3, 0) = 0.0;
  // invCovMat_(3, 1) = 40984663.024091452;
  // invCovMat_(3, 2) = -231227676.48426008;
  // invCovMat_(3, 3) = 616643833.28414011;
  // invCovMat_(3, 4) = -1583045922.3074424;
  // invCovMat_(3, 5) = 2001978657.2773747;
  // invCovMat_(3, 6) = -1125775522.6668854;
  // invCovMat_(3, 7) = 317010377.39433080;
  // invCovMat_(3, 8) = -36621042.135975793;

  // invCovMat_(4, 0) = 0.0;
  // invCovMat_(4, 1) = -84638945.013663828;
  // invCovMat_(4, 2) = 523320149.45882702;
  // invCovMat_(4, 3) = -1583045921.2283897;
  // invCovMat_(4, 4) = 4937255245.5470905;
  // invCovMat_(4, 5) = -7261416322.2394047;
  // invCovMat_(4, 6) = 5049585965.2160759;
  // invCovMat_(4, 7) = -1902780247.4704247;
  // invCovMat_(4, 8) = 322015652.31520826;

  // invCovMat_(5, 0) = 0.0;
  // invCovMat_(5, 1) = 97316900.522174478;
  // invCovMat_(5, 2) = -619639669.56130409;
  // invCovMat_(5, 3) = 2001978653.6812286;
  // invCovMat_(5, 4) = -7261416316.0957413;
  // invCovMat_(5, 5) = 12200521218.568886;
  // invCovMat_(5, 6) = -10108354297.680387;
  // invCovMat_(5, 7) = 4627749765.6307192;
  // invCovMat_(5, 8) = -938765886.65570676;

  // invCovMat_(6, 0) = 0.0;
  // invCovMat_(6, 1) = -49343953.430073619;
  // invCovMat_(6, 2) = 321400102.73960257;
  // invCovMat_(6, 3) = -1125775517.7485409;
  // invCovMat_(6, 4) = 5049585953.6595955;
  // invCovMat_(6, 5) = -10108354287.130436;
  // invCovMat_(6, 6) = 10220252556.225853;
  // invCovMat_(6, 7) = -5633885753.7398357;
  // invCovMat_(6, 8) = 1326788375.8412738;

  // invCovMat_(7, 0) = 0.0;
  // invCovMat_(7, 1) = 11674972.064385563;
  // invCovMat_(7, 2) = -79056189.117130280;
  // invCovMat_(7, 3) = 317010374.19440556;
  // invCovMat_(7, 4) = -1902780239.1231403;
  // invCovMat_(7, 5) = 4627749755.9647923;
  // invCovMat_(7, 6) = -5633885749.8849869;
  // invCovMat_(7, 7) = 3610024411.9811840;
  // invCovMat_(7, 8) = -951159600.84332108;

  // invCovMat_(8, 0) = 0.0;
  // invCovMat_(8, 1) = -847642.67718991451;
  // invCovMat_(8, 2) = 6598929.7827849686;
  // invCovMat_(8, 3) = -36621041.315511465;
  // invCovMat_(8, 4) = 322015650.09003186;
  // invCovMat_(8, 5) = -938765883.88217449;
  // invCovMat_(8, 6) = 1326788374.4404285;
  // invCovMat_(8, 7) = -951159600.57284117;
  // invCovMat_(8, 8) = 272107474.58081555;
}


double ExpMoms::piFac() const {
  return 24.*pow(M_PI*vud_*fPi_, 2)*SEW_;
}
double ExpMoms::dPiFac() const {
  return piFac()*sqrt(
    4.*pow(dVud_/vud_, 2)
    + pow(dSEW_/SEW_, 2)
    + 4.*pow(dFPi_/fPi_, 2)
  );
}
