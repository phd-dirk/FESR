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
  pionMinusMass_ = config.kPionMinusMass;
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
  const double &pionMinusMass,
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
  pionMinusMass_ = pionMinusMass;

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
  axialMoment += piFac()/s0*w.wR(pow(pionMinusMass_, 2)/s0).real();
  pseudoMoment += axialMoment*(-2.*pow(pionMinusMass_, 2)/(sTau_ + 2.*pow(pionMinusMass_, 2)));
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

  double sum = 0.0;
  std::cout << "mat" << std::endl << std::setprecision(15);
  for(int i = 0; i < covMat.size1(); i++) {
    for(int j = 0; j < covMat.size2(); j++) {
      sum += covMat(i, j);
    }
  }
  std::cout << "sum : " << sum << std::endl;
   // Numerics::invertMatrix(covMat, invCovMat_);
  Numerics::invMat(covMat, invCovMat_);
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
