#include "./experimentalMoments.hpp"


ExpMoms::ExpMoms(const string &filename, const Configuration &config)
  :config_(config), data_(Data(filename, config.RVANormalization_)),
  covMat(config_.momCount_, config_.momCount_)
{
  inputs_ = config.inputs_;
  sTau_ = config.sTau_;
  be_ = config.be_;
  pionMinusMass_ = config.kPionMinusMass;

  momCount_ = config.momCount_;

  // cache experimental moments
  initExpMoms();

  // cache covariance matrix
  initCovMat();
}
// ExpMoms::ExpMoms(
//   const string &filename,
//   const std::vector<Input> &inputs,
//   const double &sTau,
//   const double &be,
//   const double &RVANormalization
// ) : data_(Data(filename, RVANormalization)) {
//   inputs_ = inputs;
//   sTau_ = sTau;
//   be_ = be;

//   // cache experimental moments
//   initExpMoms();

//   // cache covariance matrix
//   initCovMat();
// }

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
      *wRatio(s0, w, data_.sbins[j], data_.dsbins[j]);
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
  const double &s0, const Weight &w, const double &sbin,
  const double &dsbin
) const {
  double binRight = sbin+dsbin/2.;
  double binLeft = sbin-dsbin/2.;
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
  axialMoment += kPiFac()/s0*w.wR(pow(pionMinusMass_, 2)/s0).real();
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
  errMat(data_.binCount, data_.binCount) = pow(config_.dBe_, 2);
  errMat(data_.binCount+1, data_.binCount+1) = pow(kDPiFac(), 2);
  return errMat;
}

mat ExpMoms::jacMat() const {
  mat jac(data_.binCount+2, momCount_);
  int xMom = 0;
  for(auto const &input: inputs_) {
    vec s0s = input.s0s;
    Weight w = input.weight;
    for(auto const &s0: s0s) {
      for (int j = 0; j < data_.binCount+2; j++) {
        if (j <= closestBinToS0(s0)) {
          jac(j, xMom) = sTau_/s0/be_*wRatio(s0, w, data_.sbins[j], data_.dsbins[j]);
        } else {
          jac(j, xMom) = 0.;
        }
        jac(data_.binCount, xMom) = (
          pionPoleMoment(s0, w) - expPlusPionMom(s0, w)
        )/be_;
        jac(data_.binCount+1, xMom) = pionPoleMoment(s0, w)/kPiFac();
      }
      xMom++;
    }
  }
  return jac;
}

void ExpMoms::initCovMat() {
  mat jac = jacMat();
  mat err = errMat();
  for(uint i = 0; i < momCount_; i++) {
    for(uint j = 0; j < momCount_; j++) {
      covMat(i, j) = 0.0;
      for (int k = 0; k < data_.binCount+2; k++) {
        for (int l = 0; l < data_.binCount+2; l++) {
          covMat(i, j) += jac(k, i)*err(k, l)*jac(l, j);
        }
      }
    }
  }
}

double ExpMoms::kPiFac() const {
  return 24.*pow(M_PI*config_.vud_*config_.kFPi, 2)*config_.kSEW;
}
double ExpMoms::kDPiFac() const {
  return kPiFac()*sqrt(
    4.*pow(config_.dVud_/config_.vud_, 2)
    + pow(config_.kDSEW/config_.kSEW, 2)
    + 4.*pow(config_.kDFPi/config_.kFPi, 2)
  );
}
