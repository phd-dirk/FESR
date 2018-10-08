#include "./experimentalMoments.hpp"


ExperimentalMoments::ExperimentalMoments(const string &filename, const Configuration &config)
  :config_(config), data_(Data(filename, config.RVANormalization_)),
  covMat(config_.momCount_, config_.momCount_)
{
  // cache experimental moments
  initExpMoms();

  // cache covariance matrix
  initCovMat();
}

vec ExperimentalMoments::operator () () const {
  return expMoms;
}

void ExperimentalMoments::initExpMoms() {
  for(auto const &input: config_.inputs_) {
    vec s0s = input.s0s;
    Weight w = input.weight;
    for(auto const &s0: s0s) {
      expMoms.push_back(expPlusPionMom(s0, w));
    }
  }
}

double ExperimentalMoments::expMom(const double &s0, const Weight &w) const {
  double mom = 0;
  for(int j = 0; j <= closestBinToS0(s0); j++) {
    mom += config_.sTau_/s0/config_.be_*data_.sfm2s[j]
      *wRatio(s0, w, data_.sbins[j], data_.dsbins[j]);
  }
  return mom;
}

int ExperimentalMoments::closestBinToS0(const double &s0) const {
  int i = data_.dsbins.size()-1;
  if (s0 > data_.sbins.back())
    return i;

  // -1.e-6: if exactly between two bins choose smaller one as closest
  // e.g. s0 = 2.1; sbins[71] = 2.05; sbins[72] = 2.15 => closest bin: 71
  while(abs(s0-data_.sbins[i]-1.e-6) > data_.dsbins[i]/2.)
    i--;

  return i;
}

double ExperimentalMoments::wRatio(
  const double &s0, const Weight &w, const double &sbin,
  const double &dsbin
) const {
  double binRight = sbin+dsbin/2.;
  double binLeft = sbin-dsbin/2.;
  return (
    s0/config_.sTau_*(
      (w.wD(binLeft/s0) - w.wD(binRight/s0))/
      (w.wTau(binLeft/config_.sTau_) - w.wTau(binRight/config_.sTau_))
    )
  ).real();
}


double ExperimentalMoments::pionPoleMoment(
  const double &s0, const Weight &w
) const
{
  double axialMoment = 0;
  double pseudoMoment = 0;
  axialMoment += kPiFac()/s0*w.wR(pow(config_.kPionMinusMass, 2)/s0).real();
  pseudoMoment += axialMoment*(-2.*pow(config_.kPionMinusMass, 2)/(config_.sTau_ + 2.*pow(config_.kPionMinusMass, 2)));
  return axialMoment + pseudoMoment;
}
