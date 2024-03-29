#ifndef SRC_EXPERIMENTALMOMENTS_HPP
#define SRC_EXPERIMENTALMOMENTS_HPP

#include "./types.hpp"
#include "./configuration.hpp"
#include "./data.hpp"

class ExpMoms {
 public:
  // Init. data, vector of all wRatios, exp. Spectral Moments, error Matrix
  // and Covariance matrix
  //
  // The Spectral Moments and Covariance matrix can then bet exported via the
  // public getter functions
  ExpMoms(const string &filename, const Configuration &config);
  ExpMoms(
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
  );

  vec operator ()() const;

  void initExpMoms();

  // return the experimental spectral moment
  double expMom(const double &s0, const Weight &w) const;
  double pionPoleMoment(const double &s0, const Weight &w) const;
  double expPlusPionMom(const double &s0, const Weight &w) const;

  // Selects the closest bin number from s0
  // If s0 is exactly between two bins we select the smaller one
  int closestBinToS0(const double &s0) const;

  // returns weight ratio
  double wRatio(
    const double &s0, const Weight &w, const int &bin
  ) const;

  // returns the error matrix
  mat errMat() const;

  // returns the Jacobian Matrix
  mat jacMat() const;

  // returns the covariance matrix
  void initCovMat();
  void initInvCovMat(mat covMat);

  double piFac() const;
  double dPiFac() const;

  std::vector<Input> inputs_;
  double sTau_;
  double be_;
  double dBe_;
  double vud_;
  double dVud_;
  double SEW_;
  double dSEW_;
  double fPi_;
  double dFPi_;
  double mPiM_;
  int momCount_;
  const Data data_;
  vec expMoms;
  mat covMat_;
  mat invCovMat_;
}; // END ExpMoms

#endif
