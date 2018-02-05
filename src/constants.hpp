#ifndef SRC_CONSTANTS_H
#define SRC_CONSTANTS_H

#include <vector>
#include <cmath>

using std::pow;
using std::vector;

class Constants {
 public:
  Constants() : zeta_(8) {
    zeta_[3] = 1.2020569031595942;  // coefficients.nb
    zeta_[5] = 1.036927755143369926;  // coefficients.nb
    zeta_[7] = 1.008349277381922827;
  }

  vector<double> zeta_;

  constexpr static double kBe = 17.827; // HFAG 2011
  constexpr static double kDBe = 0.04; // HFAG 2011
  constexpr static double kFPi = 92.21e-3; // PDG 2010
  constexpr static double kDFPi = 0.14e-3;
  constexpr static double kTauMass = 1.77682; // PDF 2012
  constexpr static double kPionMinusMass = 0.13957018; // M_pi^-
  constexpr static double kSEW = 1.0198; // EW radiative corr.
  constexpr static double kDSEW = 0.0006;
  constexpr static double kSTau = 3.1570893124; // kTauMass^2
  constexpr static double kAlphaTau = 0.31927; // kTauMass^2
  constexpr static double kPi = 3.141592653589793;
  constexpr static double kVud = 0.97425; // Towner, Hardy 2009
  constexpr static double kDVud = 0.00022;
  constexpr static double kDRTauVex = 0.0;
};

#endif
