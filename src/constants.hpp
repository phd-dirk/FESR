#ifndef SRC_CONSTANTS_H
#define SRC_CONSTANTS_H

#include <cmath>

using std::pow;

class Constants {
 public:
  constexpr static double kBe = 17.827; // HFAG 2011
  constexpr static double kDBe = 0.04; // HFAG 2011
  constexpr static double kFPi = 92.21e-3; // PDG 2010
  constexpr static double kTauMass = 1.77682; // PDF 2012
  constexpr static double kPionMinusMass = 0.13957018; // M_pi^-
  constexpr static double kSEW = 1.0198; // EW radiative corr.
  constexpr static double kSTauMass = 3.1570893124; // kTauMass^2
  constexpr static double kPi = 3.141592653589793;
  constexpr static double kVud = 0.97425; // Towner, Hardy 2009
  constexpr static double kDRTauVex = 0.0;
};

#endif
