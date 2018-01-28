#ifndef SRC_PI_POLE_H
#define SRC_PI_POLE_H

#include "./constants.hpp"
#include "./data.hpp"
#include <cmath>
#include <functional>
#include <complex>
#include <vector>

using std::function;
using std::complex;
using std::vector;

const double kPiFac = 24.*pow(kPi*kVud*kFPi, 2)*kSEW;

inline double pionPoleSpectralMoment(double s0, double pionMinusMass, double sTauMass,
                              function<complex<double>(complex<double>)> weight,
                              vector<double> sbins, vector<double> dsbins
                              ) {
  double axialPionPionPoleSpectralMoment = 0;
  double pseudoscalarPionPionPoleSpectralMoment = 0;
  axialPionPionPoleSpectralMoment += kPiFac/s0*weight(pow(kPionMinusMass, 2)/s0).real();
  pseudoscalarPionPionPoleSpectralMoment += axialPionPionPoleSpectralMoment*
    (-2.*pow(pionMinusMass, 2)/(sTauMass + 2.*pow(pionMinusMass, 2)));
  return axialPionPionPoleSpectralMoment + pseudoscalarPionPionPoleSpectralMoment;
}

#endif
