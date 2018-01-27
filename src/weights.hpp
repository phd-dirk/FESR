#include <cmath>
#include <complex>

using std::pow;
using std::complex;

complex<double> wR00(complex<double> x) {
  return pow((1. - x), 2)*(1. + 2.*x);
}

complex<double> wD00(complex<double> x) {
  return pow((1. - x), 3)*(1. + x);
}
