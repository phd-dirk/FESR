#ifndef SRC_NUMERICS_H
#define SRC_NUMERICS_H

#include "constants.hpp"
#include <gsl/gsl_integration.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <functional>
#include <complex>
#include <iostream>

namespace ublas = boost::numeric::ublas;
using namespace std::complex_literals;
using std::cout;
using std::endl;
using std::function;
using std::complex;

class Numerics {
 public:
  Numerics(const double &epsabs, const double &epsrel, Constants constants)
    : const_(constants), w_(gsl_integration_workspace_alloc(1200)), epsrel_(epsrel), epsabs_(epsabs) {}

  double gaussLegendre(function<double(double)> func, double from, double to) {
    double result;
    gsl_integration_fixed_workspace * q;
    const gsl_integration_fixed_type * T = gsl_integration_fixed_legendre;
    q = gsl_integration_fixed_alloc(T, 1200, from, to, 0.0, 0.0);

    gsl_function F;
    F = {
      [](double d, void* vf) -> double {
        auto& f = *static_cast<std::function<double(double)>*>(vf);
        return f(d);
      },
      &func
    };

    gsl_integration_fixed(&F, &result, q);

    return result;
  }

  double adaptiveIntegrate(function<double(double)> func, double from, double to) {
    double result, error;
    gsl_function F;
    F = {
      [](double d, void* vf) -> double {
        auto& f = *static_cast<std::function<double(double)>*>(vf);
        return f(d);
      },
      &func
    };
    gsl_integration_qag(&F, from, to, epsabs_, epsrel_, 1200, 6, w_, &result, &error);
    // cout << "error \t" << error << endl;
    return result;
  }

  complex<double> integrateComplex(function<complex<double>(double)> func, double from, double to) {
    auto funcReal = [func](double t) {
      return func(t).real();
    };
    auto funcImag = [func](double t) {
      return func(t).imag();
    };

    double cintReal = adaptiveIntegrate(funcReal, from, to);
    double cintImag = adaptiveIntegrate(funcImag, from, to);

    return complex<double>(cintReal, cintImag);
  }

  complex<double> complexContourIntegral(function<complex<double>(complex<double>)> f) {
    auto gamma = [](double t) {
      return exp(1i*t);
    };

    auto func = [&](double t) {
      return f(gamma(t));
    };

    return integrateComplex(func, -const_.kPi, const_.kPi);
  }

  // from https://gist.github.com/lilac/2464434
  template<class T>
  bool invertMatrix (const ublas::matrix<T>& input, ublas::matrix<T>& inverse) {
    using namespace boost::numeric::ublas;
    typedef permutation_matrix<std::size_t> pmatrix;
    // create a working copy of the input
    matrix<T> A(input);
    // create a permutation matrix for the LU-factorization
    pmatrix pm(A.size1());

    // perform LU-factorization
    int res = lu_factorize(A,pm);
    if( res != 0 ) return false;

    // create identity matrix of "inverse"
    inverse.assign(ublas::identity_matrix<T>(A.size1()));

    // backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);

    return true;
  }



 private:
  Constants const_;
  gsl_integration_workspace * w_;
  double epsrel_; // relative error
  double epsabs_; // absolute error
};

#endif
