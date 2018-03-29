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
#include <cmath>

namespace ublas = boost::numeric::ublas;
using std::cout;
using std::endl;
using std::function;
using std::complex;
using std::cos;
using std::abs;

class Numerics {
 public:
  Numerics(Constants constants)
    : const_(constants), w_(gsl_integration_workspace_alloc(1200)),
      gaulegX(1201), gaulegW(1201) {
    // init Gaussian quadratures from Numerical recepies
    gauleg(-const_.kPi, const_.kPi, gaulegX, gaulegW, 1201);
  }

  // find root of function
  complex<double> newtonRaphson(
                                function<complex<double>(complex<double>)> f,
                                function<complex<double>(complex<double>)> df,
                                complex<double> x0,
                                double acc ) {
    complex<double> xNext = x0;
    while( abs(f(xNext)) > acc ) {
      xNext = xNext - f(xNext)/df(xNext);
    }
    return xNext;
  }


  void gauleg(const double &x1, const double &x2, vector<double> &x,
                vector<double> &w, const int &n) {
    double z1, z, pp, p3, p2, p1;
    int m = (n + 1)/2;
    double xm = 0.5*(x2+x1);
    double xl = 0.5*(x2-x1);

    for(int i = 0; i < m; i++) {
      z = cos(const_.kPi*(i + 0.75)/(n + 0.5));

      do {
        p1 = 1.;
        p2 = 0.;
        for(int j = 0; j < n; j++){
          p3 = p2;
          p2 = p1;
          p1 = ((2.*j + 1.)*z*p2 - j*p3)/(j+1);
        }
        pp = n*(z*p1 - p2)/(z*z - 1.);
        z1 = z;
        z = z1 - p1/pp;
      } while(abs(z - z1) > epsabs_);
      x[i] = xm - xl*z;
      x[n-1-i] = xm+xl*z;
      w[i] = 2.*xl/((1. - z*z)*pp*pp);
      w[n-1-i] = w[i];
    }
  }

  complex<double> gaussIntegration(function<complex<double>(complex<double>)> func) {
    complex<double> I(0., 1.);
    complex<double> sum(0.0, 0.0);
    for( int i = 0; i < 1201; i++) {
      complex<double> x = -exp(I*gaulegX[i]);
      sum += func(x)*gaulegW[i];
    }
    return sum;
  }

  double gslFixedPointLegendre(function<double(double)> func, double from, double to) {
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
      complex<double> I(0., 1.);
      return exp(I*t);
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
  const double epsrel_ = 0.; // relative error
  const double epsabs_ = 1e-13; // absolute error
  vector<double> gaulegX;
  vector<double> gaulegW;
};

#endif
