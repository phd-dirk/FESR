#ifndef SRC_NUMERICS_H
#define SRC_NUMERICS_H

#include "./types.hpp"
#include <gsl/gsl_integration.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/spirit/include/karma.hpp>
#include <functional>
#include <complex>
#include <iostream>
#include <cmath>

namespace ublas = boost::numeric::ublas;
namespace karma = boost::spirit::karma;
using std::cout;
using std::endl;
using std::function;
using std::complex;
using std::cos;
using std::abs;


class Numerics {
 public:
  // find root of function
  complex<double> newtonRaphson(
    function<complex<double>(complex<double>)> f,
    function<complex<double>(complex<double>)> df,
    complex<double> x0,
    double acc
  ) const {
    complex<double> xNext = x0;
    while( abs(f(xNext)) > acc ) {
      xNext = xNext - f(xNext)/df(xNext);
      cout << "N \t" << xNext << endl;
    }
    return xNext;
  }

  static void gauleg(const double &x1, const double &x2, vec &x, vec &w, const int &n);
  static complex<double> gaussInt(
    function<complex<double>(complex<double>)> func
  );

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

  static double nonAdaptiveIntegrate(func f, double from, double to);
  static double adaptiveIntegrate(func f, double from, double to);
  static double gaussLegendre(func f, double from, double to);

  double semiInfInt(function<double(double)> func, double from) const
  {
    double result, error;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(5000);
    gsl_function F;
    F = {
      [](double d, void* vf) -> double {
        auto& f = *static_cast<std::function<double(double)>*>(vf);
        return f(d);
      },
      &func
    };
    gsl_integration_qagiu(&F, from, epsabs_, epsrel_, 5000, w, &result, &error);
    return result;
  }

  double qawf(func func, cDbl &omega, cDbl &from) const {
    double result, error;
    gsl_integration_workspace * wSpace = gsl_integration_workspace_alloc(1200);
    gsl_integration_workspace * cycleSpace = gsl_integration_workspace_alloc(1200);

    gsl_integration_qawo_table *qawoTable = gsl_integration_qawo_table_alloc(omega, 1.0, GSL_INTEG_SINE, 1200);
    gsl_function F;
    F = {
      [](double d, void* vf) -> double {
        auto& f = *static_cast<std::function<double(double)>*>(vf);
        return f(d);
      },
      &func
    };

    gsl_integration_qawf(&F, from, epsabs_, 1200, wSpace, cycleSpace, qawoTable, &result, &error);

    return result;
  }

  static std::complex<double> integrateComplex(
    function<complex<double>(double)> func, double from, double to
  );

  static complex<double> complexContourIntegral(function<complex<double>(complex<double>)> f) {
    auto gamma = [](double t) {
      complex<double> I(0., 1.);
      return exp(I*t);
    };

    auto func = [&](double t) {
      return f(gamma(t));
    };

    return integrateComplex(func, -M_PI, M_PI);
  }

  // from https://gist.github.com/lilac/2464434
  static bool invertMatrix (
    const ublas::matrix<double>& input,
    ublas::matrix<double>& inverse
  );
  static bool invMat(
    const ublas::matrix<double> &mat,
    ublas::matrix<double> &invMat
  );
  static void testInvMat(
    const ublas::matrix<double> &mat,
    const ublas::matrix<double> &invMat
  );

  static const vec zeta_;
 private:
  const double epsrel_ = 0.; // relative error
  const double epsabs_ = 1e-10; // absolute error
};

#endif
