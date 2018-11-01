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


namespace Numerics {
  // find root of function
  extern complex<double> newtonRaphson(
    function<complex<double>(complex<double>)> f,
    function<complex<double>(complex<double>)> df,
    complex<double> x0,
    double acc
  );
  extern void gauleg(const double &x1, const double &x2, vec &x, vec &w, const int &n);
  extern std::vector<double> gaulegX;
  extern std::vector<double> gaulegW;
  extern complex<double> gaussInt(
    function<complex<double>(complex<double>)> func
  );

  extern double gslFixedPointLegendre(function<double(double)> func, double from, double to);

  extern double nonAdaptiveIntegrate(func f, double from, double to);
  extern double adaptiveIntegrate(func f, double from, double to);
  extern double gaussLegendre(func f, double from, double to);


  extern std::complex<double> integrateComplex(
    function<complex<double>(double)> func, double from, double to
  );

  extern complex<double> complexContourIntegral(function<complex<double>(complex<double>)> f);

  // from https://gist.github.com/lilac/2464434
  extern bool invertMatrix (
    const ublas::matrix<double>& input,
    ublas::matrix<double>& inverse
  );
  extern bool invMat(
    const ublas::matrix<double> &mat,
    ublas::matrix<double> &invMat
  );
  extern void testInvMat(
    const ublas::matrix<double> &mat,
    const ublas::matrix<double> &invMat
  );

  extern const std::vector<double> zeta_;
}

#endif
