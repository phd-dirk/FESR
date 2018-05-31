#ifndef SRC_TYPES_H
#define SRC_TYPES_H

#include <functional>
#include <complex>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>

typedef std::vector<double> vec;
typedef boost::numeric::ublas::matrix<double> mat;
typedef std::complex<double> cmplx;
typedef std::function<cmplx(cmplx)> cmplxFunc;

#endif
