#ifndef SRC_TYPES_H
#define SRC_TYPES_H

#include <functional>
#include <complex>

typedef std::complex<double> cmplx;
typedef std::function<cmplx(cmplx)> cmplxFunc;

#endif
