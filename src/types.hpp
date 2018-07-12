#ifndef SRC_TYPES_H
#define SRC_TYPES_H

#include <math.h> // defines M_PI
#include <functional>
#include <complex>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include "json.hpp"
#include <fstream>
#include <string>


// vectors
typedef std::vector<double> vec;
typedef std::vector<int> intVec;

typedef boost::numeric::ublas::matrix<double> mat;
typedef std::complex<double> cmplx;
typedef std::function<cmplx(cmplx)> cmplxFunc;
typedef nlohmann::json json;
typedef std::string string;
typedef std::ifstream ifstream;


#endif
