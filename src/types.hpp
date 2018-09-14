#ifndef SRC_TYPES_H
#define SRC_TYPES_H

#include <math.h> // defines M_PI
#include <functional>
#include <complex>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <nlohmann/json.hpp>
#include <fstream>
#include <string>


// vectors
typedef std::vector<double> vec;
typedef std::vector<int> intVec;

typedef const double cDbl;
typedef boost::numeric::ublas::matrix<double> mat;
typedef std::complex<double> cmplx;
typedef std::function<cmplx(cmplx)> cmplxFunc;
typedef std::function<double(double)> func;
typedef nlohmann::json json;
typedef std::string string;
typedef std::ifstream ifstream;


#endif
