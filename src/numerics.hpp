#ifndef SRC_NUMERICS_H
#define SRC_NUMERICS_H

#include <gsl/gsl_integration.h>
#include <functional>

using std::function;

class Numerics {
 public:
  Numerics(double epsabs, double epsrel) : w_(gsl_integration_workspace_alloc(1000)), epsrel_(epsrel), epsabs_(epsabs) {}

  static double test(double x) {
    return x*x;
  }

  double integrate(function<double(double)> func) {
    double result, error;
    gsl_function F;
    F = {
      [](double d, void* vf) -> double {
        auto& f = *static_cast<std::function<double(double)>*>(vf);
        return f(d);
      },
      &func
    };
    gsl_integration_qag(&F, 0, 1, epsabs_, epsrel_, 1000, 6, w_, &result, &error);
    return result;
  }
 private:
  gsl_integration_workspace * w_;
  double epsrel_; // relative error
  double epsabs_; // absolute error
};

#endif
