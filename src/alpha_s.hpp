#ifndef SRC_ALPHAS_H
#define SRC_ALPHAS_H

#include "./types.hpp"
#include "./configuration.hpp"
#include "numerics.hpp"
#include <iostream>
#include <cmath>

// Mutlidimensional Root-finding
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>

// test multidimensional root-finding

class AlphaS {
 public:
  static cmplx runAlpha(const cmplx &mup, const cmplx &muq, const cmplx &aq);
  // Computes ap(mup) from integrating the \beta function in the complex plane in
  // Mathematica (runAlpha.nb). The initial condition a(muq)=aq is used.
  static cmplx run(
    const cmplx &q2, const cmplx &p2, const cmplx &ap
  );

 private:
  static int alpha_f(const gsl_vector *x, void *params, gsl_vector *f);
  static int alpha_df(const gsl_vector *x, void *params, gsl_matrix *J);
  static int alpha_fdf(
    const gsl_vector *x,
    void *params,
    gsl_vector *f,
    gsl_matrix *J
  );

  struct rparams {
    cmplx mup;
    cmplx muq;
    cmplx aq;
  };
};

#endif
