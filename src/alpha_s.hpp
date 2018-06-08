#ifndef SRC_ALPHAS_H
#define SRC_ALPHAS_H

#include "./constants.hpp"
#include "./types.hpp"
#include "numerics.hpp"
#include <iostream>
#include <cmath>

// Mutlidimensional Root-finding
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>

// test multidimensional root-finding

class AlphaS: public Numerics {
 public:
  AlphaS(const Constants &constants) : Numerics(constants), const_(constants) {}

  complex<double> operator ()(const cmplx &mup, const cmplx &muq, const cmplx &aq) const {
    const gsl_multiroot_fdfsolver_type *T;
    gsl_multiroot_fdfsolver *s;

    int status;
    size_t i, iter = 0;

    const size_t n = 2;
    struct rparams p = {mup, muq, aq};
    gsl_multiroot_function_fdf f = {AlphaS::alpha_f, AlphaS::alpha_df, AlphaS::alpha_fdf, n, &p};

    double x_init[2] = { 0.07, -0.02 };
    gsl_vector *x = gsl_vector_alloc(n);
    gsl_vector_set(x, 0, x_init[0]);
    gsl_vector_set(x, 1, x_init[1]);

    T = gsl_multiroot_fdfsolver_hybridsj;
    s = gsl_multiroot_fdfsolver_alloc(T, n);
    gsl_multiroot_fdfsolver_set(s, &f, x);

    // print_state(iter, s);

    do {
      iter++;
      status = gsl_multiroot_fdfsolver_iterate(s);
      // print_state(iter, s);

      if(status)
        break;

      status = gsl_multiroot_test_residual (s->f, 1e-15);
    } while( status == GSL_CONTINUE && iter < 1000);

    // printf ("status = %s\n", gsl_strerror(status));

    const cmplx ap(gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1));
    gsl_multiroot_fdfsolver_free(s);
    gsl_vector_free(x);

    return ap;
  }

  struct rparams {
    cmplx mup;
    cmplx muq;
    cmplx aq;
  };
  static int alpha_f(const gsl_vector *x, void *params, gsl_vector *f) {
    const cmplx mup = ((rparams *) params)->mup;
    const cmplx muq = ((rparams *) params)->muq;
    const cmplx aq = ((rparams *) params)->aq;

    const double mup1 = mup.real();
    const double mup2 = mup.imag();
    const double muq1 = muq.real();
    const double muq2 = muq.imag();
    const double aq1 = aq.real();
    const double aq2 = aq.imag();

    const double ap1 = gsl_vector_get (x, 0);
    const double ap2 = gsl_vector_get (x, 1);

    using std::arg;

    const double y0 = (-0.22222222222222224*ap1)/(pow(ap1,2) + pow(ap2,2)) +
      (0.22222222222222224*aq1)/(pow(aq1,2) + pow(aq2,2)) -
      0.16227475612833714*arg(cmplx(-0.14395661700847226,-0.30687590223497346) + ap1 + cmplx(0,1)*ap2) + 
      0.16227475612833714*arg(cmplx(-0.14395661700847226,0.30687590223497346) + ap1 + cmplx(0,1)*ap2) - 
      0.03722833922971071*arg(cmplx(0.329423287332271,-0.21280498219449714) + ap1 + cmplx(0,1)*ap2) + 
      0.03722833922971071*arg(cmplx(0.329423287332271,0.21280498219449714) + ap1 + cmplx(0,1)*ap2) + 
      0.16227475612833714*arg(cmplx(-0.14395661700847226,-0.30687590223497346) + aq1 + cmplx(0,1)*aq2) - 
      0.16227475612833714*arg(cmplx(-0.14395661700847226,0.30687590223497346) + aq1 + cmplx(0,1)*aq2) + 
      0.03722833922971071*arg(cmplx(0.329423287332271,-0.21280498219449714) + aq1 + cmplx(0,1)*aq2) - 
      0.03722833922971071*arg(cmplx(0.329423287332271,0.21280498219449714) + aq1 + cmplx(0,1)*aq2) + 
      0.012337599968934259*log(pow(-0.14395661700847226 + ap1,2) + pow(-0.30687590223497346 + ap2,2)) + 
      0.08642783212983117*log(pow(0.329423287332271 + ap1,2) + pow(-0.21280498219449714 + ap2,2)) - 
      0.19753086419753088*log(pow(ap1,2) + pow(ap2,2)) + 
      0.08642783212983117*log(pow(0.329423287332271 + ap1,2) + pow(0.21280498219449714 + ap2,2)) + 
      0.012337599968934259*log(pow(-0.14395661700847226 + ap1,2) + pow(0.30687590223497346 + ap2,2)) - 
      0.012337599968934259*log(pow(-0.14395661700847226 + aq1,2) + pow(-0.30687590223497346 + aq2,2)) - 
      0.08642783212983117*log(pow(0.329423287332271 + aq1,2) + pow(-0.21280498219449714 + aq2,2)) + 
      0.19753086419753088*log(pow(aq1,2) + pow(aq2,2)) - 
      0.08642783212983117*log(pow(0.329423287332271 + aq1,2) + pow(0.21280498219449714 + aq2,2)) - 
      0.012337599968934259*log(pow(-0.14395661700847226 + aq1,2) + pow(0.30687590223497346 + aq2,2)) - 
      log(sqrt(pow(muq1,2) + pow(muq2,2))/sqrt(pow(mup1,2) + pow(mup2,2)));
    const double y1 = (0.22222222222222224*ap2)/(pow(ap1,2) + pow(ap2,2)) - (0.22222222222222224*aq2)/(pow(aq1,2) + pow(aq2,2)) + 
      0.024675199937868517*arg(cmplx(-0.14395661700847226,-0.30687590223497346) + ap1 + cmplx(0,1)*ap2) + 
      0.024675199937868517*arg(cmplx(-0.14395661700847226,0.30687590223497346) + ap1 + cmplx(0,1)*ap2) - 
      0.39506172839506176*arg(ap1 + cmplx(0,1)*ap2) + 
      0.17285566425966234*arg(cmplx(0.329423287332271,-0.21280498219449714) + ap1 + cmplx(0,1)*ap2) + 
      0.17285566425966234*arg(cmplx(0.329423287332271,0.21280498219449714) + ap1 + cmplx(0,1)*ap2) - 
      0.024675199937868517*arg(cmplx(-0.14395661700847226,-0.30687590223497346) + aq1 + cmplx(0,1)*aq2) - 
      0.024675199937868517*arg(cmplx(-0.14395661700847226,0.30687590223497346) + aq1 + cmplx(0,1)*aq2) + 
      0.39506172839506176*arg(aq1 + cmplx(0,1)*aq2) - 
      0.17285566425966234*arg(cmplx(0.329423287332271,-0.21280498219449714) + aq1 + cmplx(0,1)*aq2) - 
      0.17285566425966234*arg(cmplx(0.329423287332271,0.21280498219449714) + aq1 + cmplx(0,1)*aq2) - 
      arg((muq1 + cmplx(0,1)*muq2)/(mup1 + cmplx(0,1)*mup2)) + 
      0.08113737806416857*log(pow(-0.14395661700847226 + ap1,2) + pow(-0.30687590223497346 + ap2,2)) + 
      0.018614169614855354*log(pow(0.329423287332271 + ap1,2) + pow(-0.21280498219449714 + ap2,2)) - 
      0.018614169614855354*log(pow(0.329423287332271 + ap1,2) + pow(0.21280498219449714 + ap2,2)) - 
      0.08113737806416857*log(pow(-0.14395661700847226 + ap1,2) + pow(0.30687590223497346 + ap2,2)) - 
      0.08113737806416857*log(pow(-0.14395661700847226 + aq1,2) + pow(-0.30687590223497346 + aq2,2)) - 
      0.018614169614855354*log(pow(0.329423287332271 + aq1,2) + pow(-0.21280498219449714 + aq2,2)) + 
      0.018614169614855354*log(pow(0.329423287332271 + aq1,2) + pow(0.21280498219449714 + aq2,2)) + 
      0.08113737806416857*log(pow(-0.14395661700847226 + aq1,2) + pow(0.30687590223497346 + aq2,2));

    gsl_vector_set (f, 0, y0);
    gsl_vector_set (f, 1, y1);

    return GSL_SUCCESS;
  }
  static int alpha_df(const gsl_vector *x, void *params, gsl_matrix *J) {
    const double ap1 = gsl_vector_get (x, 0);
    const double ap2 = gsl_vector_get (x, 1);

    const double df00 = -0.053350370503906966/(pow(-0.14395661700847226 + ap1,2) + pow(-0.30687590223497346 + ap2,2)) + 
      (0.024675199937868517*ap1)/(pow(-0.14395661700847226 + ap1,2) + pow(-0.30687590223497346 + ap2,2)) + 
      0.049020305087512026/(pow(0.329423287332271 + ap1,2) + pow(-0.21280498219449714 + ap2,2)) + 
      (0.17285566425966234*ap1)/(pow(0.329423287332271 + ap1,2) + pow(-0.21280498219449714 + ap2,2)) + 
      (0.16227475612833714*ap2)/(pow(-0.14395661700847226 + ap1,2) + pow(-0.30687590223497346 + ap2,2)) + 
      (0.03722833922971071*ap2)/(pow(0.329423287332271 + ap1,2) + pow(-0.21280498219449714 + ap2,2)) + 
      (0.22222222222222224*pow(ap1,2))/pow(pow(ap1,2) + pow(ap2,2),2) - 
      (0.22222222222222224*pow(ap2,2))/pow(pow(ap1,2) + pow(ap2,2),2) - 
      (0.39506172839506176*ap1)/(pow(ap1,2) + pow(ap2,2)) + 
      0.049020305087512026/(pow(0.329423287332271 + ap1,2) + pow(0.21280498219449714 + ap2,2)) + 
      (0.17285566425966234*ap1)/(pow(0.329423287332271 + ap1,2) + pow(0.21280498219449714 + ap2,2)) - 
      (0.03722833922971071*ap2)/(pow(0.329423287332271 + ap1,2) + pow(0.21280498219449714 + ap2,2)) - 
      0.053350370503906966/(pow(-0.14395661700847226 + ap1,2) + pow(0.30687590223497346 + ap2,2)) + 
      (0.024675199937868517*ap1)/(pow(-0.14395661700847226 + ap1,2) + pow(0.30687590223497346 + ap2,2)) - 
      (0.16227475612833714*ap2)/(pow(-0.14395661700847226 + ap1,2) + pow(0.30687590223497346 + ap2,2));
    const double df01 = -0.015788300674348502/(pow(-0.14395661700847226 + ap1,2) + pow(-0.30687590223497346 + ap2,2)) + 
      (0.16227475612833714*ap1)/(pow(-0.14395661700847226 + ap1,2) + pow(-0.30687590223497346 + ap2,2)) + 
      0.049048428445967664/(pow(0.329423287332271 + ap1,2) + pow(-0.21280498219449714 + ap2,2)) + 
      (0.03722833922971071*ap1)/(pow(0.329423287332271 + ap1,2) + pow(-0.21280498219449714 + ap2,2)) - 
      (0.024675199937868517*ap2)/(pow(-0.14395661700847226 + ap1,2) + pow(-0.30687590223497346 + ap2,2)) - 
      (0.17285566425966234*ap2)/(pow(0.329423287332271 + ap1,2) + pow(-0.21280498219449714 + ap2,2)) - 
      (0.4444444444444445*ap1*ap2)/pow(pow(ap1,2) + pow(ap2,2),2) +
      (0.39506172839506176*ap2)/(pow(ap1,2) + pow(ap2,2)) -
      0.049048428445967664/(pow(0.329423287332271 + ap1,2) + pow(0.21280498219449714 + ap2,2)) - 
      (0.03722833922971071*ap1)/(pow(0.329423287332271 + ap1,2) + pow(0.21280498219449714 + ap2,2)) - 
      (0.17285566425966234*ap2)/(pow(0.329423287332271 + ap1,2) + pow(0.21280498219449714 + ap2,2)) + 
      0.015788300674348502/(pow(-0.14395661700847226 + ap1,2) + pow(0.30687590223497346 + ap2,2)) - 
      (0.16227475612833714*ap1)/(pow(-0.14395661700847226 + ap1,2) + pow(0.30687590223497346 + ap2,2)) - 
      (0.024675199937868517*ap2)/(pow(-0.14395661700847226 + ap1,2) + pow(0.30687590223497346 + ap2,2));
    const double df10 = -df01;
    const double df11 = df00;

    gsl_matrix_set (J, 0, 0, df00);
    gsl_matrix_set (J, 0, 1, df01);
    gsl_matrix_set (J, 1, 0, df10);
    gsl_matrix_set (J, 1, 1, df11);

    return GSL_SUCCESS;
  }
  static int alpha_fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J) {
    alpha_f(x, params, f);
    alpha_df(x, params, J);

    return GSL_SUCCESS;
  }
  void print_state(unsigned iter, gsl_multiroot_fdfsolver *s) const {
    printf ("iter = %3u x = % .10f % .10f "
            "f(x) = % .3e % .3e\n",
            iter,
            gsl_vector_get (s->x, 0),
            gsl_vector_get (s->x, 1),
            gsl_vector_get (s->f, 0),
            gsl_vector_get (s->f, 1));
  }

  // complex<double> operator ()(const complex<double> &q2,const complex<double> &p2, const complex<double> &ap) const {
  //   return 0.);
  // }

private:
  Constants const_;
};

#endif
