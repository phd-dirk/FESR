#include "./src/data.hpp"
#include "./src/s0_sets.hpp"
#include "./src/state.hpp"
#include "./src/weights.hpp"
#include "./src/theoretical_moments.hpp"
#include "./src/numerics.hpp"
#include "./src/alphas.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <gsl/gsl_integration.h>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <string>
#include <complex>

#include <boost/numeric/odeint/stepper/generation/make_controlled.hpp>
#include <boost/numeric/odeint/stepper/generation/generation_runge_kutta_cash_karp54.hpp> // neccessary for make_controlled( stepper )
#include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
using boost::numeric::odeint::make_controlled;
using boost::numeric::odeint::runge_kutta_cash_karp54;
using boost::numeric::odeint::integrate_adaptive;


#include "./lib/cRunDec/CRunDec.h"
using std::cout;
using std::endl;
using std::complex;
using std::log;
using std::sqrt;
using std::printf;
using std::ostream;

std::string projectRoot = "/Users/knowledge/Developer/PhD/FESR";

inline double test(double x) {
  return x*x;
}

typedef complex< double > state_type;

struct stuart_landau {
  double m_beta;

  stuart_landau( double beta = 9./2. ): m_beta ( beta ) {}

  void operator()( const state_type &x, state_type &dxdt, double t) const {
    dxdt = -1./t * m_beta * pow(x, 2);
  }
};

struct streaming_observer {
  ostream &m_out;

  streaming_observer( ostream &out ) : m_out( out ) { }

  template< class State >
  void operator()( const State &x, double t ) const {
    m_out << t;
    m_out << "\t" << x.real() << "\t" << x.imag();
    m_out << "\n";
  }
};

int main () {
  cout.precision(17);

  state_type x = complex< double > ( 0.101627, 0.0 );

  const double dt = 0.01;

  // typedef runge_kutta_cash_karp54< state_type > stepper_type;
  auto stepper = make_controlled( 1.0e-3, 1.0e-3, runge_kutta_cash_karp54< state_type >() );

  integrate_adaptive( stepper, stuart_landau( 9./2. ), x, 1.77686, 2.0, dt, streaming_observer( cout ) );

  // ExperimentalMoments expMom(projectRoot+"/aleph.json", 0.99363, s0Set, wD00, wD00);
  // State state(s0Set, projectRoot+"/aleph.json", wR00);
  // renormalizeState(state, 0.99363);
  // cout << state.weight(complex<double>(1., 2.)) << endl;

  // AdlerFunction adler(3, 3, 4);

  // cout << adler.contourIntegral(3., wD00) << endl;
  // cout << adler.contourIntegral(5., wD00) << endl;


  // gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

  // double result, error;
  // double expected = -4.0;
  // double alpha = 1.0;

  // gsl_function F;
  // F.function = &f;
  // F.params = &alpha;

  // gsl_integration_qag(&F, 0, 1, 1e-13, 0, 1000, 6, w, &result, &error);

  // printf("result = % .18f\n", result);

  // Numerics numerics(1e-7, 1e-7);
  // cout << numerics.integrate(Numerics::test) << endl;


  return 0;
}
