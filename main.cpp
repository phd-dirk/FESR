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

#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/integrate/integrate_const.hpp>
using boost::numeric::odeint::runge_kutta4;
using boost::numeric::odeint::integrate_const;


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
  double m_eta;
  double m_alpha;

  stuart_landau( double eta = 1.0, double alpha = 1.0 )
      : m_eta ( eta ) , m_alpha( alpha ) {}

  void operator()( const state_type &x, state_type &dxdt, double t) const {
    const complex< double > I( 0.0, 1.0 );
    dxdt = ( 1.0 + m_eta * I ) * x - ( 1.0 + m_alpha * I ) * norm( x ) * x;
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

  state_type x = complex< double > ( 1.0, 0.0 );

  const double dt = 0.1;

  typedef runge_kutta4< state_type > stepper_type;

  integrate_const( stepper_type(), stuart_landau( 2.0, 1.0 ), x, 0.0, 10.0, dt, streaming_observer( cout ) );

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
