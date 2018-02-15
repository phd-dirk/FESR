#ifndef SRC_ALPHAS_H
#define SRC_ALPHAS_H

// #include <complex>
// #include <iostream>

// using std::complex;
// using std::ostream;
// typedef complex< double > state_type;

// struct stuart_landau {
//   double m_eta;
//   double m_alpha;

//   stuart_landau( double eta = 1.0, double alpha = 1.0 )
//       : m_eta ( eta ) , m_alpha( alpha ) {}

//   void operator()( const state_type &x, state_type &dxdt, double t) const {
//     const complex< double > I( 0.0, 1.0 );
//     dxdt = ( 1.0 + m_eta * I ) * x - ( 1.0 + m_alpha * I ) * norm( x ) * x;
//   }
// };

// struct streaming_observer {
//   ostream &m_out;

//   streaming_observer( ostream &out ) : m_out( out ) { }

//   template< class State >
//   void operator()( const State &x, double t ) const {
//     m_out << t;
//     m_out << "\t" << x.real() << "\t" << x.imag();
//     m_out << "\n";
//   }
// };

#endif
