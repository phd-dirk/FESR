#ifndef SRC_VA_DELTAS_H
#define SRC_VA_DELTAS_H
const int &nc, const int &nf, const int &order,
  const vector<double> &s0s, function<complex<double>(complex<double>)> weight) :

class Deltas: public TheoreticalMoments {
 public:
  Deltas(const int &nc, const int &nf, const int &order,
         const vector<double> &s0s,
         function<complex<double>(complex<double>)> weight)
    : TheoreticalMoments(nc, nf, order, s0s, weight) {
  }

  double VpAD0() {
    return  
  }
};

#endif
