#import "./data.hpp"
#import "./s0_sets.hpp"
#import <vector>
#import <functional>
#import <complex>

using std::vector;
using std::string;
using std::function;
using std::complex;

struct State {
public:
  vector<double> s0s;
  Data data;
  function<complex<double>(complex<double>)> weight;
  State(vector<double> s0s, int dataSize, string dataFile,
        function<complex<double>(complex<double>)> weight) :
    s0s(s0s), data(dataSize, dataFile), weight(weight) {}
};

void renormalizeState(State &state, const double &factor) {
  state.data.sfm2s = renormalize(factor, state.data.sfm2s);
  state.data.derrs = renormalize(factor, state.data.derrs);
}
