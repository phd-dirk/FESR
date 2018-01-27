#import "./data.hpp"
#import "./s0_sets.hpp"
#import <vector>

using std::vector;
using std::string;

struct State {
public:
  vector<double> s0s;
  Data data;
  State(vector<double> s0s, int dataSize, string dataFile) : s0s(s0s), data(dataSize, dataFile) {}
};

void renormalizeState(State &state, const double &factor) {
  state.data.sfm2 = renormalize(factor, state.data.sfm2);
  state.data.derr = renormalize(factor, state.data.derr);
}
