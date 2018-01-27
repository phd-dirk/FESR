#import "data.hpp"
#import "s0_sets.hpp"
#import <vector>

using std::vector;

struct State {
public:
  vector<double> s0s;
  // Data data;
  State(vector<double> s0s) : s0s(s0s) {}
};
