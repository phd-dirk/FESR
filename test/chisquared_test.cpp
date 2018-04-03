#include "../src/constants.hpp"
#include "../src/weights.hpp"
#include "../src/chisquared.hpp"
#include "../src/s0_sets.hpp"
#include "json.hpp"
#include <gtest/gtest.h>
#include <complex>
#include <fstream>
using json = nlohmann::json;

class ChisquaredTest : public ::testing::Test {
 protected:
  Chisquared *chi_;
  Constants *const_;

  virtual void SetUp() {
    const uint order = 5;
    const uint nc = 3, nf = 3;
    std::ifstream configFile("./test/configuration_test.json");
    json config;
    configFile >> config;

    const_ = new Constants(config);
    Weight weight(1);
    chi_ = new Chisquared(order, s0Set, weight, config, *const_);
  }

  virtual void TearDown() {
    delete chi_;
    delete const_;
  }
};

TEST_F(ChisquaredTest, chi2) {
  double xx[4] = { 0.31921, 0.21e-1, -0.1894, 0.16315 };
  Chisquared chi = *chi_;
  // EXPECT_NEAR(chi(xx), 9.3758915658581827, 1e-13);
}
