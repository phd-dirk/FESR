#include "../src/types.hpp"
#include "../src/configuration.hpp"
#include "../src/weights.hpp"
#include "../src/chisquared.hpp"
#include "../src/s0_sets.hpp"
#include "json.hpp"
#include <gtest/gtest.h>
#include <fstream>
using json = nlohmann::json;

class ChisquaredTest : public ::testing::Test {
 protected:
  Chisquared *chi_;

  virtual void SetUp() {
    const uint order = 5;
    const uint nc = 3, nf = 3;
    std::ifstream configFile("./test/configuration_test.json");
    json jsonConfig;
    configFile >> jsonConfig;

    Configuration config(jsonConfig);

    chi_ = new Chisquared(config);
  }

  virtual void TearDown() {
    delete chi_;
  }
};

TEST_F(ChisquaredTest, chi2) {
  double xx[4] = { 0.31921, 0.21e-1, -0.1894, 0.16315 };
  Chisquared chi = *chi_;
  // EXPECT_NEAR(chi(xx), 9.3758915658581827, 1e-13);
}
