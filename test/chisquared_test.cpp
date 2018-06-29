#include "../src/types.hpp"
#include "../src/configuration.hpp"
#include "../src/weights.hpp"
#include "../src/chisquared.hpp"
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
    std::ifstream configFile("./test/configuration_test2.json");
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
  Chisquared chi = *chi_;
  EXPECT_NEAR(chi({ 3.0 }, 0.31927, 0.021, -0.1894, 0.16315), 0.37595084146692792 , 1e-12);
  EXPECT_NEAR(chi({ 3.0, 2.0 }, 0.31927, 0.021, -0.1894, 0.16315), 1.4619737304634208, 1e-11);
}
