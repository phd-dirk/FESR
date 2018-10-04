#include "../src/types.hpp"
#include "../src/configuration.hpp"
#include "../src/weights.hpp"
#include "../src/chisquared.hpp"
#include <gtest/gtest.h>

class ChisquaredTest : public ::testing::Test {
 protected:
  Chisquared *chi_;
  Configuration *config;

  virtual void SetUp() {
    const uint order = 5;
    const uint nc = 3, nf = 3;

    config = new Configuration("./test/configuration_test.json");
    chi_ = new Chisquared(*config);
  }

  virtual void TearDown() {
    delete chi_;
  }
};

TEST_F(ChisquaredTest, inverseCovarianceMatrix) {
  mat invCov = chi_->invCovMat_;
  EXPECT_NEAR(invCov(0, 0), 19290.123456790123 , 1e-15);
  EXPECT_NEAR(invCov(1, 1), 310415.88885578304 , 1e-7);
  EXPECT_NEAR(invCov(4, 4), 2392765.2394170612, 1e-6);
}

TEST_F(ChisquaredTest, chi2) {
  Chisquared chi = *chi_;
  EXPECT_NEAR(chi(0.31927, 0.021, -0.1894, 0.16315, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0), 31701.49710037827, 1e-5);
  // EXPECT_NEAR(chi({ 3.0, 2.0 }, 0.31927, 0.021, -0.1894, 0.16315), 1.4619737304634208, 1e-11);
}
