#include "../src/types.hpp"
#include "../src/configuration.hpp"
#include "../src/weights.hpp"
#include "../src/chisquared.hpp"
#include <gtest/gtest.h>

class ChisquaredTest : public ::testing::Test {
 protected:
  Chisquared *chi_;
  Chisquared *chi2_;
  Configuration *config;

  virtual void SetUp() {
    config = new Configuration("./test/configuration_test.json");
    chi_ = new Chisquared(*config);

    std::vector<Input> inputs;
    Weight w(1);
    std::vector<double> s0s = {
      3.15709, 3.0, 2.800, 2.600, 2.400, 2.300, 2.200, 2.100, 2.000
    };
    inputs.push_back({ w, s0s });
    chi2_ = new Chisquared(
      Configuration(
        17.815, //be
        0.023, //dBe
        inputs,
        1.77682, // mTau
        3, // nc
        3, //nf
        5, // order
        0.99743669, // RVANormalization
        "FO" // scheme
     ));
  }

  virtual void TearDown() {
    delete chi_;
    delete chi2_;
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
  EXPECT_NEAR(
    chi(
      config->inputs_, 0.32307, 0.021, -0.30949, -0.030869,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )
    , 31701.49710037827, 1e-5
  );

  std::vector<Input> inputs;
  Weight w(1);
  std::vector<double> s0s = {
    3.15709, 3.0, 2.800, 2.600, 2.400, 2.300, 2.200, 2.100, 2.000
  };
  inputs.push_back({ w, s0s });
  EXPECT_NEAR(
    (*chi2_)(
      inputs, 0.31927, 0.021, -0.1894, 0.16315,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    ), 1.4619737304634208, 1e-11
  );
}
