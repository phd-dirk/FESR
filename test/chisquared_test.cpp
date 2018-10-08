#include "../src/types.hpp"
#include "../src/configuration.hpp"
#include "../src/weights.hpp"
#include "../src/chisquared.hpp"
#include <gtest/gtest.h>

class ChisquaredTest : public ::testing::Test {
 protected:
  Chisquared *chi_;
  Chisquared *chi2_;
  Chisquared *chi3_;

  virtual void SetUp() {
    chi_ = new Chisquared(Configuration("./test/configuration_test.json"));

    std::vector<Input> inputs = {{
        Weight(1),
        {
          3.15709, 3.0, 2.800, 2.600, 2.400,
            2.300, 2.200, 2.100, 2.000
        }
      }};
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

    inputs = {{
        Weight(6),
        {
          3.15709, 3.0, 2.800, 2.600, 2.400,
          2.300, 2.200, 2.100, 2.000
        }
      }};
    chi3_ = new Chisquared(
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
      0.32307, 0.021, -0.30949, -0.030869,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )
    , 31701.49710037827, 1e-5
  );

  EXPECT_NEAR(
    (*chi2_)(
      0.31927, 0.021, -0.1894, 0.16315,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    ), 9.484259311895366, 1e-5
  );

  EXPECT_NEAR(
    (*chi3_)(
      0.31927, 0.021, -0.1894, 0.16315,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    ), 18603.923830229789, 1e-5
  );
}
