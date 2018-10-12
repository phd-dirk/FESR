#include "../src/configuration.hpp"
#include "../src/weights.hpp"
#include "../src/chisquared.hpp"
#include <gtest/gtest.h>

class ChisquaredTest : public ::testing::Test {
 protected:
  Chisquared *chi_;

  virtual void SetUp() {
    chi_ = new Chisquared(Configuration("./test/configuration_test.json"));
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


// TEST_F(ChisquaredTest, chi2) {
//   Chisquared chi = *chi_;
//   EXPECT_NEAR(
//     chi(
//       0.32307, 0.021, -0.30949, -0.030869,
//       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
//     )
//     , 31701.49710037827, 1e-5
//   );
// }

// TEST_F(ChisquaredTest, chi2DV) {
//   std::vector<Input> inputs = {{
//       Weight(1),
//       {
//         3.15709, 3.0, 2.800, 2.600, 2.400,
//         2.300, 2.200, 2.100, 2.000
//       }
//     }};
//   ThMomContribs thMomContribs = { "FO", false, false, false, true, false };
//   Chisquared chi = Chisquared(
//     Configuration(
//       17.815, //be
//       0.023, //dBe
//       inputs,
//       1.77682, // mTau
//       3, // nc
//       3, //nf
//       5, // order
//       0.99743669, // RVANormalization
//       thMomContribs
//    ));

//   EXPECT_NEAR(
//     chi(
//       0.32307, 0.021, -0.30949, -0.030869,
//       3.56, 0.58, -1.92, 4.07,
//       1.68, 1.41, 5.16, 2.13
//     )
//     , 397332.77146924846, 1e-5
//   );
// }
