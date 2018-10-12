#include "../src/duality_violations.hpp"
#include "../src/weights.hpp"
#include <gtest/gtest.h>

#include <complex>

class DualityViolationTest : public ::testing::Test {
};

// TEST_F(DualityViolationTest, DVMomentVpA) {
//   DualityViolations dv;
//   EXPECT_NEAR(
//     dv.DVMomentVpA(
//       3.15709, Weight(1),
//       0.0, 1.0, 1.0, 1.0,
//       0.0, 1.0, 1.0, 1.0
//     )
//     , -0.40033129402047263, 1e-15
//   );
// }

// TEST_F(DualityViolationTest, cintDV) {
//   DualityViolations dv;
//   EXPECT_NEAR(dv.cintDVp_VA(1.0, Weight(1), 0.0, 1.0, 1.0, 1.0), -1.73494103008802, 1e-15);
//   EXPECT_NEAR(dv.cintDVp_VA(3.15709, Weight(1), -log(3), 3.7, 4.2, 7.1), -0.00572861, 1e-15);
// }

TEST_F(DualityViolationTest, intDVp0) {
  DualityViolations dv;
  EXPECT_NEAR(dv.intDVp0(1.0, 0.0, 1.0, 1.0, 1.0), 0.0907099817825180, 1e-15);
  EXPECT_NEAR(dv.intDVp0(3.15709, 0.0, 3.7, 4.2, 7.1), 5.68364e-7, 1e-13);
}

TEST_F(DualityViolationTest, intDVp1) {
  DualityViolations dv;
  EXPECT_NEAR(dv.intDVp1(1.0, 0.0, 1.0, 1.0, 1.0), 0.0141640489454048, 1e-15);
  EXPECT_NEAR(dv.intDVp1(3.15709, -log(3), 3.7, 4.2, 7.1), 5.18599e-6, 1e-10);
}

TEST_F(DualityViolationTest, intDVp2) {
  DualityViolations dv;
  EXPECT_NEAR(dv.intDVp2(1.0, 0.0, 1.0, 1.0, 1.0), -0.306183731348453, 1e-15);
  EXPECT_NEAR(dv.intDVp2(3.15709, -log(3), 3.7, 4.2, 7.1), 0.000015651, 1e-9);
}

TEST_F(DualityViolationTest, intDVp3) {
  DualityViolations dv;
  EXPECT_NEAR(dv.intDVp3(1.0, 0.0, 1.0, 1.0, 1.0), -1.37210110295795, 1e-14);
  EXPECT_NEAR(dv.intDVp3(3.15709, 0.0, 3.7, 4.2, 7.1), 0.0000156031, 1e-11);
}
