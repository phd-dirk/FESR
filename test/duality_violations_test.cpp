#include "../src/duality_violations.hpp"
#include "../src/weights.hpp"
#include <gtest/gtest.h>

#include <complex>

TEST(DVTest, cIntVA) {
  EXPECT_NEAR(DV::cIntVA(3.0, Weight(1), 1.1, 1.2, 1.3, 1.4), -4.3482045718910615E-002, 1e-15);
}

TEST(DVTest, intP0) {
  EXPECT_NEAR(DV::intP0(3.0, 1.1, 1.2, 1.3, 1.4), -1.5364831049097595e-2, 1e-15);
}

TEST(DVTest, intP1) {
  EXPECT_NEAR(DV::intP1(3.0, 1.1, 1.2, 1.3, 1.4), -4.3816982693323563e-2, 1e-15);
}

TEST(DVTest, intP2) {
  EXPECT_NEAR(DV::intP2(3.0, 1.1, 1.2, 1.3, 1.4), -8.4654144846656038e-2, 1e-15);
}

TEST(DVTest, intP3) {
  EXPECT_NEAR(DV::intP3(3.0, 1.1, 1.2, 1.3, 1.4), -0.14103982460489065, 1e-15);
}
