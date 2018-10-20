#include "../src/types.hpp"
#include <gtest/gtest.h>
#include "../src/alpha_s.hpp"
#include <cmath>

class AlphaSTest : public ::testing::Test {
};

TEST_F(AlphaSTest, complex) {
  // BETA 5th order --------------------------------------------------
  cmplx q(3.0, 0.0);
  EXPECT_NEAR(AlphaS::run(q, 3.1, 0.1028).real(), 0.10379258455665125, 1e-14);

  // q = cmplx(3.0, 1.5);
  // EXPECT_NEAR((*amu_)(q, 3.1, 0.1028).real(), 0.09858208182028837, 1e-14);
  // EXPECT_NEAR((*amu_)(q, 3.1, 0.1028).imag(), -0.01286649420615981, 1e-14);

  // EXPECT_NEAR((*amu_)(cmplx(3.0, 2.0), 3.1, 0.1028).real(), 0.09567038211680548, 1e-14);
  // EXPECT_NEAR((*amu_)(cmplx(3.0, 2.0), 3.1, 0.1028).imag(), -0.015397224860415053, 1e-14);
};
