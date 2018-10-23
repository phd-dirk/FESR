#include "../src/pseudoscalar_pheno.hpp"
#include <gtest/gtest.h>

TEST(PSPheno_Test, deltaP) {
  // deltaP(s0, weight)
  EXPECT_NEAR(
    PSPheno::deltaP(
      3.0, Weight(1),
      pow(1.77682, 2), 0.13957018, 92.21e-3,
      2.2e-3, 1.3, 0.4, 0.19e-3, 1.8, 0.21
    ), -2.7752063630383842e-3, 1e-15
  );

  EXPECT_NEAR(
    PSPheno::deltaP(
      2.0, Weight(1),
      pow(1.77682, 2), 0.13957018, 92.21e-3,
      2.2e-3, 1.3, 0.4, 0.19e-3, 1.8, 0.21
    ), -4.138431838090359e-3, 1e-15
  );

  EXPECT_NEAR(
    PSPheno::deltaP(
      2.0, Weight(6),
      pow(1.77682, 2), 0.13957018, 92.21e-3,
      2.2e-3, 1.3, 0.4, 0.19e-3, 1.8, 0.21
    ), -4.809780763096577e-5, 1e-15
  );

  EXPECT_NEAR(
    PSPheno::deltaP(
      3.0, Weight(1),
      pow(1.77686, 2), 0.13957018, 92.21e-3,
      2.2e-3, 1.3, 0.4, 0.19e-3, 1.8, 0.21
    ), -2.7750836295103356e-3, 1e-15
  );
}
