#include "../src/configuration.hpp"
#include "../src/weights.hpp"
#include "../src/chisquared.hpp"
#include <gtest/gtest.h>

TEST(Chisquared_Test, chi2) {
  const Configuration config("./test/configuration_test2.json");
  const Chi2 chi2(config);
  EXPECT_NEAR(
    chi2(
      0.3179, 0.021, -0.15, 0.24, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5
    ), 10.520253962854611, 1e-8
  );

  const Chi2 chi2_2(
    3, 3, 5,
    pow(1.77686, 2),
    0.13957018,
    92.21e-3,
    0.14e-3,
    2.2e-3, 1.3, 0.4, 0.19e-3, 1.8, 0.21,
    { 2.8e-3, 5.0e-3, 97.0e-3 },
    17.815,
    0.023,
    0.97420,
    0.00021,
    1.0198,
    0.0006,
    0.3156,
    { -pow(0.272, 3), -pow(0.272, 3), 0.8*-pow(0.272, 3) },
    {
      {
        Weight(1),
        {
          pow(1.77686, 2)
        }
      }
    },
    { "FO", true, true, true, false, false, false, true },
    0.99743669
  );

  EXPECT_NEAR(
    chi2_2(
      0.3179, 0.021, -0.15, 0.24, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5
    ), 2.2405957912637615, 1e-12
  );

}
