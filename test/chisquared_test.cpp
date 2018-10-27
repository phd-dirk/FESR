#include "../src/configuration.hpp"
#include "../src/weights.hpp"
#include "../src/chisquared.hpp"
#include <gtest/gtest.h>

TEST(Chisquared_Test, chi2) {
  const Configuration config("./test/configuration_test2.json");
  const Chi2 chi2(config);
  EXPECT_NEAR(
    chi2(
      0.3179, 0.021, -0.15, 0.24,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5
    ), 10.520253962854611, 1e-15
  );

  // std::vector<Input> inputs = {
  //   {
  //     Weight(1),
  //     {
  //       3.15723, 3.0, 2.800, 2.600, 2.400,
  //       2.300, 2.200, 2.100, 2.000
  //     }
  //   }
  // };

  // EXPECT_NEAR(
  //   Chi2::chi2(
  //     0.3179, 0.021, -0.15, 0.24,
  //     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  //     5, pow(1.77686, 2), 0.13957018, 92.21e-3,
  //     2.2e-3, 1.3, 0.4, 0.19e-3, 1.8, 0.21,
  //     3, 3, { 2.8e-3, 5.0e-3, 97.0e-3 }, 0.9742, 1.0198,
  //     Condensates(
  //       0.3156,
  //       { -pow(0.272, 3), -pow(0.272, 3), 0.8*(-pow(0.272, 3)) },
  //       { 2.8e-3, 5.0e-3, 97.0e-3 }
  //     ),
  //     inputs,
  //     { "FO", true, true, true, false, true },
  //     ExpMoms(
  //       "/Users/knowledge/Developer/PhD/FESR/aleph.json",
  //       inputs,
  //       pow(0.177686, 2), // sTau
  //       17.815, // be
  //       0.023, // dBe
  //       0.97420, // Vud
  //       0.00021, // dVud
  //       1.0198, // SEW
  //       0.0006, // dSEW
  //       92.21e-3, // fPi
  //       0.14e-3, // dFPi
  //       0.13957018, // mPiM
  //       0.99743669 // RVANormalization
  //     )
  //   ), 564.01914479368497, 1e-15
  // );
}
