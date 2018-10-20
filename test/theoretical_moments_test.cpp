#include "../src/theoretical_moments.hpp"
#include "../src/weights.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <complex>

class TheoreticalMomentsTest : public ::testing::Test {
protected:
  TheoreticalMoments *thMom_;
  const uint order_ = 4;
  const double sTau_ = 3.1570893124;
  const double astau_ = 0.31927;
  const double aGGinv_ = 2.1e-2;
  const double rhoVpA_ = -0.1894;
  const double c8VpA_ = 0.16315;
  virtual void SetUp() {

    std::vector<Input> inputs = {
      { Weight(1), { 3.0, 2.0, 1.0 } },
      { Weight(6), { 2.8, 1.8 } },
    };
    ThMomContribs thMomContribs = { "FO", true, true, true, false, true };

    thMom_ = new TheoreticalMoments(
      Configuration(
        0.3156,
        17.815,
        0.023,
        0.97425,
        0.00022,
        1.0198,
        0.13957018,
        92.21e-3,
        2.2e-3,
        1.3,
        0.4,
        0.19e-3,
        1.8,
        0.21,
        { 2.8e-3, 5.0e-3, 97.0e-3 },
        { 0.0201236, 0.021236, 0.0160989 },
        inputs,
        1.77682,
        3,
        3,
        5,
        0.99743669,
        thMomContribs
      )
    );
  }

  virtual void TearDown() {
    delete thMom_;
  }
};

// TEST_F(TheoreticalMomentsTest, ThMom) {
//   TheoreticalMoments thMoms = TheoreticalMoments(
//     3,
//     3,
//     { 2.8e-3, 5e-3, 97e-3 },
//     Condensates(
//       0.3156,
//       { 0.0201236, 0.0201236, 0.0160989 },
//       { 2.8e-3, 5.0e-3, 97.0e-3 }
//     ),
//     3.1572314596,
//     0.13957018,
//     92.21e-3,
//     2.2e-3,
//     1.3,
//     0.4,
//     0.19e-3,
//     1.8,
//     0.21,
//     {
//       {
//         Weight(1),
//         {
//           3.1572314596000002, 3.0, 2.800, 2.600, 2.400, 2.300,
//             2.200, 2.100, 2.000
//         }
//       }
//     },
//     { "FO", true, true, true, false, true },
//     0.97425,
//     1.0198
//   );

//   EXPECT_NEAR(
//       thMoms.thMom(
//           3.0, Weight(1), , aGGinv_, rhoVpA_, c8VpA_, 5,
//           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
//       ),
//       3.4634999665533375, 1e-14
//   );

// }

// TEST_F(TheoreticalMomentsTest, IntegralMoment) {
//   // thMom(const int &i, const double &astau, const double &aGGinv,
//   // const double &rhoVpA, const double &c8VpA, const double &order)
//   TheoreticalMoments th = *thMom_;
//   EXPECT_NEAR(
//       th.thMom(
//           3.1570893124, Weight(1), astau_, aGGinv_, rhoVpA_, c8VpA_, 5,
//           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
//       ),
//       3.4634999665533375, 1e-14
//   );
//   EXPECT_NEAR(
//     th.thMom(3.0, Weight(1), 0.31927, 0.021, -0.1894, -0.161315, 5,
//              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
//     3.4848911406246970, 1e-14
//   );
//   EXPECT_NEAR(
//     th.thMom(
//       3.0, Weight(1), 0.2, 0.021, -0.1894, -0.161315, 5,
//       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
//     ),
//     3.1542708728473938, 1e-14
//   );
//   EXPECT_NEAR(
//     th.thMom(
//       3.0, Weight(1), 0.2, 0.1, -0.1894, -0.161315, 5,
//       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
//     ),
//     3.1571219772351791, 1e-14
//   );
//   EXPECT_NEAR(
//     th.thMom(
//       3.0, Weight(1), 0.2, 0.1, -0.3, 0.9, 5,
//       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
//     ),
//     3.1129921630758908, 1e-14
//   );
// }

// TEST_F(TheoreticalMomentsTest, ThMomDV) {
//   std::vector<Input> inputs = {{
//       Weight(1),
//       {
//         3.15709, 3.0, 2.800, 2.600, 2.400,
//         2.300, 2.200, 2.100, 2.000
//       }
//     }};
//   ThMomContribs thMomContribs = { "FO", false, false, false, true, false };
//   TheoreticalMoments th = TheoreticalMoments(
//     Configuration(
//       17.815, //be
//       0.023, //dBe
//       0.97425, // Vud
//       0.0022, // dVud
//       inputs,
//       1.77682, // mTau
//       3, // nc
//       3, // nf
//       5, // order
//       0.99743669, // RVANormalization
//       thMomContribs
//     ));

//   EXPECT_NEAR(
//       th.thMom(
//         3.15709, Weight(1),
//         0.31927, 0.021, -0.1894, 0.16315, 5,
//         0.0, 1.0, 1.0, 1.0,
//         0.0, 1.0, 1.0, 1.0
//       ),
//       -1.749, 1e-15
//   );
// }

TEST_F(TheoreticalMomentsTest, Delta0) {
  const double astau = 0.31927;
  const int order = 5;
  EXPECT_NEAR(thMom_->del0(sTau_, Weight(1), sTau_, astau, order), 0.20298958142552484, 1e-14);
}

// TEST_F(TheoreticalMomentsTest, Delta4) {
//   EXPECT_NEAR(thMom_->del4(sTau_, Weight(1), sTau_, astau_, aGGinv_), 7.93483938348380651e-4, 1.e-14);
//   EXPECT_NEAR(thMom_->del4(3.0, Weight(1), 3.0,  0.32307, 0.021), 2.41777114013837419e-4, 1.e-14);
// }

TEST_F(TheoreticalMomentsTest, Delta68) {
  const double rho = -0.1893979224795759;
  const double c8 = 0.16314594513667133;
  // Test delta_V+A^(8)
  EXPECT_NEAR(thMom_->del68(sTau_, Weight(1), 0., c8), -1.2966374009228992e-3, 1e-14);
  // Test delta_V+A^(6)
  EXPECT_NEAR(thMom_->del68(sTau_, Weight(1), rho, 0.), -7.1284580508113966e-3, 1e-14);
}
// TEST_F(TheoreticalMomentsTest, DeltaP) {
//   EXPECT_NEAR(thMom_->deltaP(sTau_, Weight(1)), -2.63897241291510083e-3, 1e-14);
// }
