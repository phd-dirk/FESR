#include "../src/theoretical_moments.hpp"
#include "../src/weights.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <complex>

class TheoreticalMomentsTest : public ::testing::Test {
protected:
  TheoreticalMoments *thMom_;
  Weight *weight_;
  Configuration *config;
  const uint order_ = 4;
  const double sTau_ = 3.1570893124;
  const double astau_ = 0.31927;
  const double aGGinv_ = 2.1e-2;
  const double rhoVpA_ = -0.1894;
  const double c8VpA_ = 0.16315;
  virtual void SetUp() {
    config = new Configuration("./test/configuration_test.json");
    weight_ = new Weight(1);
    thMom_ = new TheoreticalMoments(*config);
  }

  virtual void TearDown() {
    delete thMom_;
  }
};

TEST_F(TheoreticalMomentsTest, IntegralMoment) {
  // thMom(const int &i, const double &astau, const double &aGGinv,
  // const double &rhoVpA, const double &c8VpA, const double &order)
  TheoreticalMoments th = *thMom_;
  Weight w(1);
  EXPECT_NEAR(
      th.thMom(
          3.1570893124, w, astau_, aGGinv_, rhoVpA_, c8VpA_, 5,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      ),
      3.4634999665533375, config->tolerance
  );
  EXPECT_NEAR(th.thMom(3.0, w, 0.31927, 0.021, -0.1894, -0.161315, 5,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
              3.4848911406246970, config->tolerance);
  EXPECT_NEAR(th.thMom(3.0, w, 0.2, 0.021, -0.1894, -0.161315, 5,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
              3.1542708728473938, config->tolerance);
  EXPECT_NEAR(th.thMom(3.0, w, 0.2, 0.1, -0.1894, -0.161315, 5,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
              3.1571219772351791, config->tolerance);
  EXPECT_NEAR(th.thMom(3.0, w, 0.2, 0.1, -0.3, 0.9, 5,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
              3.1129921630758908, config->tolerance);
}

TEST_F(TheoreticalMomentsTest, ThMomDV) {
  std::vector<Input> inputs = {{
      Weight(1),
      {
        3.15709, 3.0, 2.800, 2.600, 2.400,
        2.300, 2.200, 2.100, 2.000
      }
    }};
  ThMomContribs thMomContribs = { "FO", false, false, false, true, false };
  TheoreticalMoments th = TheoreticalMoments(
    Configuration(
      17.815, //be
      0.023, //dBe
      0.97425, // Vud
      0.0022, // dVud
      inputs,
      1.77682, // mTau
      3, // nc
      3, // nf
      5, // order
      0.99743669, // RVANormalization
      thMomContribs
    ));

  EXPECT_NEAR(
      th.thMom(
        3.15709, Weight(1),
        0.31927, 0.021, -0.1894, 0.16315, 5,
        0.0, 1.0, 1.0, 1.0,
        0.0, 1.0, 1.0, 1.0
      ),
      -1.749, 1e-15
  );
}

TEST_F(TheoreticalMomentsTest, Delta0) {
  const double astau = 0.31927;
  const int order = 5;
  EXPECT_NEAR(thMom_->del0(sTau_, *weight_, sTau_, astau, order), 0.20298958142552484, 1e-13);
}

TEST_F(TheoreticalMomentsTest, Delta4) {
  EXPECT_NEAR(thMom_->del4(sTau_, *weight_, sTau_, astau_, aGGinv_), 7.93483938348380651e-4, 1.e-11);
  EXPECT_NEAR(thMom_->del4(3.1570893124000001, Weight(1), 3.1570893124000001,  0.32307, 0.021), 8.14191092757126224e-4, 1.e-6);
}

TEST_F(TheoreticalMomentsTest, Delta68) {
  const double rho = -0.1893979224795759;
  const double c8 = 0.16314594513667133;
  // Test delta_V+A^(8)
  EXPECT_NEAR(thMom_->del68(sTau_, *weight_, 0., c8), -1.2966374009228992e-3, 1e-13);
  // Test delta_V+A^(6)
  EXPECT_NEAR(thMom_->del68(sTau_, *weight_, rho, 0.), -7.1284580508113966e-3, 1e-13);
}

TEST_F(TheoreticalMomentsTest, DeltaP) {
  EXPECT_NEAR(thMom_->deltaP(sTau_, *weight_), -2.63897241291510083e-3, 1e-13);
}
