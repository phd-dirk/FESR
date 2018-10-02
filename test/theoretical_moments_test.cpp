#include "../src/theoretical_moments.hpp"
#include "../src/weights.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <complex>

class TheoreticalMomentsTest : public ::testing::Test {
protected:
  TheoreticalMoments *thMom_;
  Weight *weight_;
  const uint order_ = 4;
  const double sTau_ = 3.1570893124;
  const double astau_ = 0.31927;
  const double aGGinv_ = 2.1e-2;
  const double rhoVpA_ = -0.1894;
  const double c8VpA_ = 0.16315;
  virtual void SetUp() {
    Configuration config("./test/configuration_test.json");
    weight_ = new Weight(1);
    thMom_ = new TheoreticalMoments(config);
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
  EXPECT_NEAR(th.thMom(3.1570893124, w, astau_, aGGinv_, rhoVpA_, c8VpA_, 5,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
              3.4634999665533375, 1.e-14);
  EXPECT_NEAR(th.thMom(3.0, w, 0.31927, 0.021, -0.1894, -0.161315, 5,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
              3.4848911406246970, 1.e-14);
  EXPECT_NEAR(th.thMom(3.0, w, 0.2, 0.021, -0.1894, -0.161315, 5,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
              3.1542708728473938, 1.e-14);
  EXPECT_NEAR(th.thMom(3.0, w, 0.2, 0.1, -0.1894, -0.161315, 5,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
              3.1571219772351791, 1.e-14);
  EXPECT_NEAR(th.thMom(3.0, w, 0.2, 0.1, -0.3, 0.9, 5,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
              3.1129921630758908, 1.e-14);
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

