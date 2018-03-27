#include "../src/constants.hpp"
#include "../src/theoretical_moments.hpp"
#include "../src/weights.hpp"
#include "../src/s0_sets.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <complex>

using std::pow;
using std::complex;

class TheoreticalMomentsTest : public ::testing::Test {
protected:
  TheoreticalMoments *thMom_;
  Constants *const_;
  virtual void SetUp() {
    int order = 5;
    const_ = new Constants(3, 3);
    thMom_ = new TheoreticalMoments(order, s0Set, wD00, *const_);
  }

  virtual void TearDown() {
    delete thMom_;
    delete const_;
  }
};

TEST_F(TheoreticalMomentsTest, Delta0) {
  const double s0 = const_->kSTau;
  const double astau = 0.31927;
  const int order = 5;
  EXPECT_NEAR(thMom_->del0(s0, wD00, astau, order), 0.20298958142552484, 1e-13);
}

TEST_F(TheoreticalMomentsTest, Delta2) {
  const double s0 = const_->kSTau;
  const double astau = 0.31927;
  const int order = 2;
  // EXPECT_NEAR(thMom_->del2(s0, wD00, astau, order), 3.90335378985309371e-5, 1e-13);
}

TEST_F(TheoreticalMomentsTest, Delta4) {
  const double s0 = const_->kSTau;
  const double aGGinv = 2.1e-2;
  const double astau = 0.31927;
  EXPECT_NEAR(thMom_->del4(s0, wD00, astau, aGGinv), 7.93483938348380651e-4, 1e-11);
}

TEST_F(TheoreticalMomentsTest, Delta68) {
  const double s0 = const_->kSTau;
  const double rho = -0.1893979224795759;
  const double c8 = 0.16314594513667133;
  // Test delta_V+A^(8)
  EXPECT_NEAR(thMom_->del68(s0, wD00, 0., c8), -1.2966374009228992e-3, 1e-13);
  // Test delta_V+A^(6)
  EXPECT_NEAR(thMom_->del68(s0, wD00, rho, 0.), -7.1284580508113966e-3, 1e-13);
}

TEST_F(TheoreticalMomentsTest, DeltaP) {
  const double s0 = const_->kSTau;
  EXPECT_NEAR(thMom_->deltaP(const_->kSTau, wR00), -2.63897241291510083e-3, 1e-13);
}

