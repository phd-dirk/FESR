#include "../src/constants.hpp"
#include "../src/theoretical_moments.hpp"
#include "../src/weights.hpp"
#include "../src/s0_sets.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <complex>

using std::pow;
using std::complex;

class AdlerFunctionTest : public ::testing::Test {
 protected:
  TheoreticalMoments *thMom_;
  AdlerFunction *adler;
  Constants *const_;
  virtual void SetUp() {
    int order = 5;
    const_ = new Constants(3, 3);
    adler = new AdlerFunction(5, *const_);
    thMom_ = new TheoreticalMoments(order, s0Set, wD00, *const_);
  }

  virtual void TearDown() {
    delete adler;
    delete const_;
    delete thMom_;
  }
};

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

TEST_F(AdlerFunctionTest, D0) {
  complex<double> s(3., 3.);
  complex<double> mu2(3.0, 0);
  int order = 5;
  EXPECT_NEAR(adler->D0(s, mu2, const_->kAlphaTau, order).real(), 2.6878635987293748e-2, 1e-13);
  EXPECT_NEAR(adler->D0(s, mu2, const_->kAlphaTau, order).imag(), 1.5454361164073294e-3, 1e-13);

  s = complex<double>(7.0, 2.0);
  mu2 = complex<double>(1.5, 2.2);
  EXPECT_NEAR(adler->D0(s, mu2, const_->kAlphaTau, order).real(), 2.6732096654905575e-2, 1e-13);
  EXPECT_NEAR(adler->D0(s, mu2, const_->kAlphaTau, order).imag(), -1.0995430321549232e-3, 1e-13);
}

TEST_F(AdlerFunctionTest, CIntD0) {
  double s0 = 3.;
  int order = 5;
  EXPECT_NEAR(adler->D0CInt(s0, wD00, const_->kAlphaTau, order), 1.8129891021138491, 1e-13);
  // EXPECT_NEAR(adler->contourIntegral(s0, wD00).imag(), 0., 1e-13);

  s0 = 2.6;
  EXPECT_NEAR(adler->D0CInt(s0, wD00, const_->kAlphaTau, order), 1.8396044152966260, 1e-13);
  // EXPECT_NEAR(adler->contourIntegral(s0, wD00).imag(), 0., 1e-13);
}

TEST_F(AdlerFunctionTest, D68) {
  const complex<double> s(-1.7585656673884618, 1.9150579086499120);
  const double rho = -0.1894;
  const double c8 = 0.16315;
  EXPECT_NEAR(adler->D68(s, rho, c8).real(), -3.9659215963967076e-4, 1e-13);
  EXPECT_NEAR(adler->D68(s, rho, c8).imag(), 1.7341429761413442e-4, 1e-13);
}

TEST_F(AdlerFunctionTest, D68CInt) {
  const double s0 = 2.4;
  const double rho = -0.1894;
  const double c8 = 0.16315;
  EXPECT_NEAR(adler->D68CInt(s0, wD00, rho, c8), -6.0327814112243174e-2, 1e-13);
}

TEST_F(TheoreticalMomentsTest, Delta0) {
  const double s0 = const_->kSTau;
  const double astau = 0.31927;
  const int order = 5;
  EXPECT_NEAR(thMom_->del0(s0, wD00, astau, order), 0.20298958142552484, 1e-13);
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

