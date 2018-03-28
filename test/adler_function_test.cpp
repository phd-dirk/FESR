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

TEST_F(AdlerFunctionTest, D2) {
  const complex<double> s(3.1570893123374919, -1.9866720523795795e-5);
  const complex<double> mu2(3.1570893124, 0.);
  const double astau = 0.31927;
  const int order = 2;
  const uint r = 1;
  EXPECT_NEAR(adler->D2(s, mu2, astau, order, r).real(), 3.4324915657182825e-7, 1e-09);
}

TEST_F(AdlerFunctionTest, D4) {
  const complex<double> s(3.1570893123374919, 1.9866720523795795e-5);
  const complex<double> mu2(3.1570893124, 0.);
  const double astau = 0.31927;
  const double aGGinv = 2.1e-2;
  const int order = 5;
  EXPECT_NEAR(adler->D4(s, mu2, astau, aGGinv, order, 1).real(), 2.7591458364939887e-4, 1e-13);
}
TEST_F(AdlerFunctionTest, D4CInt) {
  const double s0 = 3.1570893124;
  const double astau = 0.31927;
  const double aGGinv = 2.1e-2;
  const uint r = 1;
  const int order = 5;
  EXPECT_NEAR(adler->D4CInt(s0, wD00, astau, aGGinv, order, r), 1.6944347548019322e-3, 1e-11);
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
