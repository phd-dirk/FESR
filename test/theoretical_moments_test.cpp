#include "../src/constants.hpp"
#include "../src/theoretical_moments.hpp"
#include "../src/weights.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <complex>

using std::pow;
using std::complex;

class AdlerFunctionTest : public ::testing::Test {
 protected:
  AdlerFunction *adler;
  Constants *const_;
  virtual void SetUp() {
    const_ = new Constants(3, 3);
    adler = new AdlerFunction(5, *const_);
  }

  virtual void TearDown() {
    delete adler;
    delete const_;
  }
};


TEST_F(AdlerFunctionTest, D0) {
  complex<double> s(3., 3.);
  complex<double> mu2(3.0, 0);
  int order = 3;
  EXPECT_NEAR(adler->D0(s, mu2, const_->kAlphaTau, order).real(), 2.6878635987293748e-2, 1e-13);
  EXPECT_NEAR(adler->D0(s, mu2, const_->kAlphaTau, order).imag(), 1.5454361164073294e-3, 1e-13);

  s = complex<double>(7.0, 2.0);
  mu2 = complex<double>(1.5, 2.2);
  EXPECT_NEAR(adler->D0(s, mu2, const_->kAlphaTau, order).real(), 2.6732096654905575e-2, 1e-13);
  EXPECT_NEAR(adler->D0(s, mu2, const_->kAlphaTau, order).imag(), -1.0995430321549232e-3, 1e-13);
}

TEST_F(AdlerFunctionTest, CIntD0) {
  double s0 = 3.;
  int order = 3;
  EXPECT_NEAR(adler->D0CInt(s0, wD00, const_->kAlphaTau, order), 1.8129891021138491, 1e-13);
  // EXPECT_NEAR(adler->contourIntegral(s0, wD00).imag(), 0., 1e-13);

  s0 = 2.6;
  EXPECT_NEAR(adler->D0CInt(s0, wD00, const_->kAlphaTau, order), 1.8396044152966260, 1e-13);
  // EXPECT_NEAR(adler->contourIntegral(s0, wD00).imag(), 0., 1e-13);
}

TEST_F(AdlerFunctionTest, DeltaP) {
  
}
