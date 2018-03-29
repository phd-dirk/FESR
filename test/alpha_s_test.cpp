#include <gtest/gtest.h>
#include "../src/constants.hpp"
#include "../src/alpha_s.hpp"
#include <complex>
#include <cmath>

using std::complex;
using std::sqrt;


class AlphaSTest : public ::testing::Test {
 protected:
  Constants *const_;
  AlphaS *amu_;
  const double astau_ = 0.31927;

  virtual void SetUp() {
    const int nc = 3;
    const int nf = 3;
    const_ = new Constants(nc, nf);
    amu_ = new AlphaS(*const_);
  }
};

TEST_F(AlphaSTest, complex) {
  complex<double> q2(3., 3.);
  EXPECT_NEAR((*amu_)(q2, const_->kSTau, astau_/const_->kPi).real(),
              8.9725951245272245e-2, 1e-14);
  EXPECT_NEAR((*amu_)(q2, const_->kSTau, astau_/const_->kPi).imag(),
              -1.8025823222497021e-2, 1e-14);

  q2 = complex<double>(1.5, 2.);
  EXPECT_NEAR((*amu_)(q2, const_->kSTau, astau_/const_->kPi).real(),
              9.9728959436126022e-2, 1e-14);
  EXPECT_NEAR((*amu_)(q2, const_->kSTau, astau_/const_->kPi).imag(),
              -2.8233893471810996e-2, 1e-14);

  q2 = complex<double>(2., 0.);
  EXPECT_NEAR((*amu_)(q2, const_->kSTau, astau_/const_->kPi).real(),
              0.11721885219529733, 1e-14);

  q2 = complex<double>(2.6, 0.);
  EXPECT_NEAR((*amu_)(q2, const_->kSTau, astau_/const_->kPi).real(),
              0.10764199710166374, 1e-14);
}
