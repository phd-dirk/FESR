#include <gtest/gtest.h>
#include "../src/constants.hpp"
#include "../src/alpha_s.hpp"
#include <complex>
#include <cmath>
#include "json.hpp"
#include <fstream>
using json = nlohmann::json;

using std::complex;
using std::sqrt;


class AlphaSTest : public ::testing::Test {
 protected:
  Constants *const_;
  AlphaS *amu_;
  const double astau_ = 0.31927;

  virtual void SetUp() {
    std::ifstream configFile("./test/configuration_test.json");
    json config;
    configFile >> config;
    const_ = new Constants(config);
    amu_ = new AlphaS(*const_);
  }
};

TEST_F(AlphaSTest, complex) {
  // BETA 4th order --------------------------------------------------
  // complex<double> q2(3., 3.);
  // EXPECT_NEAR((*amu_)(q2, const_->kSTau, astau_/const_->kPi).real(),
  //             8.9725951245272245e-2, 1e-14);
  // EXPECT_NEAR((*amu_)(q2, const_->kSTau, astau_/const_->kPi).imag(),
  //             -1.8025823222497021e-2, 1e-14);

  // q2 = complex<double>(1.5, 2.);
  // EXPECT_NEAR((*amu_)(q2, const_->kSTau, astau_/const_->kPi).real(),
  //             9.9728959436126022e-2, 1e-14);
  // EXPECT_NEAR((*amu_)(q2, const_->kSTau, astau_/const_->kPi).imag(),
  //             -2.8233893471810996e-2, 1e-14);

  // q2 = complex<double>(2., 0.);
  // EXPECT_NEAR((*amu_)(q2, const_->kSTau, astau_/const_->kPi).real(),
  //             0.11721885219529733, 1e-14);

  // q2 = complex<double>(2.6, 0.);
  // EXPECT_NEAR((*amu_)(q2, const_->kSTau, astau_/const_->kPi).real(),
  //             0.10764199710166374, 1e-14);

  // BETA 5th order --------------------------------------------------
  complex<double> q(3.0, 0.0);
  EXPECT_NEAR((*amu_)(q, 3.1572314596, 0.10162679736189885).real(),
              0.10470997092165853, 1e-14);
}
