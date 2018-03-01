#include "../src/theoretical_moments.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <complex>

using std::pow;
using std::complex;

class AdlerFunctionTest : public ::testing::Test {
 protected:
  AdlerFunction *adler;
  virtual void SetUp() {
    adler = new AdlerFunction(3, 3, 5);
  }

  virtual void TearDown() {
    delete adler;
  }
};

TEST_F(AdlerFunctionTest, AdlerCoefficients) {
  // EXPECT_NEAR(adler->getC(0, 0), 1., 1e-14);
  EXPECT_NEAR(adler->getC(0, 1), 1., 1e-14);

  EXPECT_NEAR(adler->getC(1, 1), 1., 1e-14);

  EXPECT_NEAR(adler->getC(2, 1), 1.6398212048969847, 1e-14);
  EXPECT_NEAR(adler->getC(2, 2), -1.125,  1e-14);

  EXPECT_NEAR(adler->getC(3, 1), 6.3710144831009403, 1e-13);
  EXPECT_NEAR(adler->getC(3, 2), -5.6895977110182159, 1e-14);
  EXPECT_NEAR(adler->getC(3, 3), 1.6875, 1e-14);

  EXPECT_NEAR(adler->getC(4, 1), 49.075700002947983, 1e-13);
  EXPECT_NEAR(adler->getC(4, 2), -33.091406616720278, 1e-13);
  EXPECT_NEAR(adler->getC(4, 3), 15.801594849790986, 1e-14);
  EXPECT_NEAR(adler->getC(4, 4), -2.84765625, 1e-14);

  EXPECT_NEAR(adler->getC(5, 1), 283., 1e-14);
  EXPECT_NEAR(adler->getC(5, 2), -299.17718720515239, 1e-12);
  EXPECT_NEAR(adler->getC(5, 3), 129.57753256923368, 1e-12);
  EXPECT_NEAR(adler->getC(5, 4), -40.616088412029718, 1e-13);
  EXPECT_NEAR(adler->getC(5, 5), 5.1257812500000002, 1e-14);
}

TEST_F(AdlerFunctionTest, D0) {
  complex<double> s(3., 3.);
  complex<double> mu2(3.0, 0);
  EXPECT_NEAR(adler->D0(s, mu2).real(), 2.6878635987293748e-2, 1e-13);
  EXPECT_NEAR(adler->D0(s, mu2).imag(), 1.5454361164073294e-3, 1e-13);

  s = complex<double>(7.0, 2.0);
  mu2 = complex<double>(1.5, 2.2);
  EXPECT_NEAR(adler->D0(s, mu2).real(), 2.6732096654905575e-2, 1e-13);
  EXPECT_NEAR(adler->D0(s, mu2).imag(), -1.0995430321549232e-3, 1e-13);

  // complex<double> mu2_2(3., 3.);
  // EXPECT_NEAR(adler->D0(s, mu2_2).real(), 2.6878635987293748e-2, 1e-13);
}
