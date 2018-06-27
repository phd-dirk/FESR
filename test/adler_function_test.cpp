#include "../src/types.hpp"
#include "../src/configuration.hpp"
#include "../src/theoretical_moments.hpp"
#include "../src/weights.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <complex>
#include "json.hpp"
#include <fstream>

using json = nlohmann::json;
using std::pow;
using std::complex;

class AdlerFunctionTest : public ::testing::Test {
protected:
  TheoreticalMoments *thMom_;
  AdlerFunction *adler;
  const double astau_ = 0.31927;
  const double aGGinv_ = 2.1e-2;
  const double sTau_ = 3.1570893124000001;
  virtual void SetUp() {
    std::ifstream configFile("./test/configuration_test.json");
    json jsonConfig;
    configFile >> jsonConfig;
    Configuration config(jsonConfig);

    adler = new AdlerFunction(config);
    thMom_ = new TheoreticalMoments(config);
  }
};

TEST_F(AdlerFunctionTest, D0) {
  cmplx s(3.0, 3.0);
  cmplx mu2(3.0, 0);
  int order = 5;
  EXPECT_NEAR(adler->D0(s, mu2, 3.1570893124000001, 0.3, 5).real(), 2.6951237510689965e-2, 1e-15);
  EXPECT_NEAR(adler->D0(s, mu2, 3.1570893124000001, 0.3, 5).imag(), 1.3861376488833205e-3, 1e-15);

  EXPECT_NEAR(adler->D0(cmplx(3.0, 3.0), cmplx(3.0, 1.0), 3.1570893124000001, 0.3, 5).real(), 2.689493917535999e-2, 1e-15);
  EXPECT_NEAR(adler->D0(cmplx(3.0, 3.0), cmplx(3.0, 1.0), 3.1570893124000001, 0.3, 5).imag(), 1.4938194348006232e-3, 1e-15);
}

TEST_F(AdlerFunctionTest, CIntD0) {
  //D0CIntFO(s0, weight, sTau, aStau, order)
  EXPECT_NEAR(adler->D0CIntFO(3.0, Weight(1), 3.1570893124000001, 0.31927, 5), 1.8130325533146323, 1e-14);
  EXPECT_NEAR(adler->D0CIntFO(2.6, Weight(1), 3.1570893124000001, 0.31927, 5), 1.8398156114670496 , 1e-14);
  EXPECT_NEAR(adler->D0CIntFO(2.6, Weight(6), 3.1570893124000001, 0.31927, 5), 0.28724759227498259, 1e-14);
  EXPECT_NEAR(adler->D0CIntFO(2.6, Weight(6), 3.1570893124000001, 0.29, 5), 0.28350497970822786, 1e-14);

  // s0 = 2.;
  // EXPECT_NEAR(adler->D0CIntFO(s0, *weight_, sTau_, astau_, order_), 1.9023472322728381, 1e-13);
}

// TEST_F(AdlerFunctionTest, D2) {
//   const complex<double> s(3.1570893123374919, -1.9866720523795795e-5);
//   const complex<double> mu2(3.1570893124, 0.);
//   const double astau = 0.31927;
//   const int order = 2;
//   const uint r = 1;
//   EXPECT_NEAR(adler->D2(s, mu2, astau, order, r).real(), 3.4324915657182825e-7, 1e-09);
// }

// TEST_F(AdlerFunctionTest, D4) {
//   complex<double> s(3.1570893123374919, 1.9866720523795795e-5);
//   complex<double> mu2(3.1570893124, 0.);
//   EXPECT_NEAR(adler->D4(s, mu2, astau_, aGGinv_, order_, 1).real(), 2.7591458364939887e-4, 1e-15);
//   s = complex<double>(2.8735226351854122, -1.3077004976474509);
//   mu2 = complex<double>(3.1570893124000001, 0.);
//   EXPECT_NEAR(adler->D4(s, mu2, astau_, aGGinv_, order_, 1).real(), 1.6118530368558565e-4, 1e-15);
//   EXPECT_NEAR(adler->D4(s, mu2, astau_, aGGinv_, order_, 1).imag(), 2.2444197547483058e-4, 1e-15);
//   s = complex<double>(-2.1537305702568648, 2.3083885195545712);
//   mu2 = complex<double>(3.1570893124000001, 0.);
//   EXPECT_NEAR(adler->D4(s, mu2, astau_, aGGinv_, order_, 1).real(), -1.1095466575342775e-5, 1e-15);
//   EXPECT_NEAR(adler->D4(s, mu2, astau_, aGGinv_, order_, 1).imag(), 2.7283500628353362e-4, 1e-15);
//   EXPECT_NEAR(adler->D4(s, mu2, astau_, aGGinv_, order_, -1).real(), -1.6819854236697827e-5, 1e-15);
//   EXPECT_NEAR(adler->D4(s, mu2, astau_, aGGinv_, order_, -1).imag(), 2.8808404431198147e-4, 1e-15);

//   s = complex<double>(-1.3643773470694827, 1.4623523702595214);
//   mu2 = complex<double>(2., 0.);
//   EXPECT_NEAR(adler->D4(s, mu2, astau_, aGGinv_, order_, 1).real(), -1.8100662040598564e-5, 1e-15);

// }
// TEST_F(AdlerFunctionTest, D4CInt) {
//   double s0 = 3.1570893124;
//   EXPECT_NEAR(adler->D4CInt(s0, *weight_, sTau_, astau_, aGGinv_, 1), 1.6944347548019322e-3, 1e-15);
//   // EXPECT_NEAR(adler->D4CInt(s0, *weight_, astau_, aGGinv_, -1), 6.8601706024321041e-4, 1e-15);
//   // s0 = 2.6;
//   // EXPECT_NEAR(adler->D4CInt(s0, *weight_, astau_, aGGinv_, 1), 2.8575812002800097e-3, 1e-15);
//   // s0 = 2.0;
//   // EXPECT_NEAR(adler->D4CInt(s0, *weight_, astau_, aGGinv_, 1), 5.9016068520545409e-3, 1e-15);
//   // EXPECT_NEAR(adler->D4CInt(s0, *weight_, astau_, aGGinv_, -1), 2.2220948788460476e-3, 1e-15);
// }

// TEST_F(AdlerFunctionTest, D68) {
//   const complex<double> s(-1.7585656673884618, 1.9150579086499120);
//   const double rho = -0.1894;
//   const double c8 = 0.16315;
//   EXPECT_NEAR(adler->D68(s, rho, c8).real(), -3.9659215963967076e-4, 1e-13);
//   EXPECT_NEAR(adler->D68(s, rho, c8).imag(), 1.7341429761413442e-4, 1e-13);
// }

// TEST_F(AdlerFunctionTest, D68CInt) {
//   const double s0 = 2.4;
//   const double rho = -0.1894;
//   const double c8 = 0.16315;
//   EXPECT_NEAR(adler->D68CInt(s0, *weight_, rho, c8), -6.0327814112243174e-2, 1e-13);
// }
