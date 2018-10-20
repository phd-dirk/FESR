#include "../src/types.hpp"
#include "../src/configuration.hpp"
#include "../src/theoretical_moments.hpp"
#include "../src/weights.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <complex>
#include <boost/numeric/ublas/matrix.hpp>

using std::pow;
using std::complex;
using boost::numeric::ublas::matrix;

class OPETest : public ::testing::Test {
protected:
  TheoreticalMoments *thMom_;
  OPE *ope_;
  Configuration *config;
  matrix<double> c;

  virtual void SetUp() {
    config = new Configuration("./test/configuration_test.json");
    ope_ = new OPE(*config);
    thMom_ = new TheoreticalMoments(*config);
    c = Configuration::adlerCoefficients(3, Configuration::betaCoefficients(3, 3));
  }
};

TEST_F(OPETest, D0) {
  cmplx s(3.0, 3.0);
  cmplx mu2(3.0, 0);
  int order = 5;
  EXPECT_NEAR(OPE::D0(s, mu2, 3.1570893124000001, 0.3, c, 5).real(), 2.6951237510689965e-2, 1e-15);
  EXPECT_NEAR(OPE::D0(s, mu2, 3.1570893124000001, 0.3, c, 5).imag(), 1.3861376488833205e-3, 1e-15);

  EXPECT_NEAR(OPE::D0(cmplx(3.0, 3.0), cmplx(3.0, 1.0), 3.1570893124000001, 0.3, c, 5).real(), 2.689493917535999e-2, 1e-15);
  EXPECT_NEAR(OPE::D0(cmplx(3.0, 3.0), cmplx(3.0, 1.0), 3.1570893124000001, 0.3, c, 5).imag(), 1.4938194348006232e-3, 1e-15);

  EXPECT_NEAR(OPE::D0(3.0, -3.0, 3.1570893124000001, 0.31927, c, 3).real(), 2.6828728340843492e-2, 1e-15);
  EXPECT_NEAR(OPE::D0(3.0, -3.0, 3.1570893124000001, 0.31927, c, 3).imag(), -1.4849720954483447e-3, 1e-15);
}

TEST_F(OPETest, D0CInt) {
  //D0CIntFO(s0, weight, sTau, aStau, order)
  EXPECT_NEAR(
    OPE::D0CIntFO(
      3.0, Weight(1), 3.1570893124000001, 0.31927, c, 5
    ), 1.8130325533146323, 1e-14
  );
  EXPECT_NEAR(
    OPE::D0CIntFO(
      2.6, Weight(1), 3.1570893124000001, 0.31927, c, 5
    ), 1.8398156114670496 , 1e-14
  );
  EXPECT_NEAR(
    OPE::D0CIntFO(
      2.6, Weight(6), 3.1570893124000001, 0.31927, c, 5
    ), 0.28724759227498259, 1e-14
  );
  EXPECT_NEAR(
    OPE::D0CIntFO(
      2.6, Weight(6), 3.1570893124000001, 0.29, c, 5
    ), 0.28350497970822786, 1e-14
  );

  // CIPT
  EXPECT_NEAR(ope_->D0CIntCI(2.6, Weight(1), 3.1570893124000001, 0.31927, 1), 1.7203914169060759, 1e-13);
  EXPECT_NEAR(ope_->D0CIntCI(2.6, Weight(1), 3.1570893124000001, 0.31927, 2), 1.7643779920232532, 1e-13);
  EXPECT_NEAR(ope_->D0CIntCI(2.6, Weight(1), 3.1570893124000001, 0.31927, 5), 1.8005358682689818, 1e-13);

  // OPE ope = OPE(
  //   3,
  //   3,
  //   { 2.8e-3, 5.0e-3, 97.0e-3 },
  //   Condensates(
  //     0.3156,
  //     { 0.0201236, 0.0201236, 0.0160989 },
  //     { 2.8e-3, 5.0e-3, 97.0e-3 }
  //   ),
  //   3.1572314596,
  //   0.13957018,
  //   92.21e-3,
  //   2.2e-3,
  //   1.3,
  //   0.4,
  //   0.19e-3,
  //   1.8,
  //   0.21
  // );

  // EXPECT_NEAR(ope.D0CIntCI(3.0, Weight(1), 3.1570893124000001, 0.31927, 5), 1.8005358682689818, 1e-13);
}

// TEST_F(OPETest, D2) {
//   const complex<double> s(3.1570893123374919, -1.9866720523795795e-5);
//   const complex<double> mu2(3.1570893124, 0.);
//   const double astau = 0.31927;
//   const int order = 2;
//   const uint r = 1;
//   EXPECT_NEAR(ope_->D2(s, mu2, astau, order, r).real(), 3.4324915657182825e-7, 1e-09);
// }

// TEST_F(OPETest, D4) {
//   // D4(s, mu, astau, aGGinv, order, r)
//   EXPECT_NEAR(ope_->D4(3.0, 3.0, 3.1570893124000001, 0.31927, 2.1e-2, 5, 1).real(), 3.0462940737588433e-4, 1e-15);
// }

// TEST_F(OPETest, D4CInt) {
//   // D4CInt(s0, weight, sTau, astau, aGGinv, r)
//   EXPECT_NEAR(ope_->D4CInt(3.0, Weight(1), 3.1570893124000001, 0.31927, 2.1e-2, 1),
//               1.942350658196885e-3, config->tolerance);
//   EXPECT_NEAR(ope_->D4CInt(3.0, Weight(6), 3.1570893124000001, 0.31927, 2.1e-2, 1),
//               -1.8595234733280692e-2 , config->tolerance);
//   EXPECT_NEAR(ope_->D4CInt(3.0, Weight(6), 3.1570893124000001, 0.28, 2.1e-2, 1),
//               -1.8945624330502554e-2, config->tolerance);
// }

TEST_F(OPETest, D68) {
  // D68(s, rho, c8)
   EXPECT_NEAR(ope_->D68(3.0, -0.1894, 0.16315).real(), -1.2987654320987656e-4, 1e-15);
   EXPECT_NEAR(ope_->D68(2.0, -0.1894, 0.16315).real(), -3.02375e-4, 1e-15);
   EXPECT_NEAR(ope_->D68(2.0, -0.5, 0.16315).real(), -1.467124999999999999e-3, 1e-15);
   EXPECT_NEAR(ope_->D68(2.0, -0.5, 0.1).real(), -1.6249999999999999999999e-3, 1e-15);
}

TEST_F(OPETest, D68CInt) {
  // D68CInt(s0, weight, rho, c8)
  EXPECT_NEAR(ope_->D68CInt(3.0, Weight(1), -0.1894, 0.16315), -2.9695080856551679e-2, config->tolerance);
  EXPECT_NEAR(ope_->D68CInt(2.0, Weight(1), -0.1894, 0.16315), -0.10827202768105049, config->tolerance);
  EXPECT_NEAR(ope_->D68CInt(2.0, Weight(6), -0.1894, 0.16315), 8.1905379523540371e-3, config->tolerance);
  EXPECT_NEAR(ope_->D68CInt(2.0, Weight(6), -0.7, 0.16315), -6.7400762155589364e-2, config->tolerance);
  EXPECT_NEAR(ope_->D68CInt(2.0, Weight(6), -0.7, 0.9), 9.6228642910621415e-2, config->tolerance);
}

TEST_F(OPETest, deltaP) {
  // deltaP(s0, weight)
  EXPECT_NEAR(ope_->deltaP(3.0, Weight(1)), -2.7752063630383842e-3, 1e-15);
  EXPECT_NEAR(ope_->deltaP(2.0, Weight(1)), -4.138431838090359e-3, 1e-15);
  EXPECT_NEAR(ope_->deltaP(2.0, Weight(6)), -4.809780763096577e-5, 1e-15);
}
