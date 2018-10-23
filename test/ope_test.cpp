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
  ThMoms *thMom_;
  OPE *ope_;
  Configuration *config;
  matrix<double> c;

  virtual void SetUp() {
    config = new Configuration("./test/configuration_test.json");
    ope_ = new OPE(*config);
    thMom_ = new ThMoms(*config);
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

  EXPECT_NEAR(OPE::D0(3.0, 3.0, pow(1.77686, 2), 0.3156, c, 5).real(), 2.5734381131592217E-002, 1e-15);
}

TEST_F(OPETest, D0CInt) {
  //D0CIntFO(s0, weight, sTau, aStau, c, order)
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
  EXPECT_NEAR(OPE::D0CIntFO(3.0, Weight(1), pow(1.77686, 2), 0.3156, c, 5), 1.8062025328897611, 1e-14);

  // CIPT
  EXPECT_NEAR(ope_->D0CIntCI(2.6, Weight(1), 3.1570893124000001, 0.31927, c, 1), 1.7203914169060759, 1e-13);
  EXPECT_NEAR(ope_->D0CIntCI(2.6, Weight(1), 3.1570893124000001, 0.31927, c, 2), 1.7643779920232532, 1e-13);
  EXPECT_NEAR(ope_->D0CIntCI(2.6, Weight(1), 3.1570893124000001, 0.31927, c, 5), 1.8005358682689818, 1e-13);
}

TEST_F(OPETest, D4) {
  // D4(s, mu, astau, aGGinv, order, r)
  EXPECT_NEAR(
    OPE::D4(
      3.0, 3.0, pow(1.77686, 2), 0.3156, 2.1e-2, 5, 1, { 2.8e-3, 5.0e-3, 97.0e-3 },
      Condensates(
        0.3156,
        { -pow(0.272, 3), -pow(0.272, 3), 0.8*(-pow(0.272, 3)) },
        { 2.8e-3, 5.0e-3, 97.0e-3 }
      )
    ).real(), 3.0537543015491460e-4, 1e-15
  );
  EXPECT_NEAR(
    OPE::D4(
      3.0, 3.0, pow(1.77686, 2), 0.3156, 2.1e-2, 5, 1, { 2.8e-3, 5.0e-3, 97.0e-3 },
      Condensates(
        0.3156,
        { -pow(0.272, 3), -pow(0.272, 3), 0.8*(-pow(0.272, 3)) },
        { 2.8e-3, 5.0e-3, 97.0e-3 }
      )
    ).imag(), 3.3420812507851665e-5, 1e-15
  ); // imag has minus sign?
}

TEST_F(OPETest, D4CInt) {
  // D4CInt(s0, weight, sTau, astau, aGGinv, r)
  EXPECT_NEAR(
    OPE::D4CInt(
      3.0, Weight(1), pow(1.77686, 2), 0.3156, 2.1e-2, 1, { 2.8e-3, 5.0e-3, 97.0e-3 },
      Condensates(
        0.3156,
        { -pow(0.272, 3), -pow(0.272, 3), 0.8*(-pow(0.272, 3)) },
        { 2.8e-3, 5.0e-3, 97.0e-3 }
      )
    ), 1.889902422930036e-3, 1e-15
  );
}

TEST_F(OPETest, D68) {
  // D68(s, rho, c8)
  EXPECT_NEAR(OPE::D68(3.0, -0.1894, 0.16315).real(), -1.2987654320987656e-4, 1e-15);
  EXPECT_NEAR(OPE::D68(2.0, -0.1894, 0.16315).real(), -3.02375e-4, 1e-15);
  EXPECT_NEAR(OPE::D68(2.0, -0.5, 0.16315).real(), -1.467124999999999999e-3, 1e-15);
  EXPECT_NEAR(OPE::D68(2.0, -0.5, 0.1).real(), -1.6249999999999999999999e-3, 1e-15);

  EXPECT_NEAR(OPE::D68(3.0, -0.5, 0.1).real(), -5.0617283950617285E-004, 1e-15);
  EXPECT_NEAR(OPE::D68(3.0, -0.5, 0.1).imag(), 0.0, 1e-15);
}

TEST_F(OPETest, D68CInt) {
  // D68CInt(s0, weight, rho, c8)
  EXPECT_NEAR(OPE::D68CInt(3.0, Weight(1), -0.1894, 0.16315), -2.9695080856551679e-2, 1e-15);
  EXPECT_NEAR(OPE::D68CInt(2.0, Weight(1), -0.1894, 0.16315), -0.10827202768105049, 1e-15);
  EXPECT_NEAR(OPE::D68CInt(2.0, Weight(6), -0.1894, 0.16315), 8.1905379523540371e-3, 1e-15);
  EXPECT_NEAR(OPE::D68CInt(2.0, Weight(6), -0.7, 0.16315), -6.7400762155589364e-2, 1e-15);
  EXPECT_NEAR(OPE::D68CInt(2.0, Weight(6), -0.7, 0.9), 9.6228642910621415e-2, 1e-15);

  EXPECT_NEAR(OPE::D68CInt(3.0, Weight(1), -0.5, 0.1), -6.8721689903881375E-002, 1e-15);
}

