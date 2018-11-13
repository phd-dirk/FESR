#include <gtest/gtest.h>
#include "../src/numerics.hpp"
#include <boost/numeric/ublas/matrix.hpp>

using boost::numeric::ublas::matrix;


TEST(NumericsTest, matrixInvert) {
  matrix<double> m(2, 2);
  matrix<double> mR(2, 2);
  m(0, 0) = 1.;
  m(0, 1) = 2.;
  m(1, 0) = 3.;
  m(1, 1) = 4.;
  Numerics::invertMatrix(m, mR);
  EXPECT_DOUBLE_EQ(mR(0, 0), -2);
  EXPECT_DOUBLE_EQ(mR(0, 1), 1.);
  EXPECT_DOUBLE_EQ(mR(1, 0), 3./2.);
  EXPECT_DOUBLE_EQ(mR(1, 1), -1./2.);

  m(0, 0) = 0.56;
  m(0, 1) = 0.123;
  m(1, 0) = 123;
  m(1, 1) = 0.002211;
  Numerics::invertMatrix(m, mR);
  EXPECT_NEAR(mR(0, 0), -0.000146155, 1e-9);
  EXPECT_NEAR(mR(0, 1), 0.00813075, 1e-6);
}
