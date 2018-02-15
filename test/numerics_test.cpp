#include <gtest/gtest.h>
#include "../src/numerics.hpp"
#include <boost/numeric/ublas/matrix.hpp>

using boost::numeric::ublas::matrix;

class NumericsTest : public ::testing::Test {
 protected:
  Numerics * numerics;
  virtual void SetUp() {
    numerics = new Numerics(1e-14, 1e-14);
  }

  virtual void TearDown() {
    delete numerics;
  }
};

TEST_F(NumericsTest, matrixInvert) {
  matrix<double> m(2, 2);
  matrix<double> mR(2, 2);
  m(0, 0) = 1.;
  m(0, 1) = 2.;
  m(1, 0) = 3.;
  m(1, 1) = 4.;
  numerics->invertMatrix(m, mR);
  EXPECT_DOUBLE_EQ(mR(0, 0), -2);
  EXPECT_DOUBLE_EQ(mR(0, 1), 1.);
  EXPECT_DOUBLE_EQ(mR(1, 0), 3./2.);
  EXPECT_DOUBLE_EQ(mR(1, 1), -1./2.);
}
