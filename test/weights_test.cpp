#include <gtest/gtest.h>
#include "../src/weights.hpp"
#include <complex>

using std::complex;

TEST (weights_test, wDs) {
  EXPECT_DOUBLE_EQ(Weight::wD00(complex<double>(2., 0.)).real(), -3.);
}
