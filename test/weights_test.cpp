#include <gtest/gtest.h>
#include "../src/weights.hpp"
#include <complex>


TEST(weights_test, wRs) {
  EXPECT_DOUBLE_EQ(Weight::wR02(std::complex<double>(3., 0.)).real(), 252.);
}

TEST (weights_test, wDs) {
  EXPECT_DOUBLE_EQ(Weight::wD00(std::complex<double>(2., 0.)).real(), -3.);
  EXPECT_DOUBLE_EQ(Weight::wD02(std::complex<double>(2., 0.)).real(), -142.0/15.0);
}
