#include <gtest/gtest.h>
#include "../src/constants.hpp"

// compare constants to Matthias params
// from fesr_aleph_2015/params.f90
TEST (constants_test, various_constants) {
  EXPECT_EQ(kBe, 17.827);
  EXPECT_EQ(kFPi, 92.21e-3);
  EXPECT_EQ(kSEW, 1.0198);
  EXPECT_EQ(kSTauMass, 3.1570893124000001);
  EXPECT_EQ(kPionMinusMass, 0.13957018);
  EXPECT_EQ(kVud, 0.97425);
}
