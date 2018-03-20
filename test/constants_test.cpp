#include <gtest/gtest.h>
#include "../src/constants.hpp"

class ConstantsTest : public ::testing::Test {
protected:
  Constants * constants;
  virtual void SetUp() {
    constants = new Constants();
  }

  virtual void TearDown() {
    delete constants;
  }
};

// compare constants to Matthias params
// from fesr_aleph_2015/params.f90
TEST_F (ConstantsTest, various_constants) {
  EXPECT_EQ(Constants::kBe, 17.827);
  EXPECT_EQ(Constants::kFPi, 92.21e-3);
  EXPECT_EQ(constants->kSEW, 1.0198);
  EXPECT_EQ(Constants::kSTau, 3.1570893124000001);
  EXPECT_EQ(Constants::kPionMinusMass, 0.13957018);
  EXPECT_EQ(constants->kVud, 0.97425);
}
