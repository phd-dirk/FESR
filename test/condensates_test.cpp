#include "../src/condensates.hpp"
#include <gtest/gtest.h>


class CondensatesTest : public ::testing::Test {
protected:
  Condensates condensates_;

  virtual void SetUp() {
    condensates_ = Condensates(
      0.3156,
      { -pow(0.272, 3), -pow(0.272, 3), 0.8*(-pow(0.272, 3)) },
      { 2.8e-3, 5.0e-3, 97.0e-3 }
    );
  }
};

TEST_F(CondensatesTest, condensates) {
  EXPECT_NEAR(condensates_.qqInv_[0], -2.0123640616267518E-002, 1e-15);
}
