#include <gtest/gtest.h>
#include "../src/constants.hpp"

// compare constants to Matthias params
// from fesr_aleph_2015/params.f90
TEST (constants_test, various_constants) {
  EXPECT_EQ(be, 17.827);
  EXPECT_EQ(fpi, 92.21e-3);
  EXPECT_EQ(sew, 1.0198);
  EXPECT_EQ(sTau, 3.1570893124000001);
  EXPECT_EQ(vud, 0.97425);
}
