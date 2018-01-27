#include <gtest/gtest.h>
#include "../src/constants.hpp"

// compare constants to Matthias params
// from fesr_aleph_2015/params.f90
TEST (constants_test, various_constants) {
  EXPECT_EQ(sTau, 3.1570893124000001);
  EXPECT_EQ(be, 17.827);
}
