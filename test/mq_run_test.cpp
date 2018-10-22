#include "../src/types.hpp"
#include "../src/configuration.hpp"
#include "../src/mq_run.hpp"
#include <gtest/gtest.h>

TEST(MQRunTest, mq) {
  const double sTau = 3.1570893124;
  complex<double> q2(3.0, 0);
  cout.precision(17);
  EXPECT_NEAR(MQ::run(q2, sTau, sTau, 0.31927/M_PI).real(), 1.0082569029387689, 1e-6);
  q2 = complex<double>(3.0, 1.5);
  EXPECT_NEAR(MQ::run(q2, sTau, sTau, 0.31927/M_PI).real(), 0.98122985881717040, 1e-6);
}
