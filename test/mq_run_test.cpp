#include "../src/types.hpp"
#include "../src/configuration.hpp"
#include "../src/mq_run.hpp"
#include <gtest/gtest.h>

class MQRunTest : public ::testing::Test {
protected:
  MQRun *mq_;
  const double sTau_ = 3.1570893124;
  virtual void SetUp() {
    mq_ = new MQRun(sTau_);
  }
};

TEST_F(MQRunTest, mq) {
  complex<double> q2(3.0, 0);
  cout.precision(17);
  EXPECT_NEAR((*mq_)(q2, sTau_, 0.31927/M_PI).real(), 1.0082569029387689, 1e-6);
  q2 = complex<double>(3.0, 1.5);
  EXPECT_NEAR((*mq_)(q2, sTau_, 0.31927/M_PI).real(), 0.98122985881717040, 1e-6);
}
