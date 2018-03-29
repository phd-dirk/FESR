#include "../src/constants.hpp"
#include "../src/mq_run.hpp"
#include <gtest/gtest.h>

class MQRunTest : public ::testing::Test {
protected:
  Constants *const_;
  MQRun *mq_;
  const double astau_ = 0.31927;
  double atau_;
  virtual void SetUp() {
    const int nc = 3;
    const int nf = 3;
    const_ = new Constants(nc, nf);
    mq_ = new MQRun(*const_);
    atau_ = astau_/const_->kPi;
  }
};

TEST_F(MQRunTest, mq) {
  complex<double> q2(0., 1.5);
  EXPECT_NEAR((*mq_)(q2, const_->kSTau, atau_).real(), 0.97537991922188128, 1e-13);
  q2 = complex<double>(2., 0);
  EXPECT_NEAR((*mq_)(q2, const_->kSTau, atau_).real(), 1.0839550043211856, 1e-13);
}
