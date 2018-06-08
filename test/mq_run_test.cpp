#include "../src/constants.hpp"
#include "../src/mq_run.hpp"
#include <gtest/gtest.h>
#include "json.hpp"
#include <fstream>
using json = nlohmann::json;

class MQRunTest : public ::testing::Test {
protected:
  Constants *const_;
  MQRun *mq_;
  const double astau_ = 0.31927;
  double atau_;
  virtual void SetUp() {
    std::ifstream configFile("./test/configuration_test.json");
    json config;
    configFile >> config;
    const_ = new Constants(config);
    int order = config["parameters"]["order"];
    mq_ = new MQRun(*const_);
    atau_ = astau_/const_->kPi;
  }
};

TEST_F(MQRunTest, mq) {
  complex<double> q2(0., 1.5);
  EXPECT_NEAR((*mq_)(q2, const_->kSTau, atau_).real(), 0.97537991922188128, 1e-13);
  q2 = complex<double>(3., 0);
  EXPECT_NEAR((*mq_)(q2, const_->kSTau, atau_).real(), 1.0082152936637017, 1e-13);
  q2 = complex<double>(2., 0);
  EXPECT_NEAR((*mq_)(q2, const_->kSTau, atau_).real(), 1.0839550043211856, 1e-13);
}
