#include <gtest/gtest.h>
#include "../src/constants.hpp"
#include "json.hpp"
#include <fstream>
using json = nlohmann::json;

class ConstantsTest : public ::testing::Test {
protected:
  Constants * constants;
  virtual void SetUp() {
    std::ifstream configFile("./test/configuration_test.json");
    json config;
    configFile >> config;

    constants = new Constants(config);
  }

  virtual void TearDown() {
    delete constants;
  }
};

// compare constants to Matthias params
// from fesr_aleph_2015/params.f90
TEST_F (ConstantsTest, various_constants) {
  EXPECT_EQ(constants->kBe, 17.827);
  EXPECT_EQ(constants->kFPi, 92.21e-3);
  EXPECT_EQ(constants->kSEW, 1.0198);
  EXPECT_EQ(constants->kMTau, 1.77682);
  EXPECT_EQ(constants->kPionMinusMass, 0.13957018);
  EXPECT_EQ(constants->kVud, 0.97425);
}

TEST_F(ConstantsTest, AdlerCoefficients) {
  // EXPECT_NEAR(adler->getC(0, 0), 1., 1e-14);
  EXPECT_NEAR(constants->c_[0][1], 1., 1e-14);

  EXPECT_NEAR(constants->c_[1][1], 1., 1e-14);

  EXPECT_NEAR(constants->c_[2][1], 1.6398212048969847, 1e-14);
  EXPECT_NEAR(constants->c_[2][2], -1.125,  1e-14);

  EXPECT_NEAR(constants->c_[3][1], 6.3710144831009403, 1e-13);
  EXPECT_NEAR(constants->c_[3][2], -5.6895977110182159, 1e-14);
  EXPECT_NEAR(constants->c_[3][3], 1.6875, 1e-14);

  EXPECT_NEAR(constants->c_[4][1], 49.075700002947983, 1e-13);
  EXPECT_NEAR(constants->c_[4][2], -33.091406616720278, 1e-13);
  EXPECT_NEAR(constants->c_[4][3], 15.801594849790986, 1e-14);
  EXPECT_NEAR(constants->c_[4][4], -2.84765625, 1e-14);

  EXPECT_NEAR(constants->c_[5][1], 283., 1e-14);
  EXPECT_NEAR(constants->c_[5][2], -299.17718720515239, 1e-12);
  EXPECT_NEAR(constants->c_[5][3], 129.57753256923368, 1e-12);
  EXPECT_NEAR(constants->c_[5][4], -40.616088412029718, 1e-13);
  EXPECT_NEAR(constants->c_[5][5], 5.1257812500000002, 1e-14);
}
