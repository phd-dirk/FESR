#include <gtest/gtest.h>
#include "../src/configuration.hpp"

class ConfigurationTest : public ::testing::Test {
protected:
  Configuration *config;
  virtual void SetUp() {
    config = new Configuration("./test/configuration_test.json");
  }

  virtual void TearDown() {
    delete config;
  }
};

// compare constants to Matthias params
// from fesr_aleph_2015/params.f90
TEST_F (ConfigurationTest, various_constants) {
  EXPECT_EQ(config->be, 17.815);
  EXPECT_EQ(config->kFPi, 92.21e-3);
  EXPECT_EQ(config->kSEW, 1.0198);
  EXPECT_EQ(config->mTau, 1.77682);
  EXPECT_EQ(config->kPionMinusMass, 0.13957018);
  EXPECT_EQ(config->kVud, 0.97425);
}

TEST_F(ConfigurationTest, AdlerCoefficients) {
  // EXPECT_NEAR(adler->getC(0, 0), 1., 1e-14);
  EXPECT_NEAR(config->c[0][1], 1., 1e-14);

  EXPECT_NEAR(config->c[1][1], 1., 1e-14);

  EXPECT_NEAR(config->c[2][1], 1.6398212048969847, 1e-14);
  EXPECT_NEAR(config->c[2][2], -1.125,  1e-14);

  EXPECT_NEAR(config->c[3][1], 6.3710144831009403, 1e-13);
  EXPECT_NEAR(config->c[3][2], -5.6895977110182159, 1e-14);
  EXPECT_NEAR(config->c[3][3], 1.6875, 1e-14);

  EXPECT_NEAR(config->c[4][1], 49.075700002947983, 1e-13);
  EXPECT_NEAR(config->c[4][2], -33.091406616720278, 1e-13);
  EXPECT_NEAR(config->c[4][3], 15.801594849790986, 1e-14);
  EXPECT_NEAR(config->c[4][4], -2.84765625, 1e-14);

  EXPECT_NEAR(config->c[5][1], 283., 1e-14);
  EXPECT_NEAR(config->c[5][2], -299.17718720515239, 1e-12);
  EXPECT_NEAR(config->c[5][3], 129.57753256923368, 1e-12);
  EXPECT_NEAR(config->c[5][4], -40.616088412029718, 1e-13);
  EXPECT_NEAR(config->c[5][5], 5.1257812500000002, 1e-14);
}
