#include <gtest/gtest.h>
#include "../src/experimentalMoments.hpp"
#include "../src/weights.hpp"
#include "../src/constants.hpp"
#include "json.hpp"
using json = nlohmann::json;
#include <fstream>

using std::vector;

class ExperimentalMomentsTest : public ::testing::Test {
 protected:
  ExperimentalMoments * expMom;
  virtual void SetUp() {

    Weight weight(1);
    std::ifstream configFile("./test/configuration_test.json");
    json config;
    configFile >> config;

    Constants constants(config);
    expMom = new ExperimentalMoments("/Users/knowledge/Developer/PhD/FESR/aleph.json",
                                     1, config["parameters"]["s0Set"], weight, constants);
  }

  virtual void TearDown() {
    delete expMom;
  }
};


// compare closest bin to s0 with Matthias
// from fesr_aleph_2015/VAmomAleph2014.f90
TEST_F(ExperimentalMomentsTest, closestBinToS0) {
  // ! Caveat: The bin number corresponds to the array index ( Fortran arrays start at 1, CPP at 0)
  // e.g. Matthias maxBin: 67 => CPP maxbin: 66
  EXPECT_EQ(expMom->closestBinToS0(pow(1.77682, 2)), 78); // Matthias 79
  EXPECT_EQ(expMom->closestBinToS0(3.), 77); // Matthias 78
  EXPECT_EQ(expMom->closestBinToS0(2.1), 71); // Matthias 72
}

TEST_F(ExperimentalMomentsTest, weightRatios) {
  matrix<double> weightRatios = expMom->getWeightRatios();
  EXPECT_NEAR(weightRatios(0, 0), 1., 1.e-15);
  EXPECT_NEAR(weightRatios(1, 1), 0.99968000816184477, 1.e-14);
  EXPECT_NEAR(weightRatios(2, 2), 0.99851573355014112, 1.e-14);
  EXPECT_NEAR(weightRatios(1, 0), 0.99994042624850787, 1.e-14);
}

TEST_F(ExperimentalMomentsTest, ExperimentalMoments) {
  vector<double> expMoments = expMom->getExperimentalMoments();
  EXPECT_NEAR(expMoments[0], 2.8673485909014409, 1.e-15);
  EXPECT_NEAR(expMoments[3], 2.7726476698371396, 1.e-15);
  EXPECT_NEAR(expMoments[7], 2.6421766712865162, 1.e-14);
}

TEST_F(ExperimentalMomentsTest, experimentalMomentPlusPionPoleMoments) {
  vector<double> expPlusPionPoleMoments = expMom->getExpPlusPionPoleMoments();
  EXPECT_NEAR(expPlusPionPoleMoments[1], 3.4855000824156286, 1.e-14);
  EXPECT_NEAR(expPlusPionPoleMoments[4], 3.5298822719575593, 1.e-14);
  EXPECT_NEAR(expPlusPionPoleMoments[8], 3.5694369485452664, 1.e-14);
}

TEST_F(ExperimentalMomentsTest, kPiFac) {
  EXPECT_NEAR(expMom->kDPiFac(), 6.0937786116003652e-3  ,1.e-14);
}

TEST_F(ExperimentalMomentsTest, errorMatrix) {
  EXPECT_NEAR(expMom->getErrorMatrix(0, 0), 2.2268960999999995e-7, 1.e-15);
  EXPECT_NEAR(expMom->getErrorMatrix(21, 48), -1.6821367980824591e-5, 1.e-15);
  EXPECT_NEAR(expMom->getErrorMatrix(0, 79), -8.9303148579275616e-11, 1.e-15);
  EXPECT_NEAR(expMom->getErrorMatrix(80, 80), 1.6e-3, 1.e-15);
  EXPECT_NEAR(expMom->getErrorMatrix(81, 81), 3.7134137767198075e-5, 1.e-15);
}

TEST_F(ExperimentalMomentsTest, jacobianMatrix) {
  EXPECT_NEAR(expMom->getJacobianMatrix(0, 0), 5.6094687833062200e-2, 1.e-15);
  EXPECT_NEAR(expMom->getJacobianMatrix(21, 4), 6.8565471537930633e-2, 1.e-14);
  EXPECT_NEAR(expMom->getJacobianMatrix(79, 0), 0., 1.e-14);
  EXPECT_NEAR(expMom->getJacobianMatrix(80, 0), -0.16084302411518711, 1.e-14);
  EXPECT_NEAR(expMom->getJacobianMatrix(81, 2), 0.35273831585987608, 1.e-14);
}

TEST_F(ExperimentalMomentsTest, covarianceMatrix) {
  EXPECT_NEAR(expMom->covarianceMatrix(0, 0), 1.38387768233595e-4, 1.e-12);
  EXPECT_NEAR(expMom->covarianceMatrix(1, 1), 1.1211551248845397e-4 , 1.e-14);
  EXPECT_NEAR(expMom->covarianceMatrix(2, 6), 9.7822218694526206e-5, 1.e-14);
  EXPECT_NEAR(expMom->covarianceMatrix(7, 8), 1.2914710292441723e-4, 1.e-14);
}

TEST_F(ExperimentalMomentsTest, inverseCovarianceMatrix) {
  EXPECT_NEAR(expMom->inverseCovarianceMatrix(0, 0), 19290.123456790123, 1.e-12);
  EXPECT_NEAR(expMom->inverseCovarianceMatrix(1, 1), 4358547.9558245243, 1.e-2);
  EXPECT_NEAR(expMom->inverseCovarianceMatrix(6, 7), -5610057155.6616163, 1e1);
}
