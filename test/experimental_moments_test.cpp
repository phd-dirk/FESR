#include <gtest/gtest.h>
#include "../src/experimentalMoments.hpp"
#include "../src/s0_sets.hpp"
#include "../src/weights.hpp"

class ExperimentalMomentsTest : public ::testing::Test {
 protected:
  ExperimentalMoments * expMom;
  virtual void SetUp() {
    expMom = new ExperimentalMoments("/Users/knowledge/Developer/PhD/FESR/aleph.json",
                                     0.99363, s0Set, wD00, wD00);
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
  EXPECT_NEAR(expMoments[1], 2.8255554004717451, 1.e-14);
  EXPECT_NEAR(expMoments[7], 2.6253460058904214, 1.e-14);
}

TEST_F(ExperimentalMomentsTest, experimentalMomentPlusPionPoleMoments) {
  vector<double> expPlusPionPoleMoments = expMom->getExpPlusPionPoleMoments();
  EXPECT_NEAR(expPlusPionPoleMoments[1], 3.4673859072186186, 1.e-14);
}

TEST_F(ExperimentalMomentsTest, kPiFac) {
  EXPECT_NEAR(expMom->kDPiFac(), 6.0937786116003652e-3  ,1.e-14);
}

TEST_F(ExperimentalMomentsTest, errorMatrix) {
  EXPECT_NEAR(expMom->getErrorMatrix(0, 0), 2.1986158042263601e-7, 1.e-15);
  EXPECT_NEAR(expMom->getErrorMatrix(21, 48), -1.660774631171531e-5, 1.e-15);
  EXPECT_NEAR(expMom->getErrorMatrix(0, 79), -8.8169050111305228e-11, 1.e-15);
  EXPECT_NEAR(expMom->getErrorMatrix(80, 80), 1.6e-3, 1.e-15);
  EXPECT_NEAR(expMom->getErrorMatrix(81, 81), 3.7134137767198075e-5, 1.e-15);
}

TEST_F(ExperimentalMomentsTest, jacobianMatrix) {
  EXPECT_NEAR(expMom->getJacobianMatrix(0, 0), 5.6094687833062200e-2, 1.e-15);
  EXPECT_NEAR(expMom->getJacobianMatrix(21, 4), 6.8565471537930633e-2, 1.e-14);
  EXPECT_NEAR(expMom->getJacobianMatrix(79, 0), 0., 1.e-14);
  EXPECT_NEAR(expMom->getJacobianMatrix(80, 0), -0.15981845405157338, 1.e-14);
  EXPECT_NEAR(expMom->getJacobianMatrix(81, 2), 0.35273831585987608, 1.e-14);
}

TEST_F(ExperimentalMomentsTest, covarianceMatrix) {
  EXPECT_NEAR(expMom->getCovarianceMatrix(0, 0), 1.3667648017148091e-4, 1.e-9);
  EXPECT_NEAR(expMom->getCovarianceMatrix(1, 1), 1.1074282562391433e-4, 1.e-14);
}
