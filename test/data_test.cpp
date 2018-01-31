#include <string>
#include <gtest/gtest.h>
#include "../src/data.hpp"
#include "../src/s0_sets.hpp"
#include "../src/constants.hpp"
#include "../src/weights.hpp"
#include <iostream>


using std::cout;
using std::endl;
const Data data("/Users/knowledge/Developer/PhD/FESR/aleph.json", 1.);

// ! Caveat: Fortran array start at 1 (non at 0 like arr[0])

// compare aleph data (sbin) from Matthias
// from fesr_aleph_2015/VAmomAleph2014.f90
TEST (data_test, sbin) {
  EXPECT_EQ(data.sbins[0], 3.7499999999999999e-2); // Matthias: dsbin(1)
  EXPECT_EQ(data.sbins[16], 0.48749999999999999);  // Matthias: dsbin(17)
}

// compare aleph data (dsbin) from Matthias
// from fesr_aleph_2015/VAmomAleph2014.f90
TEST (data_test, dsbin) {
  EXPECT_EQ(data.dsbins[32], 2.5000000000000001e-2); // Matthias: dsbin(33)
  EXPECT_EQ(data.dsbins[17], 2.5000000000000001e-2); // Matthias: dsbin(18)
}

// compare aleph data (sfm2) from Matthias
// from fesr_aleph_2015/VAmomAleph2014.f90
TEST (data_test, sfm2) {
  // ! Caveat: sfm2 will be renormalized (multiplied by a factor)
  EXPECT_EQ(data.sfm2s[0], 2.6332999999999999e-4); // Matthias: sfm2(1)
  EXPECT_EQ(data.sfm2s[29], 0.59462000000000004); // Matthias: sfm2(30)
}

// compare aleph data (derr) from Matthias
// from fesr_aleph_2015/VAmomAleph2014.f90
TEST (data_test, derr) {
  EXPECT_EQ(data.derrs[18], 3.83239999999999997e-2); // Matthias: derr(19)
  EXPECT_EQ(data.derrs[53], 1.4302000000000000e-2); // Matthias: derr(54)
}

// compare aleph data (corerr) from Matthias
// from fesr_aleph_2015/VAmomAleph2014.f90
TEST (data_test, corerr) {
  EXPECT_EQ(data.corerrs(0, 0), 100 ); // Matthias: correrr(1, 1)
  EXPECT_EQ(data.corerrs(43, 11), 0.75721001625061035); // Matthias: corerr(44, 12)
}


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

TEST_F(ExperimentalMomentsTest, errorMatrix) {
  EXPECT_NEAR(expMom->getErrorMatrix(0, 0), 2.1986158042263601e-7, 1.e-15);
  EXPECT_NEAR(expMom->getErrorMatrix(21, 48), -1.660774631171531e-5, 1.e-15);
}

TEST_F(ExperimentalMomentsTest, jacobianMatrix) {
  EXPECT_NEAR(expMom->getJacobianMatrix(0, 0), 5.6094687833062200e-2, 1.e-15);
  EXPECT_NEAR(expMom->getJacobianMatrix(21, 4), 6.8565471537930633e-2, 1.e-14);
}
