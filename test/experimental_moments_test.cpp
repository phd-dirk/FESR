#include <gtest/gtest.h>
#include "../src/types.hpp"
#include "../src/experimentalMoments.hpp"
#include "../src/weights.hpp"
#include "json.hpp"
using json = nlohmann::json;
#include <fstream>

using std::vector;

class ExperimentalMomentsTest : public ::testing::Test {
 protected:
  ExperimentalMoments *expMom;
  ExperimentalMoments *expMom2;


  virtual void SetUp() {
    Weight weight(1);
    std::ifstream configFile("./test/configuration_test.json");
    json config;
    configFile >> config;

    expMom = new ExperimentalMoments("/Users/knowledge/Developer/PhD/FESR/aleph.json", config);

    std::ifstream configFile2("./test/configuration_test2.json");
    json config2;
    configFile2 >> config2;

    expMom2 = new ExperimentalMoments("/Users/knowledge/Developer/PhD/FESR/aleph.json", config2);
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

// TEST_F(ExperimentalMomentsTest, weightRatios) {
//   matrix<double> weightRatios = expMom->weightRatios;
//   EXPECT_NEAR(weightRatios(0, 0), 1., 1.e-15);
//   EXPECT_NEAR(weightRatios(1, 1), 0.99968000816184477, 1.e-14);
//   EXPECT_NEAR(weightRatios(2, 2), 0.99851573355014112, 1.e-14);
//   EXPECT_NEAR(weightRatios(1, 0), 0.99994042624850787, 1.e-14);
// }

// TEST_F(ExperimentalMomentsTest, ExperimentalMoments) {
//   vector<double> expMoments = expMom->getExperimentalMoments();
//   EXPECT_NEAR(expMoments[0], 2.8673485909014409, 1.e-15);
//   EXPECT_NEAR(expMoments[3], 2.7726476698371396, 1.e-15);
//   EXPECT_NEAR(expMoments[7], 2.6421766712865162, 1.e-14);
// }

TEST_F(ExperimentalMomentsTest, expMom) {
  vec expMoms = (*expMom)();
  EXPECT_NEAR(expMoms[0], 3.4801214322769858, 1.e-14);
  EXPECT_NEAR(expMoms[1], 3.5645062448703584, 1.e-14);
  EXPECT_NEAR(expMoms[2], 3.6634484709047848, 1.e-14);
  EXPECT_NEAR(expMoms[3], 0.55315051406537252, 1.e-14);
  EXPECT_NEAR(expMoms[4], 0.54867716084362927, 1.e-14);
}

TEST_F(ExperimentalMomentsTest, kPiFac) {
  EXPECT_NEAR(expMom->kDPiFac(), 6.0937786116003652e-3  ,1.e-14);
}

TEST_F(ExperimentalMomentsTest, errorMatrix) {
  mat errMat = expMom->errMat();
  EXPECT_NEAR(errMat(0, 0), 2.2154942818661503e-7, 1.e-15);
  EXPECT_NEAR(errMat(21, 48), -1.6735241745083317e-5, 1.e-15);
  EXPECT_NEAR(errMat(80, 80), 5.28999999999999999e-4, 1.e-15);
  EXPECT_NEAR(errMat(81, 81), 3.7134137767198075e-5, 1.e-15);
}

TEST_F(ExperimentalMomentsTest, jacobianMatrix) {
  mat jacMat = expMom->jacMat();
  EXPECT_NEAR(jacMat(0, 0), 5.9068224019943202e-2, 1.e-15);
  EXPECT_NEAR(jacMat(80, 1), -0.14605172447424614, 1.e-14);
  EXPECT_NEAR(jacMat(50, 2), 0.0, 1.e-14);
  EXPECT_NEAR(jacMat(34, 3), 1.3217232991932832e-2, 1.e-14);
  EXPECT_NEAR(jacMat(21, 4), 1.792620664346542e-2, 1.e-14);
}

TEST_F(ExperimentalMomentsTest, covarianceMatrix) {
  mat covMat = expMom->getCovMat();
  EXPECT_NEAR(covMat(0, 0), 8.4576517915615828e-5, 1.e-15);
  EXPECT_NEAR(covMat(4, 2), 1.9349658295540526e-5, 1.e-15);
  EXPECT_NEAR(covMat(1, 3), 9.5344778533946512e-6, 1.e-15);
  EXPECT_NEAR(covMat(2, 2), 2.5799227204695101e-4, 1.e-15);
}

// TEST_F(ExperimentalMomentsTest, inverseCovarianceMatrix) {
//   EXPECT_NEAR(expMom->inverseCovarianceMatrix(0, 0), 19290.123456790123, 1.e-12);
//   EXPECT_NEAR(expMom->inverseCovarianceMatrix(1, 1), 4358547.9558245243, 1.e-2);
//   EXPECT_NEAR(expMom->inverseCovarianceMatrix(6, 7), -5610057155.6616163, 1e1);
// }
