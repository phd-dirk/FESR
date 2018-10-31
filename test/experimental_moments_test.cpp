#include <gtest/gtest.h>
#include "../src/types.hpp"
#include "../src/experimentalMoments.hpp"
#include "../src/weights.hpp"

using std::vector;

class ExperimentalMomentsTest : public ::testing::Test {
 protected:
  ExpMoms *expMom;
  ExpMoms *expMom2;

  virtual void SetUp() {
    std::vector<Input> inputs = {
      { Weight(1), { 3.0, 2.0, 1.0 } },
      { Weight(6), { 2.8, 1.8 } },
    };
    ThMomContribs thMomContribs = { "FO", true, true, true, false, true };

    expMom = new ExpMoms(
      "/Users/knowledge/Developer/PhD/FESR/aleph.json",
      Configuration(
        0.3156,
        17.815,
        0.023,
        0.97425,
        0.00022,
        1.0198,
        0.0006,
        0.13957018,
        92.21e-3,
        0.14e-3,
        2.2e-3,
        1.3,
        0.4,
        0.19e-3,
        1.8,
        0.21,
        { 2.8e-3, 5.0e-3, 97.0e-3 },
        { 0.0201236, 0.021236, 0.0160989 },
        inputs,
        1.77682,
        3,
        3,
        5,
        0.99743669,
        thMomContribs
      )
    );

    std::vector<Input> inputs2 = {{
        Weight(1),
        {
          3.1572314596000002, 3.0, 2.800, 2.600, 2.400,
          2.300, 2.200, 2.100, 2.000
        }
      }
    };
    expMom2 = new ExpMoms(
      "/Users/knowledge/Developer/PhD/FESR/aleph.json",
      inputs2,
      3.1572314596000002, // sTau
      17.815, // be
      0.023, // dBe
      0.97420, // Vud
      0.00021, // dVud
      1.0198, // SEW
      0.0006, // dSEW
      92.21e-3, // fPi
      0.14e-3, // dFPi
      0.13957018, // pionMinusMass
      0.99743669 // RVANormalization
    );
  }

  virtual void TearDown() {
    delete expMom;
    delete expMom2;
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
  EXPECT_NEAR(
    expMom2->wRatio(3.1572314596000002, Weight(1), 0), 1., 1.e-15
  );
  EXPECT_NEAR(
    expMom2->wRatio(3.0, Weight(1), 53), 0.93963444985821865, 1.e-14
  );
  EXPECT_NEAR(
    expMom2->wRatio(2.0, Weight(1), 20), 0.87108121877135736, 1.e-14
  );
}

TEST_F(ExperimentalMomentsTest, expMom) {
  vec expMoms = (*expMom)();
  EXPECT_NEAR(expMoms[0], 3.4801214322769858, 1.e-14);
  EXPECT_NEAR(expMoms[1], 3.5645062448703584, 1.e-14);
  EXPECT_NEAR(expMoms[2], 3.6634484709047848, 1.e-14);
  EXPECT_NEAR(expMoms[3], 0.55315051406537252, 1.e-14);
  EXPECT_NEAR(expMoms[4], 0.54867716084362927, 1.e-14);

  // weight(1), s0 = sTau
  ExpMoms expMom1(
    "/Users/knowledge/Developer/PhD/FESR/aleph.json",
    {
      {
        Weight(1),
        {
          3.1572314596000002, 3.0, 2.800, 2.600, 2.400,
          2.300, 2.200, 2.100, 2.000
        }
      }
    },
    3.1572314596000002, // sTau
    17.815, // be
    0.023, // dBe
    0.97420, // Vud
    0.00021, // dVud
    1.0198, // SEW
    0.0006, // dSEW
    92.21e-3, // fPi
    0.14e-3, // dFPi
    0.13957018, // pionMinusMass
    0.99743669 // RVANormalization
  );

  EXPECT_NEAR(expMom1()[0], 3.4717374120768105, 1.e-14);

  // weight(3), s0 = sTau
  ExpMoms expMom2(
    "/Users/knowledge/Developer/PhD/FESR/aleph.json",
    {
      {
        Weight(3),
        {
          3.1572314596000002, 3.0, 2.800, 2.600, 2.400,
          2.300, 2.200, 2.100, 2.000
        }
      }
    },
    3.1572314596000002, // sTau
    17.815, // be
    0.023, // dBe
    0.97420, // Vud
    0.00021, // dVud
    1.0198, // SEW
    0.0006, // dSEW
    92.21e-3, // fPi
    0.14e-3, // dFPi
    0.13957018, // pionMinusMass
    0.99743669 // RVANormalization
  );

  EXPECT_NEAR(expMom2()[0], 0.44080705131582015, 1.e-14);
}

TEST_F(ExperimentalMomentsTest, piFac) {
  EXPECT_NEAR(expMom->dPiFac(), 6.0937786116003652e-3, 1.e-14);

  EXPECT_NEAR(expMom2->piFac(),  1.9492982336116860, 1.e-14);
  EXPECT_NEAR(expMom2->dPiFac(), 6.0875061675025877E-003, 1.e-14);
}

TEST_F(ExperimentalMomentsTest, errorMatrix) {
  mat errMat = expMom->errMat();
  EXPECT_NEAR(errMat(0, 0), 2.2154942818661503e-7, 1.e-15);
  EXPECT_NEAR(errMat(21, 48), -1.6735241745083317e-5, 1.e-15);
  EXPECT_NEAR(errMat(80, 80), 5.28999999999999999e-4, 1.e-15);
  EXPECT_NEAR(errMat(81, 81), 3.7134137767198075e-5, 1.e-15);

  mat errMat2 = expMom2->errMat();
  EXPECT_NEAR(errMat(0, 0), 2.2154942818661503E-007, 1.e-15);
  double errMat2Sum = 0;
  for(int i = 0; i < errMat2.size1(); i++) {
    for(int j = 0; j < errMat2.size2(); j++) {
      errMat2Sum += errMat2(i, j);
    }
  }
  EXPECT_NEAR(errMat2Sum, 3.0084246919183249E-002, 1.e-13);
}

TEST_F(ExperimentalMomentsTest, jacobianMatrix) {
  mat jacMat = expMom->jacMat();
  EXPECT_NEAR(jacMat(0, 0), 5.9068224019943202e-2, 1.e-15);
  EXPECT_NEAR(jacMat(80, 1), -0.14605172447424614, 1.e-14);
  EXPECT_NEAR(jacMat(50, 2), 0.0, 1.e-14);
  EXPECT_NEAR(jacMat(34, 3), 1.3217232991932832e-2, 1.e-14);
  EXPECT_NEAR(jacMat(21, 4), 1.792620664346542e-2, 1.e-14);

  mat jacMat2 = expMom2->jacMat();
  EXPECT_NEAR(jacMat2(0, 0), 5.6132472635419588E-002, 1.e-14);
  EXPECT_NEAR(jacMat2(27, 4), 6.5891460929989584E-002, 1.e-14);
  EXPECT_NEAR(jacMat2(80, 7), -0.14803522840184405, 1.e-14);
  EXPECT_NEAR(jacMat2(81, 3), 0.37986349105621087, 1.e-14);
  double jacMat2Sum = 0;
  for(int i = 0; i < jacMat2.size1(); i++) {
    for(int j = 0; j < jacMat2.size2(); j++) {
      jacMat2Sum += jacMat2(i, j);
    }
  }
  EXPECT_NEAR(jacMat2Sum, 39.389541422983221, 1.e-13);
}

TEST_F(ExperimentalMomentsTest, covarianceMatrix) {
  mat covMat = expMom->covMat_;
  EXPECT_NEAR(covMat(0, 0), 8.4576517915615828e-5, 1.e-15);
  EXPECT_NEAR(covMat(4, 2), 1.9349658295540526e-5, 1.e-15);
  EXPECT_NEAR(covMat(1, 3), 9.5344778533946512e-6, 1.e-15);
  EXPECT_NEAR(covMat(2, 2), 2.5799227204695101e-4, 1.e-15);

  mat covMat2 = expMom2->covMat_;
  EXPECT_NEAR(covMat2(0, 0), 1.1028649814955913e-4, 1.e-15);
  EXPECT_NEAR(covMat2(2, 6), 7.2254964242188837E-005, 1.e-15);
  double covMat2Sum = 0;
  for(int i = 0; i < covMat2.size1(); i++) {
    for(int j = 0; j < covMat2.size2(); j++) {
      covMat2Sum += covMat2(i, j);
    }
  }
  EXPECT_NEAR(covMat2Sum, 6.3547447523870154E-003, 1.e-14);
}

// TEST_F(ExperimentalMomentsTest, inverseCovarianceMatrix) {
//   mat invCov = expMom2->invCovMat_;

//   EXPECT_NEAR(invCov(0, 0), 19290.123456790123 , 1e-15);
//   EXPECT_NEAR(invCov(1, 1), 4380565.4450900145 , 1e-15);
//   // EXPECT_NEAR(invCov(4, 4), 2392765.2394170612, 1e-6);
// }
