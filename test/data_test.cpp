#include <string>
#include <gtest/gtest.h>
#include "../src/data.hpp"

std::string projectRoot = "/Users/knowledge/Developer/PhD/FESR";
const Data data = readData(80, projectRoot+"/aleph.json");


// ! Caveat: Fortran array start at 1 (non at 0 like arr[0])


// compare aleph data (sbin) from Matthias
// from fesr_aleph_2015/VAmomAleph2014.f90
TEST (data_test, sbin) {
  EXPECT_EQ(data.sbin[0], 3.7499999999999999e-2); // Matthias: dsbin(1)
  EXPECT_EQ(data.sbin[16], 0.48749999999999999);  // Matthias: dsbin(17)
}

// compare aleph data (dsbin) from Matthias
// from fesr_aleph_2015/VAmomAleph2014.f90
TEST (data_test, dsbin) {
  EXPECT_EQ(data.dsbin[32], 2.5000000000000001e-2); // Matthias: dsbin(33)
  EXPECT_EQ(data.dsbin[17], 2.5000000000000001e-2); // Matthias: dsbin(18)
}

// compare aleph data (sfm2) from Matthias
// from fesr_aleph_2015/VAmomAleph2014.f90
TEST (data_test, sfm2) {
  // ! Caveat: sfm2 will be renormalized (multiplied by a factor)
  EXPECT_EQ(data.sfm2[0], 2.6332999999999999e-4); // Matthias: sfm2(1)
  EXPECT_EQ(data.sfm2[29], 0.59462000000000004); // Matthias: sfm2(30)
}

// compare aleph data (derr) from Matthias
// from fesr_aleph_2015/VAmomAleph2014.f90
TEST (data_test, derr) {
  EXPECT_EQ(data.derr[18], 3.83239999999999997e-2); // Matthias: derr(19)
  EXPECT_EQ(data.derr[53], 1.4302000000000000e-2); // Matthias: derr(54)
}

// compare aleph data (corerr) from Matthias
// from fesr_aleph_2015/VAmomAleph2014.f90
TEST (data_test, corerr) {
  EXPECT_EQ(data.corerr(0, 0), 100 ); // Matthias: correrr(1, 1)
  EXPECT_EQ(data.corerr(43, 11), 0.75721001625061035); // Matthias: corerr(44, 12)
}




