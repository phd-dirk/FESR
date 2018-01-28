#include <gtest/gtest.h>
#include "../src/pion_pole.hpp"
#include "../src/constants.hpp"
#include "../src/weights.hpp"
#include "../src/data.hpp"

const Data data(80, "/Users/knowledge/Developer/PhD/FESR/aleph.json");

TEST (pion_pole_test, pifac) {
  EXPECT_DOUBLE_EQ(kPiFac, 1.9494983309486447);
}

TEST (pion_pole_test, pion_pole_spectral_momentum) {
  EXPECT_DOUBLE_EQ(pionPoleSpectralMoment(3., kPionMinusMass, kSTauMass, wR00,
                                          data.sbins, data.dsbins),
                  0.64183050674687347);
}
