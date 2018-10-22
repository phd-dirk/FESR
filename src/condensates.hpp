#ifndef SRC_CONDENSATES_HPP
#define SRC_CONDENSATES_HPP

#include <math.h>
#include <vector>
#include <iostream>
#include <iomanip>

class Condensates {
 public:
  Condensates(
    const double &asTau = 0.3156,
    const std::vector<double> &qqMTau = {
      -pow(0.272, 3), -pow(0.272, 3), 0.8*(-pow(0.272, 3))
    },
    const std::vector<double> &mq = {
      2.8e-3, 5.0e-3, 97.0e-3
    }
  ) {

    qqMTau_ = qqMTau;
    const double aTau = asTau/M_PI;
    qqInv_ = {
      qqMTau_[0] + 3.0*pow(mq[0], 3)/(7.0*pow(M_PI, 2))*(1.0/aTau - 53.0/24.0),
      qqMTau_[1] + 3.0*pow(mq[1], 3)/(7.0*pow(M_PI, 2))*(1.0/aTau - 53.0/24.0),
      qqMTau_[2] + 3.0*pow(mq[2], 3)/(7.0*pow(M_PI, 2))*(1.0/aTau - 53.0/24.0)
    };
  };

  std::vector<double> qqMTau_; // ( <uu>(mTau), <dd>(mTau), <ss>(Mtau) )
  std::vector<double> qqInv_;
};

#endif
