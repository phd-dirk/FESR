#ifndef SRC_CONDENSATES_HPP
#define SRC_CONDENSATES_HPP

#include <math.h>
#include <vector>

class Condensates {
 public:
  Condensates(
    const double &asTau = 0.3156,
    const std::vector<double> &qqMTau = {
      0.0201236, 0.0201236, 0.0160989
    },
    const std::vector<double> &mq = {
      2.8e-3, 5.0e-3, 97.0e-3
    }
  ) {
    qqMTau_ = qqMTau;
    qqInv_ = {
      qqMTau_[0] + 3.*pow(mq[1], 3)/(7.*pow(M_PI, 2))*(1./asTau/M_PI - 53./24.),
      qqMTau_[1] + 3.*pow(mq[2], 3)/(7.*pow(M_PI, 2))*(1./asTau/M_PI - 53./24.),
      qqMTau_[2] + 3.*pow(mq[3], 3)/(7.*pow(M_PI, 2))*(1./asTau/M_PI - 53./24.)
    };
  };

  std::vector<double> qqMTau_; // ( <uu>(mTau), <dd>(mTau), <ss>(Mtau) )
  std::vector<double> qqInv_;
};

#endif
