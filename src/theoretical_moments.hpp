#ifndef SRC_THEORETICAL_MOMENTS_H
#define SRC_THEORETICAL_MOMENTS_H

#include "./types.hpp"
#include "./configuration.hpp"
#include "./adler_function.hpp"
#include "./duality_violations.hpp"

template <typename T>
void print(T t)
{
  std::cout << t << " ";
}

template<typename T, typename... Args>
void print(T t, Args... args)
{
  std::cout << t << " ";
  print(args...);
  std::cout << std::endl;
}

class TheoreticalMoments:
    public AdlerFunction, DualityViolations
{
 public:
  TheoreticalMoments(const Configuration &config);

  double thMom(
    cDbl &s0, const Weight &w, cDbl &astau,
    cDbl &aGGinv, cDbl &rhoVpA, cDbl &c8VpA, cDbl &order,
    cDbl &deV, cDbl &gaV, cDbl &alV, cDbl &beV,
    cDbl &deA, cDbl &gaA, cDbl &alA, cDbl &beA
  ) const;

  double cIntVpAD0FO(const double &s0, const Weight &weight, const double &sTau,
                     const double &astau, const int &order) const {
    return 2.*D0CIntFO(s0, weight, sTau, astau, order);
  }
  double cIntVpAD0CI(const double &s0, const Weight &weight, const double &sTau,
                     const double &astau, const int &order) const {
    return 2.*D0CIntCI(s0, weight, sTau, astau, order);
  }


  double cIntVpAD4FO(const double &s0, const Weight &weight, const double &sTau,
                     const double &astau, const double &aGGinv) const {
    return  D4CInt(s0, weight, sTau, astau, aGGinv, 1) + D4CInt(s0, weight, sTau, astau, aGGinv, -1);
  }
  double del0(const double &s0, const Weight &weight,
              const double &sTau, const double &astau, const int &order) const {
    return (cIntVpAD0FO(s0, weight, sTau, astau, order)
            - cIntVpAD0FO(s0, weight, sTau, astau, 0)
            )/3.0;
  }

  double del2(const double &s0, const Weight &weight,
              const double &astau, const int &order) const {
    return (D2CInt(s0, weight, astau, order, 1)
            + D2CInt(s0, weight, astau, order, -1)
            )/3.;
  }

  double del4(const double &s0, const Weight &weight, const double &sTau,
              const double &astau, const double &aGGinv) const {
    return cIntVpAD4FO(s0, weight, sTau, astau, aGGinv)/3.;
  }

  double del6(const double &s0, const Weight &weight,
               const double &rhoVpA) const {
    return D68CInt(s0, weight, rhoVpA, 0.0)/3.;
  }

  double del8(const double &s0, const Weight &weight,
              const double &c8VpA) const {
    return D68CInt(s0, weight, 0.0, c8VpA)/3.;
  }

  double del68(const double &s0, const Weight &weight,
               const double &rhoVpA, const double &c8VpA) {
    return D68CInt(s0, weight, rhoVpA, c8VpA)/3.;
  }

  void log(const double &astau, const double &aGGinv, const double &rhoVpA, const double &c8VpA, const int &order) const {
    // cout << "thMom: \t" << operator() (config_.sTau, astau, aGGinv, rhoVpA, c8VpA, order) << endl;
    // cout << "Delta^(0): \t" << del0(config_.sTau, config_.weight, config_.sTau, astau, order) << endl;
    // cout << "Delta^(4): \t" << del4(config_.sTau, config_.weight, config_.sTau, astau, aGGinv) << endl;
    // cout << "Delta^(6): \t" << del6(config_.sTau, config_.weight, rhoVpA) << endl;
    // cout << "Delta^(8): \t" << del8(config_.sTau, config_.weight, c8VpA) << endl;
    // cout << "Delta_P+S: \t" << deltaP(config_.sTau, config_.weight) << endl;
  }


 private:
  Configuration config_;
  std::vector<Input> inputs_;
};

#endif
