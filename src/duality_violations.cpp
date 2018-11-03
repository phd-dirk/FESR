#import "duality_violations.hpp"

double DV::cIntVpA(
  const double &s0, const Weight &w,
  const double &deV, const double &gaV, const double &alV, const double &beV,
  const double &deA, const double &gaA, const double &alA, const double &beA
) {
  return (
    cIntVA(s0, w, deV, gaV, alV, beV)
    +cIntVA(s0, w, deA, gaA, alA, beA)
  );
}

double DV::cIntVA(
  const double &s0, const Weight &w, const double &de, const double &ga,
  const double &al, const double &be
) {
  return w.poli().x0*intP0(s0, de, ga, al, be)
    + w.poli().x1*intP1(s0, de, ga, al, be)
    + w.poli().x2*intP2(s0, de, ga, al, be)
    + w.poli().x3*intP3(s0, de, ga, al, be);
}

double DV::rhoDV(
  const double &s, const double &delta, const double &gamma,
  const double &alpha, const double &beta
) {
  return exp(-delta-gamma*s)*sin(alpha+beta*s);
}

double DV::intP0(
  const double &s0, const double &delta, const double &gamma,
  const double &alpha, const double &beta
) {
  return (
    -12.0*pow(M_PI, 2)*exp(-delta - gamma*s0)*
    (beta*cos(alpha + beta*s0) + gamma*sin(alpha + beta*s0)))/(pow(beta,2) + pow(gamma,2)
  )/s0;
}

double DV::intP1(
  const double &s0, const double &delta, const double &gamma,
  const double &alpha, const double &beta
) {
  return -12.0*pow(M_PI, 2)*(
    (exp(-delta - gamma*s0)*(beta*(2*gamma + (pow(beta,2) + pow(gamma,2))*s0)*cos(alpha + beta*s0) + (pow(gamma,2) + pow(gamma,3)*s0 + pow(beta,2)*(-1 + gamma*s0))*sin(alpha + beta*s0)))/
                             pow(pow(beta,2) + pow(gamma,2),2)
  )/pow(s0, 2);
}

double DV::intP2(
  const double &s0, const double &delta, const double &gamma,
  const double &alpha, const double &beta
) {
  return -12.0*pow(M_PI, 2)*(
    (exp(-delta - gamma*s0)*(beta*(-2*pow(beta,2) + 6*pow(gamma,2) + 4*gamma*(pow(beta,2) + pow(gamma,2))*s0 + pow(pow(beta,2) + pow(gamma,2),2)*pow(s0,2))*cos(alpha + beta*s0) +(2*gamma*(-3*pow(beta,2) + pow(gamma,2)) + 2*(-pow(beta,4) + pow(gamma,4))*s0 + gamma*pow(pow(beta,2) + pow(gamma,2),2)*pow(s0,2))*sin(alpha + beta*s0)))/
    pow(pow(beta,2) + pow(gamma,2),3)
  )/pow(s0, 3);
}

double DV::intP3(
  const double &s0, const double &delta, const double &gamma,
  const double &alpha, const double &beta
) {
  return -12.0*pow(M_PI, 2)*(exp(-delta - gamma*s0)*(
            beta*(pow(beta,6)*pow(s0,3) + 3*pow(beta,4)*s0*(-2 + 2*gamma*s0 + pow(gamma,2)*pow(s0,2)) +
                  3*pow(beta,2)*gamma*(-8 + 4*gamma*s0 + 4*pow(gamma,2)*pow(s0,2) + pow(gamma,3)*pow(s0,3)) + 
                  pow(gamma,3)*(24 + 18*gamma*s0 + 6*pow(gamma,2)*pow(s0,2) + pow(gamma,3)*pow(s0,3)))*cos(alpha + beta*s0) + 
            (pow(beta,6)*pow(s0,2)*(-3 + gamma*s0) + 3*pow(beta,4)*(2 - 6*gamma*s0 - pow(gamma,2)*pow(s0,2) + pow(gamma,3)*pow(s0,3)) + 
             3*pow(beta,2)*pow(gamma,2)*(-12 - 4*gamma*s0 + pow(gamma,2)*pow(s0,2) + pow(gamma,3)*pow(s0,3)) + 
             pow(gamma,4)*(6 + 6*gamma*s0 + 3*pow(gamma,2)*pow(s0,2) + pow(gamma,3)*pow(s0,3)))*sin(alpha + beta*s0)))/pow(pow(beta,2) + pow(gamma,2),4
             )/pow(s0, 4);
}
