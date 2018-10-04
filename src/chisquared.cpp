#include "./chisquared.hpp"

// Public
Chisquared::Chisquared(Configuration config)
  :config_(config), inputs_(config.inputs),
   expMom_(ExperimentalMoments("/Users/knowledge/Developer/PhD/FESR/aleph.json", config))
{
  invCovMat_(config.momCount_, config.momCount_);
  std::cout << "mhmm" << std::endl;
  initInvCovMat();
}

double Chisquared::operator() ( const double *xx) const {
  // init fit parameters
  double astau = xx[0];
  double aGGinv = xx[1];
  double rhoD6VpA = xx[2];
  double c8D8VpA = xx[3];
  double deV = xx[4];
  double gaV = xx[5];
  double alV = xx[6];
  double beV = xx[7];
  double deA = xx[8];
  double gaA = xx[9];
  double alA = xx[10];
  double beA = xx[11];

  return chi2(
    config_.inputs, astau, aGGinv, rhoD6VpA, c8D8VpA,
    deV, gaV, alV, beV,
    deA, gaA, alA, alA
  );
}

double Chisquared::operator ()(
  std::vector<Input> inputs,
  const double &astau,
  const double &aGGinv,
  const double &rho,
  const double &c8,
  const double &vKappa,
  const double &vGamma,
  const double &vAlpha,
  const double &vBeta,
  const double &aKappa,
  const double &aGamma,
  const double &aAlpha,
  const double &aBeta
) const
{
  return chi2(
    inputs, astau, aGGinv, rho, c8,
    vKappa, vGamma, vAlpha, vBeta, aKappa, aGamma, aAlpha, aBeta
  );
}

// Private
double Chisquared::chi2(
  std::vector<Input> inputs,
  const double &astau,
  const double &aGGinv,
  const double &rho,
  const double &c8,
  const double &deV,
  const double &gaV,
  const double &alV,
  const double &beV,
  const double &deA,
  const double &gaA,
  const double &alA,
  const double &beA
) const
{
  double chi = 0;

  vec momDiff(config_.momCount_);
  vec thMoms = calcThMoms(
    inputs, astau, aGGinv, rho, c8, config_.order_,
    deV, gaV, alV, beV, deA, gaA, alA, beA
  );

  for(uint i = 0; i < config_.momCount_; i++) {
    momDiff[i] = expMom_()[i] - thMoms[i];
  }

  for(uint k = 0; k < config_.momCount_; k++) {
    for(uint l = 0; l < config_.momCount_; l++) {
      chi += momDiff[k] * invCovMat_(k, l) * momDiff[l];
    }
  }
  return chi;
}

vec Chisquared::calcThMoms(
  std::vector<Input> inputs,
  const double &astau,
  const double &aGGinv,
  const double &rhoVpA,
  const double &c8VpA,
  const double &order,
  const double &deV,
  const double &gaV,
  const double &alV,
  const double &beV,
  const double &deA,
  const double &gaA,
  const double &alA,
  const double &beA
) const
{
  vec thMoms(config_.momCount_);
  std::vector<std::future<double>> ftrs(config_.momCount_);
  TheoreticalMoments th(config_);
  int i = 0;
  for(auto const& input: inputs_) {
    vec s0s = input.s0s;
    Weight w = input.weight;
    for(auto const& s0: s0s) {
      ftrs[i] = std::async(
        &TheoreticalMoments::thMom,
        &th, s0, w, astau, aGGinv, rhoVpA, c8VpA, order,
        deV, gaV, alV, beV, deA, gaA, alA, beA
      );
      i++;
    }
  }

  for(uint i=0; i<config_.momCount_; i++) {
    thMoms[i] = ftrs[i].get();
  }

  return thMoms;
}

void Chisquared::initInvCovMat() {
  mat covMat = expMom_.getCovMat();

  std::cout << "jo" << std::endl;
  std::cout << config_.momCount_ << endl;
  // Remove correlations with R_tau,V+A in Aleph fit
  for (uint i = 1; i < config_.momCount_; i++) {
    covMat(0, i) = 0.;
    covMat(i, 0) = 0.;
  }
  // employ uncertainity of R_VA = 3.4718(72) (HFLAV 2017)
  covMat(0, 0) = pow(0.0072, 2);

  Numerics::invertMatrix(covMat, invCovMat_);
}
