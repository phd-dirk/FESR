// MINUIT2
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "./chisquared.hpp"

class Minuit {
 public:
  static ROOT::Math::Minimizer* FESR(const Configuration &config) {
    const Chi2 chi2(config);

    // MINUIT
    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    // set tolerances
    // min->SetMaxFunctionCalls(10000000); // for Minuit2
    // min->SetMaxIteraions(10000000); // for GSL
    // min->SetTolerance(1e-13);
    min->SetStrategy(2);
    min->SetPrintLevel(3); // activate logging

    // function wrapper
    ROOT::Math::Functor chi2F(chi2, 12);
    min->SetFunction(chi2F);

    // set free variables to be minimized
    if (config.astau.isFixed) {
      min->SetFixedVariable(0, "astau", config.astau.value);
    } else {
      min->SetVariable(0, "astau", config.astau.value, config.astau.stepSize);
    }
    if (config.aGGInv.isFixed) {
      min->SetFixedVariable(1, "aGGInv", config.aGGInv.value);
    } else {
      min->SetVariable(1, "aGGInv", config.aGGInv.value, config.aGGInv.stepSize);
    }
    if (config.rhoVpA.isFixed) {
      min->SetFixedVariable(2, "rhoVpA", config.rhoVpA.value);
    } else {
      min->SetVariable(2, "rhoVpA", config.rhoVpA.value, config.rhoVpA.stepSize);
    }
    if (config.c8VpA.isFixed) {
      min->SetFixedVariable(3, "c8VpA", config.c8VpA.value);
    } else {
      min->SetVariable(3, "c8VpA", config.c8VpA.value, config.c8VpA.stepSize);
    }
    // D10
    if (config.c10.isFixed) {
      min->SetFixedVariable(4, "c10", config.c10.value);
    } else {
      min->SetVariable(4, "c10", config.c10.value, config.c10.stepSize);
    }
    // D12
    if (config.c12.isFixed) {
      min->SetFixedVariable(5, "c12", config.c12.value);
    } else {
      min->SetVariable(5, "c12", config.c12.value, config.c12.stepSize);
    }
    if (config.deltaV.isFixed) {
      min->SetFixedVariable(6, "deltaV", config.deltaV.value);
    } else {
      min->SetVariable(6, "deltaV", config.deltaV.value, config.deltaV.stepSize);
    }
    if (config.gammaV.isFixed) {
      min->SetFixedVariable(7, "gammaV", config.gammaV.value);
    } else {
      min->SetVariable(7, "gammaV", config.gammaV.value, config.gammaV.stepSize);
    }
    if (config.alphaV.isFixed) {
      min->SetFixedVariable(8, "alphaV", config.alphaV.value);
    } else {
      min->SetVariable(8, "alphaV", config.alphaV.value, config.alphaV.stepSize);
    }
    if (config.betaV.isFixed) {
      min->SetFixedVariable(9, "betaV", config.betaV.value);
    } else {
      min->SetVariable(9, "betaV", config.betaV.value, config.betaV.stepSize);
    }
    if (config.deltaA.isFixed) {
      min->SetFixedVariable(10, "deltaA", config.deltaA.value);
    } else {
      min->SetVariable(10, "deltaA", config.deltaA.value, config.deltaA.stepSize);
    }
    if (config.gammaA.isFixed) {
      min->SetFixedVariable(11, "gammaA", config.gammaA.value);
    } else {
      min->SetVariable(11, "gammaA", config.gammaA.value, config.gammaA.stepSize);
    }
    if (config.alphaA.isFixed) {
      min->SetFixedVariable(12, "alphaA", config.alphaA.value);
    } else {
      min->SetVariable(12, "alphaA", config.alphaA.value, config.alphaA.stepSize);
    }
    if (config.betaA.isFixed) {
      min->SetFixedVariable(13, "betaA", config.betaA.value);
    } else {
      min->SetVariable(13, "betaA", config.betaA.value, config.betaA.stepSize);
    }

    // minimize!
    min->Minimize();

    const double *xs = min->X();
    std::cout << "alpha: \t" << xs[0] << std::endl;
    std::cout << "aGGInv: \t" << xs[1] << std::endl;
    std::cout << "c6: \t" << xs[2] << std::endl;
    std::cout << "c8: \t" << xs[3] << std::endl;
    std::cout << "c10: \t" << xs[4] << std::endl;
    std::cout << "c12: \t" << xs[5] << std::endl;

    ThMoms::logDeltas(
      config.inputs_[0].s0s[0],
      config.inputs_[0].weight,
      config.sTau_,
      xs[0],
      xs[1],
      Configuration::adlerCoefficients(
        config.nf_, Configuration::betaCoefficients( config.nc_, config.nf_ )
      ),
      config.mq_,
      config.condensates_,
      config.order_,
      xs[2], xs[3], xs[4], xs[5]
    );

    return min;
  }
};
