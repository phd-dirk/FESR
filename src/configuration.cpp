#include "./configuration.hpp"

Configuration::Configuration(string configFilePath) {
    ifstream configFile(configFilePath);
    json jsonConfig;
    configFile >> jsonConfig;

    order = jsonConfig["parameters"]["order"];
    RVANormalization = jsonConfig["parameters"]["RVANormalization"];

    nc = jsonConfig["parameters"]["nc"];
    nf = jsonConfig["parameters"]["nf"];

    thMomContribs = {
      jsonConfig["scheme"], jsonConfig["thMomContribs"]["D0"], jsonConfig["thMomContribs"]["D4"],
      jsonConfig["thMomContribs"]["D68"], jsonConfig["thMomContribs"]["DV"], jsonConfig["thMomContribs"]["PionPole"]
    };

    astau = {
      jsonConfig["variables"]["astau"]["fixed"],
      jsonConfig["variables"]["astau"]["value"],
      jsonConfig["variables"]["astau"]["stepSize"]
    };
    aGGInv = {
      jsonConfig["variables"]["aGGInv"]["fixed"],
      jsonConfig["variables"]["aGGInv"]["value"],
      jsonConfig["variables"]["aGGInv"]["stepSize"]
    };
    rhoVpA = {
      jsonConfig["variables"]["rhoVpA"]["fixed"],
      jsonConfig["variables"]["rhoVpA"]["value"],
      jsonConfig["variables"]["rhoVpA"]["stepSize"]
    };
    c8VpA = {
      jsonConfig["variables"]["c8VpA"]["fixed"],
      jsonConfig["variables"]["c8VpA"]["value"],
      jsonConfig["variables"]["c8VpA"]["stepSize"]
    };
    deltaV = {
      jsonConfig["variables"]["deltaV"]["fixed"],
      jsonConfig["variables"]["deltaV"]["value"],
      jsonConfig["variables"]["deltaV"]["stepSize"]
    };
    gammaV = {
      jsonConfig["variables"]["gammaV"]["fixed"],
      jsonConfig["variables"]["gammaV"]["value"],
      jsonConfig["variables"]["gammaV"]["stepSize"]
    };
    alphaV = {
      jsonConfig["variables"]["alphaV"]["fixed"],
      jsonConfig["variables"]["alphaV"]["value"],
      jsonConfig["variables"]["alphaV"]["stepSize"]
    };
    betaV = {
      jsonConfig["variables"]["betaV"]["fixed"],
      jsonConfig["variables"]["betaV"]["value"],
      jsonConfig["variables"]["betaV"]["stepSize"]
    };
    deltaA = {
      jsonConfig["variables"]["deltaA"]["fixed"],
      jsonConfig["variables"]["deltaA"]["value"],
      jsonConfig["variables"]["deltaA"]["stepSize"]
    };
    gammaA = {
      jsonConfig["variables"]["gammaA"]["fixed"],
      jsonConfig["variables"]["gammaA"]["value"],
      jsonConfig["variables"]["gammaA"]["stepSize"]
    };
    alphaA = {
      jsonConfig["variables"]["alphaA"]["fixed"],
      jsonConfig["variables"]["alphaA"]["value"],
      jsonConfig["variables"]["alphaA"]["stepSize"]
    };
    betaA = {
      jsonConfig["variables"]["betaA"]["fixed"],
      jsonConfig["variables"]["betaA"]["value"],
      jsonConfig["variables"]["betaA"]["stepSize"]
    };


    mTau = jsonConfig["parameters"]["mTau"];
    sTau = pow(mTau, 2);

    be = jsonConfig["parameters"]["be"];
    dBe = jsonConfig["parameters"]["dBe"];

    // add weights & s0s
    for(auto const& input : jsonConfig["parameters"]["input"]) {
      Weight w(input["weight"].get<int>());
      vec s0s = input["s0s"].get<vec>();
      inputs.push_back({ w, s0s });
      momCount += s0s.size();
    }
    initializeBetaCoefficients();
    initializeAdlerCoefficients();

    // minuit
    tolerance = jsonConfig["tolerance"];
  }
