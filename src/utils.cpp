#include "utils.hpp"

void writeOutput(const string configFilePath, const string outputFilePath, Minimizer* min) {
  ifstream configFile;
  configFile.open(configFilePath);
  ofstream outputFile;
  outputFile.open(outputFilePath, std::ios::app);
  outputFile << std::setprecision(15);

  string line;
  if (configFile.is_open())
  {
    if(outputFile.is_open())
    {
      // write config.json to output file
      while( getline (configFile, line) )
      {
        outputFile << line << "\n";
      }

      const double *xs = min->X();
      const double *errors = min->Errors();

      outputFile << endl;
      outputFile << "astau: \t = " << xs[0] << "\t +/- \t " << errors[0] << "\n";
      outputFile << "aGG: \t = " << xs[1] << "\t +/- \t " << errors[1] << "\n";
      outputFile << "rho: \t = " << xs[2] << "\t +/- \t " << errors[2] << "\n";
      outputFile << "c8: \t = " << xs[3] << "\t +/- \t " << errors[3] << "\n";

      outputFile << endl;
      outputFile.close();
    }
    configFile.close();
  }
}

// int dof(const json config) {
//   int numS0s = config["parameters"]["s0Set"].size();
//   int numVar = 0;
//   if (!config["variables"]["astau"]["fixed"]) {
//     numVar++;
//   }
//   if (!config["variables"]["aGGInv"]["fixed"]) {
//     numVar++;
//   }
//   if (!config["variables"]["rhoVpA"]["fixed"]) {
//     numVar++;
//   }
//   if (!config["variables"]["c8VpA"]["fixed"]) {
//     numVar++;
//   }
//   return numS0s-numVar;
// }

int dof(const Configuration config) {
  int dof = 1;
  // int dof = config.s0Set.size();
  // if(!config.astau.isFixed)
  //   dof--;
  // if(!config.aGGInv.isFixed)
  //   dof--;
  // if(!config.rhoVpA.isFixed)
  //   dof--;
  // if(!config.c8VpA.isFixed)
  //   dof--;
  return dof;
}
