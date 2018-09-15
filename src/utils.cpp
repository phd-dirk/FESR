#include "utils.hpp"

void writeOutput (
    const string configFilePath,
    const string outputFilePath,
    Minimizer* min,
    const Configuration config
 )
{
  ifstream configFile;
  configFile.open(configFilePath);
  ofstream outputFile;
  outputFile.open(outputFilePath, std::ios::out);
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

      // write FVAL & Edm
      outputFile << endl;
      outputFile << "FVAL: \t = " << min->MinValue() << endl;
      outputFile << "Edm: \t = " << min->Edm() << endl;

      // write vars with errors
      const double *xs = min->X();
      const double *errors = min->Errors();
      outputFile << endl;
      outputFile << "astau: \t = " << xs[0] << "\t +/- \t " << errors[0] << endl;
      outputFile << "aGG: \t = " << xs[1] << "\t +/- \t " << errors[1] << endl;
      outputFile << "rho: \t = " << xs[2] << "\t +/- \t " << errors[2] << endl;
      outputFile << "c8: \t = " << xs[3] << "\t +/- \t " << errors[3] << endl;

      // dof
      outputFile << "dof: \t" << config.dof() << endl;
      outputFile << "Chi2/dof = \t" << min->MinValue()/config.dof() << endl;

      outputFile << endl;
      outputFile.close();
    }
    configFile.close();
  }
}
