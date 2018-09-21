#include "utils.hpp"

void addCSV(double value, ofstream &outputFile, bool last = false) {
  if(!last) {
    outputFile << value << ",";
  } else {
    outputFile << value;
  }
}

void writeOutput (
    const string outputFilePath,
    Minimizer* min,
    const Configuration config
 )
{
  ofstream outputFile;
  outputFile.open(outputFilePath, std::ios::app);
  outputFile << std::setprecision(15);

  if(outputFile.is_open())
  {
    const double *xs = min->X();
    const double *errors = min->Errors();
    // alpha_s
    addCSV(xs[0], outputFile);
    addCSV(errors[0], outputFile);
    // aGGInv
    addCSV(xs[1], outputFile);
    addCSV(errors[1], outputFile);
    // O6
    addCSV(xs[2], outputFile);
    addCSV(errors[2], outputFile);
    // O8
    addCSV(xs[2], outputFile);
    addCSV(errors[2], outputFile);

    // write FVAL
    addCSV(min->MinValue(), outputFile);
    // dof
    addCSV(config.dof(), outputFile);
    addCSV(min->MinValue()/config.dof(), outputFile);
    addCSV(min->Edm()/config.dof(), outputFile, true);

    outputFile << endl;
    outputFile.close();
  }
}

void writeOutput(const string text, const string outputFilePath)
{
  ofstream outputFile;
  outputFile.open(outputFilePath, std::ios::app);
  outputFile << text << endl;
}
