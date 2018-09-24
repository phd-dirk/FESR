#include "utils.hpp"

void addCSV(double value, ofstream &outputFile, bool last = false) {
  if(!last) {
    outputFile << value << ",";
  } else {
    outputFile << value;
  }
}

void writeToCSV(
    string &csvFilePath,
    Minimizer* min,
    const Configuration config
) {
  ofstream csvFile;
  csvFile.open(csvFilePath, std::ios::app);
  csvFile << std::setprecision(15);

  if(csvFile.is_open())
  {
    const double *xs = min->X();
    const double *errors = min->Errors();
    // min status
    addCSV(min->Status(), csvFile);
    // alpha_s
    addCSV(xs[0], csvFile);
    addCSV(errors[0], csvFile);
    // aGGInv
    addCSV(xs[1], csvFile);
    addCSV(errors[1], csvFile);
    // O6
    addCSV(xs[2], csvFile);
    addCSV(errors[2], csvFile);
    // O8
    addCSV(xs[2], csvFile);
    addCSV(errors[2], csvFile);

    // write FVAL
    addCSV(min->MinValue(), csvFile);
    // dof
    addCSV(config.dof(), csvFile);
    addCSV(min->MinValue()/config.dof(), csvFile);
    addCSV(min->Edm(), csvFile, true);

    csvFile << endl;
    csvFile.close();
  }
}

void writeToOutputFile(
    string &outputFilePath,
    Minimizer* min
) {
  ofstream outputFile;
  outputFile.open(outputFilePath);
  outputFile << std::setprecision(15);

  if(outputFile.is_open())
  {
    outputFile << "Minuit2Minimizer\t Status:\t" << min->Status() << endl;
    outputFile << "FVAL\t=\t" << min->MinValue() << endl;
    outputFile << "Edm\t=\t" << min->Edm() << endl;
    outputFile << "Nfcn\t=\t" << min->NCalls() << endl;

    const double *xs = min->X();
    const double *errors = min->Errors();
    outputFile << "astau\t=\t" << xs[0] << "\t+/-\t" << errors[0] << endl;
    outputFile << "aGGInv\t=\t" << xs[1] << "\t+/-\t" << errors[1] << endl;
    outputFile << "O6\t=\t" << xs[2] << "\t+/-\t" << errors[2] << endl;
    outputFile << "O8\t=\t" << xs[3] << "\t+/-\t" << errors[3] << endl;

    // DV
    outputFile << "deltaV\t=\t" << xs[4] << "\t+/-\t" << errors[4] << endl;
    outputFile << "gammaV\t=\t" << xs[5] << "\t+/-\t" << errors[5] << endl;
    outputFile << "alphaV\t=\t" << xs[6] << "\t+/-\t" << errors[6] << endl;
    outputFile << "betaV\t=\t" << xs[7] << "\t+/-\t" << errors[7] << endl;
    outputFile << "deltaA\t=\t" << xs[8] << "\t+/-\t" << errors[8] << endl;
    outputFile << "gammaA\t=\t" << xs[9] << "\t+/-\t" << errors[9] << endl;
    outputFile << "alphaA\t=\t" << xs[10] << "\t+/-\t" << errors[10] << endl;
    outputFile << "betaA\t=\t" << xs[11] << "\t+/-\t" << errors[11] << endl;

    outputFile.close();
  }
}

void writeOutput (
    const string configFilePath,
    Minimizer* min,
    const Configuration config
) {
  std::size_t pos = configFilePath.find_last_of("/\\");
  std::string csvFilePath = configFilePath.substr(0,pos) + "/fits.csv";
  pos = configFilePath.find_last_of(".");
  std::string outputFilePath = configFilePath.substr(0,pos) + ".txt";

  writeToCSV(csvFilePath, min, config);
  writeToOutputFile(outputFilePath, min);
}


void writeOutput(const string text, const string outputFilePath)
{
  ofstream outputFile;
  outputFile.open(outputFilePath, std::ios::app);
  outputFile << text << endl;
}
