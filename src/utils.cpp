#include "utils.hpp"

void addCSV(double value, ofstream &outputFile, bool last = false) {
  if(!last) {
    outputFile << value << ",";
  } else {
    outputFile << value;
  }
}

string variableOutput(Minimizer* min, const uint &i) {
  const bool isFixed = min->IsFixedVariable(i);
  const double *xs = min->X();
  const double *errors = min->Errors();
  string out = "";

  // name
  out += min->VariableName(i) + "\t=\t" + std::to_string(xs[0]);

  // error
  if(!isFixed) {
    out += "\t+/-\t" + std::to_string(errors[0]) + "\t";
  }

  // fixed
  if(isFixed) {
    out += "\t(fixed)";
  }

  out += "\n";
  return out;
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
    const string &configFilePath,
    Minimizer* min
) {
  std::size_t pos = configFilePath.find_last_of("/\\");
  pos = configFilePath.find_last_of(".");
  std::string outputFilePath = configFilePath.substr(0,pos) + ".txt";

  ofstream outputFile;
  outputFile.open(outputFilePath);
  outputFile << std::setprecision(15);

  ifstream configFile;
  configFile.open(configFilePath);

  if(outputFile.is_open() && configFile.is_open())
  {
    outputFile << "Minuit2Minimizer\t Status:\t" << min->Status() << endl;
    outputFile << "FVAL\t=\t" << min->MinValue() << endl;
    outputFile << "Edm\t=\t" << min->Edm() << endl;
    outputFile << "Nfcn\t=\t" << min->NCalls() << endl;

    const double *xs = min->X();
    const double *errors = min->Errors();

    // output variables
    for(uint i = 0; i <= 11; i++) {
      outputFile << variableOutput(min, i);
    }

    outputFile << endl;

    // add configFile
    outputFile << configFile.rdbuf();

    configFile.close();
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

  writeToCSV(csvFilePath, min, config);
  writeToOutputFile(configFilePath, min);
}


void writeOutput(const string text, const string outputFilePath)
{
  ofstream outputFile;
  outputFile.open(outputFilePath, std::ios::app);
  outputFile << text << endl;
}
